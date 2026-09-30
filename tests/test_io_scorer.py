from __future__ import annotations

from dataclasses import replace
import os
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]

from vec_scorelens.io import InputPolicyError, load_bundle
from vec_scorelens.scorer import PINNED_COMMIT, PINNED_VERSION, run_official, verify_source_pin


def _write(path: Path, X: np.ndarray, genes: list[str], labels: list[str], *, sparse_x: bool = False) -> None:
    matrix = sparse.csr_matrix(X) if sparse_x else X
    data = ad.AnnData(matrix)
    data.var_names = genes
    data.obs_names = [f"cell-{index}" for index in range(X.shape[0])]
    data.obs["celltype"] = labels
    data.write_h5ad(path)


def _inputs(tmp_path: Path, *, sparse_x: bool = False) -> tuple[Path, Path, Path]:
    rng = np.random.default_rng(22)
    genes = ["g1", "g2", "g3", "g4"]
    labels = ["a"] * 16 + ["b"] * 16
    reference = rng.uniform(0.1, 1.0, size=(32, 4)).astype(np.float32)
    target = (reference + np.array([0.35, 0.0, 0.25, 0.0], dtype=np.float32)).astype(np.float32)
    prediction = (reference + np.array([0.25, 0.05, 0.18, 0.0], dtype=np.float32)).astype(np.float32)
    pred_path, target_path, ref_path = tmp_path / "prediction.h5ad", tmp_path / "target.h5ad", tmp_path / "reference.h5ad"
    _write(pred_path, prediction, genes, labels, sparse_x=sparse_x)
    _write(target_path, target, genes, labels, sparse_x=sparse_x)
    _write(ref_path, reference, genes, labels, sparse_x=sparse_x)
    return pred_path, target_path, ref_path


def test_load_bundle_is_deterministic_readonly_and_backed_sparse(tmp_path: Path) -> None:
    pred, target, reference = _inputs(tmp_path, sparse_x=True)
    one = load_bundle(pred, target, reference, max_cells=32, seed=9)
    two = load_bundle(pred, target, reference, max_cells=32, seed=9)
    assert one.prediction.dtype == np.float32
    assert one.prediction.flags.writeable is False
    assert one.target.flags.writeable is False
    assert one.reference.flags.writeable is False
    assert one.target_labels.flags.writeable is False
    assert one.provenance["selection"] == two.provenance["selection"]
    assert one.provenance["scope"] == "sampled_subset_only"
    assert one.provenance["inputs"]["prediction"]["original_shape"] == [32, 4]
    with pytest.raises(ValueError):
        one.prediction[0, 0] = 1.0


def test_load_bundle_rejects_gene_order_mismatch(tmp_path: Path) -> None:
    pred, target, reference = _inputs(tmp_path)
    data = ad.read_h5ad(reference)
    data = data[:, ["g2", "g1", "g3", "g4"]].copy()
    data.write_h5ad(reference)
    with pytest.raises(InputPolicyError, match="exactly match"):
        load_bundle(pred, target, reference, max_cells=32)


def test_load_bundle_rejects_selected_nonfinite_values(tmp_path: Path) -> None:
    pred, target, reference = _inputs(tmp_path)
    data = ad.read_h5ad(pred)
    data.X[0, 0] = np.nan
    data.write_h5ad(pred)
    with pytest.raises(InputPolicyError, match="NaN or infinite"):
        load_bundle(pred, target, reference, max_cells=32)


def test_load_bundle_rejects_tiny_negative_before_float32_rounding(tmp_path: Path) -> None:
    pred, target, reference = _inputs(tmp_path)
    data = ad.read_h5ad(pred)
    data.X = np.asarray(data.X, dtype=np.float64)
    data.X[0, 0] = -1e-50
    data.write_h5ad(pred)
    with pytest.raises(InputPolicyError, match="negative"):
        load_bundle(pred, target, reference, max_cells=32)


def test_load_bundle_rejects_missing_expression_matrix(tmp_path: Path) -> None:
    pred, target, reference = _inputs(tmp_path)
    data = ad.AnnData(X=None, shape=(32, 4))
    data.var_names = ["g1", "g2", "g3", "g4"]
    data.obs_names = [f"cell-{index}" for index in range(32)]
    data.obs["celltype"] = ["a"] * 16 + ["b"] * 16
    data.write_h5ad(pred)
    with pytest.raises(InputPolicyError, match="X is missing"):
        load_bundle(pred, target, reference, max_cells=32)


def test_constant_target_fails_before_pinned_mmd_for_loader_and_direct_bundle(tmp_path: Path) -> None:
    pred, target, reference = _inputs(tmp_path)
    target_data = ad.read_h5ad(target)
    target_data.X[:] = 1.0
    target_data.write_h5ad(target)
    with pytest.raises(InputPolicyError, match="target-fitted MMD requires within-target variation"):
        load_bundle(pred, target, reference, max_cells=32)

    valid_bundle = load_bundle(pred, reference, reference, max_cells=32)
    constant = np.ones_like(valid_bundle.target)
    constant.setflags(write=False)
    with pytest.raises(InputPolicyError, match="target-fitted MMD requires within-target variation"):
        run_official(replace(valid_bundle, target=constant))


def test_cli_constant_target_exits_cleanly_without_partial_output(tmp_path: Path) -> None:
    pred, target, reference = _inputs(tmp_path)
    target_data = ad.read_h5ad(target)
    target_data.X[:] = 1.0
    target_data.write_h5ad(target)
    output = tmp_path / "constant-target-report"
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "vec_scorelens.cli",
            "--task", "T1",
            "--prediction", str(pred),
            "--target", str(target),
            "--reference", str(reference),
            "--out", str(output),
            "--max-cells", "32",
        ],
        capture_output=True,
        text=True,
        env=env,
        timeout=60,
    )
    assert result.returncode == 2
    assert "target-fitted MMD requires within-target variation" in result.stderr
    assert "Traceback" not in result.stderr
    assert not output.exists()


def test_run_official_matches_public_file_api_on_complete_selected_inputs(tmp_path: Path) -> None:
    pred, target, reference = _inputs(tmp_path)
    bundle = load_bundle(pred, target, reference, max_cells=32, seed=0)
    before = (bundle.prediction.copy(), bundle.target.copy(), bundle.reference.copy())
    scored = run_official(bundle, seed=3)
    import veckit

    public = veckit.score(task="T1", input=pred, target=target, reference=reference, seed=3)
    assert scored["metrics"] == public["metrics"]
    assert scored["provenance"]["scope"] == "sampled_subset_only"
    assert scored["provenance"]["commit"] == PINNED_COMMIT
    assert scored["provenance"]["version"] == PINNED_VERSION
    assert np.array_equal(bundle.prediction, before[0])
    assert np.array_equal(bundle.target, before[1])
    assert np.array_equal(bundle.reference, before[2])


def test_source_pin_is_verified() -> None:
    provenance = verify_source_pin()
    assert provenance["commit"] == PINNED_COMMIT
    assert provenance["version"] == PINNED_VERSION
    assert set(provenance["source_hashes"]) >= {"score_h5ad.py", "T1/metrics_v2.py"}
