"""Independent acceptance through public APIs; no implementation inspection.

The PI scope amendment in PLANTED_FAILURE_PROTOCOL.md distinguishes exposed
evidence from unsupported state localization, module discovery, and severity
labels. VEC_SCORELENS_SYNTHETIC_EVIDENCE optionally saves measured JSONL evidence
to a new file; it never overwrites an existing evidence file.
"""

from functools import lru_cache
import json
import os
from pathlib import Path

import numpy as np
import pytest

from vec_scorelens.diagnostics import diagnose
from vec_scorelens.io import Bundle, load_bundle
from vec_scorelens.scorer import run_official
from vec_scorelens.synthetic import create_case, create_healthy_sample, write_case


FAILURES = ("direction", "covariance", "missing", "overrepresented", "within_state",
            "scale", "magnitude", "response_magnitude", "collapse", "combined")
_EVIDENCE = None


@pytest.fixture(scope="session", autouse=True)
def evidence_file():
    global _EVIDENCE
    destination = os.environ.get("VEC_SCORELENS_SYNTHETIC_EVIDENCE")
    if destination:
        path = Path(destination)
        path.parent.mkdir(parents=True, exist_ok=True)
        _EVIDENCE = path.open("x", encoding="utf-8")
    yield
    if _EVIDENCE is not None:
        _EVIDENCE.close()
        _EVIDENCE = None


def _bundle(case):
    arrays = []
    for name in ("prediction", "target", "reference"):
        matrix = np.asarray(case[name].X, dtype=np.float32)
        matrix.flags.writeable = False
        arrays.append(matrix)
    return Bundle(prediction=arrays[0], target=arrays[1], reference=arrays[2],
                  genes=case["target"].var_names.tolist(),
                  target_labels=case["target"].obs["celltype"].astype(str).to_numpy(),
                  provenance={"synthetic_only": True, "assessed_scope": "full synthetic fixture",
                              "source_shapes": case["planted"]["shapes"]})


def _module_errors(diagnostic):
    """Test-only prespecified module summaries, not inferred product modules."""
    errors = [[] for _ in range(4)]
    variances = [[] for _ in range(4)]
    for row in diagnostic["genes"]:
        index = int(row["gene"].rsplit("_", 1)[1])
        if 96 <= index < 128:
            variances[(index - 96) // 8].append(row["variance_truth"])
    for row in diagnostic["pairs"]:
        i = int(row["gene_i"].rsplit("_", 1)[1])
        j = int(row["gene_j"].rsplit("_", 1)[1])
        if 96 <= i < 128 and 96 <= j < 128 and (i - 96) // 8 == (j - 96) // 8:
            errors[(i - 96) // 8].append(abs(row["covariance_error"]))
    return [float(np.mean(error) / max(np.mean(variance), 1e-12))
            for error, variance in zip(errors, variances)]


def _record(case, scored, diagnostic):
    if _EVIDENCE is None:
        return
    rows = sorted(diagnostic["genes"], key=lambda row: (-row["absolute_error"], row["gene"]))
    record = {"case": case["planted"]["case"], "seed": case["planted"]["seed"],
              "replicate": case["planted"].get("replicate"),
              "synthetic_only": True, "matrix_sha256": case["planted"]["matrix_sha256"],
              "official_metrics": scored["metrics"],
              "composition": diagnostic["composition"], "top_gene_errors": rows[:12],
              "max_gene_absolute_error": max(row["absolute_error"] for row in rows),
              "test_only_prespecified_module_normalized_covariance_error": _module_errors(diagnostic)}
    _EVIDENCE.write(json.dumps(record, allow_nan=False) + "\n")
    _EVIDENCE.flush()


@lru_cache(maxsize=36)
def _evaluate(name, seed):
    case = create_case(name, seed)
    bundle = _bundle(case)
    scored = run_official(bundle, seed=seed)
    diagnostic = diagnose(bundle, scored, seed=seed)
    _record(case, scored, diagnostic)
    return case, scored, diagnostic


def _gene_recovery(diagnostic, planted, count=8):
    rows = sorted(diagnostic["genes"], key=lambda row: (-row["absolute_error"], row["gene"]))
    selected = {row["gene"] for row in rows[:count]}
    hit = len(selected.intersection(planted))
    assert hit / len(planted) >= 0.875
    assert 1.0 - len(selected.difference(planted)) / (128 - len(planted)) >= 0.99


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("name", FAILURES)
def test_ten_planted_failures_recover_specific_measured_evidence(name, seed):
    case, scored, diagnostic = _evaluate(name, seed)
    genes = {row["gene"]: row for row in diagnostic["genes"]}
    composition = {row["celltype"]: row for row in diagnostic["composition"]}
    metrics = scored["metrics"]
    selected = [f"synthetic_gene_{i:03d}" for i in [64, 65, 66, 67, 80, 81, 82, 83]]
    if name in ("direction", "combined"):
        assert all(not genes[g]["literal_sign_agreement"] for g in selected)
        response_rows = {g: row for g, row in genes.items() if 64 <= int(g.rsplit("_", 1)[1]) < 96}
        top = sorted(response_rows, key=lambda g: (-response_rows[g]["absolute_error"], g))[:8]
        assert len(set(top).intersection(selected)) / 8 >= 0.875
        assert metrics["de_score"] < 0.9
    if name in ("covariance", "combined"):
        assert sum(value > 0.50 for value in _module_errors(diagnostic)) >= 3
        if name == "covariance":
            assert max(row["absolute_error"] for row in genes.values()) < 1e-4
            assert metrics["variance_ratio"] == pytest.approx(1.0, abs=0.002)
            assert metrics["de_score"] == pytest.approx(1.0, abs=0.001)
    if name in ("missing", "combined"):
        row = composition["state_3"]
        assert row["pred_soft_fraction"] <= 0.02
        assert abs(row["proportion_error"]) >= 0.075
        assert max(composition, key=lambda g: abs(composition[g]["proportion_error"])) == "state_3"
    if name == "overrepresented":
        row = composition["state_3"]
        assert row["pred_soft_fraction"] >= 0.27
        assert row["proportion_error"] >= 0.15
        assert max(composition, key=lambda g: abs(composition[g]["proportion_error"])) == "state_3"
    if name == "within_state":
        _gene_recovery(diagnostic, set(case["planted"]["planted_genes"]))
        assert all(abs(row["proportion_error"]) <= 0.02 for row in composition.values())
        # Only global gene-block evidence is exposed; subtype localization is not claimed.
        assert all(genes[g]["change_error"] == pytest.approx(0.15, abs=1e-4)
                   for g in case["planted"]["planted_genes"])
    if name == "scale":
        assert metrics["pb_rel_err"] >= 0.50
        assert metrics["pseudobulk_pearson"] == pytest.approx(1.0, abs=0.0001)
        assert all(row["mean_pred"] == pytest.approx(1.8 * row["mean_truth"] + 0.7, abs=1e-4)
                   for row in genes.values())
    if name == "magnitude":
        _gene_recovery(diagnostic, set(selected))
        assert all(genes[g]["change_pred"] / genes[g]["change_truth"] == pytest.approx(1.8, abs=0.02)
                   for g in selected)
        assert all(genes[g]["literal_sign_agreement"] for g in selected)
        for i in case["planted"]["response_gene_indices"]:
            g = f"synthetic_gene_{i:03d}"
            if g not in selected:
                assert genes[g]["change_pred"] / genes[g]["change_truth"] == pytest.approx(1.0, abs=0.02)
    if name == "response_magnitude":
        for i in case["planted"]["response_gene_indices"]:
            g = f"synthetic_gene_{i:03d}"
            assert genes[g]["change_pred"] / genes[g]["change_truth"] == pytest.approx(2.0, abs=0.02)
            assert genes[g]["literal_sign_agreement"]
        _, healthy_score, _ = _evaluate("healthy", seed)
        assert metrics["de_score"] == healthy_score["metrics"]["de_score"]
        assert metrics["de_direction"] == pytest.approx(healthy_score["metrics"]["de_direction"], abs=0.002)
    if name == "collapse":
        assert metrics["variance_ratio"] <= 0.01
        assert max(row["absolute_error"] for row in genes.values()) < 1e-4
        assert all(row["variance_ratio"] <= 0.01 for row in genes.values())
        assert metrics["de_score"] == pytest.approx(1.0, abs=0.001)
        assert sum(value > 0.50 for value in _module_errors(diagnostic)) == 4


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("name", ["healthy", "row_permutation"])
def test_exact_healthy_controls_preserve_population_evidence(name, seed):
    _, scored, diagnostic = _evaluate(name, seed)
    assert max(row["absolute_error"] for row in diagnostic["genes"]) < 1e-4
    assert max(_module_errors(diagnostic)) < 1e-4
    assert all(abs(row["proportion_error"]) < 0.02 for row in diagnostic["composition"])
    assert scored["metrics"]["variance_ratio"] == pytest.approx(1.0, abs=0.002)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_twenty_independent_healthy_controls_per_seed(seed):
    # Raw diagnostic bounds, not fabricated product severity/failure categories.
    for replicate in range(20):
        case = create_healthy_sample(seed, replicate)
        bundle = _bundle(case)
        scored = run_official(bundle, seed=seed)
        diagnostic = diagnose(bundle, scored, seed=seed)
        _record(case, scored, diagnostic)
        validation_labels = np.asarray(case["planted"]["privileged_test_only"]["prediction_state_ids"])
        assert np.mean(scored["probe"].predict(bundle.prediction) == validation_labels) >= 0.99
        assert max(row["absolute_error"] for row in diagnostic["genes"]) < 1e-4
        assert max(_module_errors(diagnostic)) < 0.50
        assert all(abs(row["proportion_error"]) < 0.02 for row in diagnostic["composition"])
        assert 0.8 < scored["metrics"]["variance_ratio"] < 1.2
        assert scored["metrics"]["de_score"] == pytest.approx(1.0, abs=0.001)


def test_default_bounded_load_still_recovers_direction_fixture(tmp_path):
    out = write_case(tmp_path / "bounded_direction", "direction", 0)
    bundle = load_bundle(out / "prediction.h5ad", out / "target.h5ad", out / "reference.h5ad")
    assert bundle.prediction.shape == (256, 128)
    scored = run_official(bundle, seed=0)
    diagnostic = diagnose(bundle, scored, seed=0)
    planted = {f"synthetic_gene_{i:03d}" for i in [64, 65, 66, 67, 80, 81, 82, 83]}
    _gene_recovery(diagnostic, planted)


def test_aligned_gene_name_order_equivariance():
    case, original_score, original_diagnostic = _evaluate("direction", 0)
    order = np.random.default_rng(918).permutation(128)
    reordered = {name: case[name][:, order].copy() for name in ("prediction", "target", "reference")}
    reordered["planted"] = case["planted"]
    bundle = _bundle(reordered)
    score = run_official(bundle, seed=0)
    diagnostic = diagnose(bundle, score, seed=0)
    before = {row["gene"]: row for row in original_diagnostic["genes"]}
    after = {row["gene"]: row for row in diagnostic["genes"]}
    assert before.keys() == after.keys()
    for gene in before:
        assert after[gene]["absolute_error"] == pytest.approx(before[gene]["absolute_error"], abs=1e-4)
        assert after[gene]["change_pred"] == pytest.approx(before[gene]["change_pred"], abs=1e-4)
    assert score["metrics"]["de_score"] == original_score["metrics"]["de_score"]
    assert score["metrics"]["pb_rel_err"] == pytest.approx(original_score["metrics"]["pb_rel_err"], abs=0.0001)
    # Random index-based sampled variogram pairs change their gene identities;
    # identical finite-sample pair traces under column permutation are not claimed.
