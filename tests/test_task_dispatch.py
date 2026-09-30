"""Pinned public-API parity and explicit participant task roles."""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pytest
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vec_scorelens.cli import main, run, run_task
from vec_scorelens.io import InputPolicyError, load_task_bundle
from vec_scorelens.scorer import SourcePinError, run_official, verify_source_pin


def _inputs(tmp_path):
    rng = np.random.default_rng(73)
    labels = ["a"] * 16 + ["b"] * 16
    reference = rng.uniform(1.5, 2.5, (32, 6)).astype(np.float32)
    response = np.array([0.7, 0.5, 0.8, -0.5, -0.7, -0.6], dtype=np.float32)
    target = reference + response
    target[16:, 0] += 1.5
    prediction = target + rng.uniform(-0.08, 0.08, target.shape).astype(np.float32)
    coords = rng.normal(size=(32, 4)).astype(np.float32)
    paths = {}
    for role, matrix in (("prediction", prediction), ("target", target), ("reference", reference)):
        data = ad.AnnData(matrix)
        data.var_names = [f"g{i}" for i in range(6)]
        data.obs_names = [f"{role}-{i}" for i in range(32)]
        if role == "target":
            data.obs["celltype"] = labels
        data.obsm["spatial_3D"] = coords.copy()
        paths[role] = tmp_path / f"{role}.h5ad"
        data.write_h5ad(paths[role])
    return paths


def _bundle(task, paths):
    contrast = {"wt" if task == "T3" else "reference": paths["reference"]}
    return load_task_bundle(task=task, prediction=paths["prediction"],
                            target=paths["target"], max_cells=32, **contrast)


def _same_metrics(actual, expected):
    assert actual.keys() == expected.keys()
    for name, value in actual.items():
        if value is None or expected[name] is None:
            assert value is expected[name], name
        else:
            np.testing.assert_equal(value, expected[name], err_msg=name)


@pytest.mark.parametrize("task", ["T1", "T2", "T3"])
def test_selected_complete_h5ad_matches_pinned_public_api(tmp_path, task):
    import veckit

    paths = _inputs(tmp_path)
    before = {role: path.read_bytes() for role, path in paths.items()}
    bundle = _bundle(task, paths)
    contrast = {"wt" if task == "T3" else "reference": paths["reference"]}
    with threadpool_limits(limits=1):
        scored = run_official(bundle, seed=3)
        public = veckit.score(task=task, input=paths["prediction"],
                              target=paths["target"], seed=3, **contrast)
    _same_metrics(scored["metrics"], public["metrics"])
    assert scored["task"] == scored["provenance"]["task"] == task
    assert scored["provenance"]["scope"] == "sampled_subset_only"
    assert scored["shape"].__name__ == "common.shape_metrics"
    if task != "T1":
        assert scored["spatial_legacy"].__name__ == "t2_metrics"
        for role in paths:
            coords = getattr(bundle, f"{role}_coords")
            assert coords.shape == (32, 3)
            assert not coords.flags.writeable
    for role, path in paths.items():
        assert path.read_bytes() == before[role]
        assert not getattr(bundle, role).flags.writeable


def test_repeated_cross_task_loader_calls_keep_their_bound_metrics(tmp_path):
    import veckit

    paths = _inputs(tmp_path)
    baselines = {}
    with threadpool_limits(limits=1):
        for task in ("T3", "T1", "T2", "T3", "T2", "T1"):
            bundle = _bundle(task, paths)
            scored = run_official(bundle, seed=11)
            if task in baselines:
                _same_metrics(scored["metrics"], baselines[task])
            baselines[task] = scored["metrics"]
            # The public loader overwrites the same generic sys.modules names.
            other = "T2" if task != "T2" else "T3"
            contrast = {"wt" if other == "T3" else "reference": paths["reference"]}
            veckit.score(task=other, input=paths["prediction"], target=paths["target"], **contrast)
    assert set(baselines) == {"T1", "T2", "T3"}


@pytest.mark.parametrize("task,extra,expected", [
    ("T1", [], "explicit --reference"),
    ("T2", ["--wt", "missing-wt.h5ad"], "--wt is for T3"),
    ("T3", [], "explicit --wt"),
    ("T3", ["--reference", "missing-ref.h5ad"], "--reference is for T1/T2"),
])
def test_cli_role_errors_precede_io_and_are_clean(tmp_path, capsys, task, extra, expected):
    out = tmp_path / "report"
    result = main(["--task", task, "--prediction", "missing-pred.h5ad",
                   "--target", "missing-target.h5ad", "--out", str(out), *extra])
    captured = capsys.readouterr()
    assert result == 2
    assert expected in captured.err
    assert "Traceback" not in captured.err
    assert "does not exist" not in captured.err
    assert not out.exists()


def test_python_role_errors_precede_loader(tmp_path, monkeypatch):
    import vec_scorelens.io as io

    def unexpected_io(**kwargs):
        raise AssertionError("role validation must happen before I/O")
    monkeypatch.setattr(io, "load_task_bundle", unexpected_io)
    with pytest.raises(ValueError, match="--reference is for T1/T2"):
        run_task(task="T3", prediction="missing", target="missing", reference="missing",
                 wt="missing", out=tmp_path / "report")


def test_old_run_forwards_to_t1(tmp_path, monkeypatch):
    import vec_scorelens.cli as cli

    seen = {}
    def capture(**kwargs):
        seen.update(kwargs)
        return "report"
    monkeypatch.setattr(cli, "run_task", capture)
    assert run("prediction", "target", "reference", tmp_path / "report", seed=9) == "report"
    assert seen["task"] == "T1"
    assert seen["reference"] == "reference"
    assert seen["seed"] == 9


def test_source_pin_includes_all_task_and_shape_semantics():
    hashes = verify_source_pin()["source_hashes"]
    assert set(hashes) >= {"T2/metrics.py", "T2/metrics_v2.py", "T3/metrics.py",
                           "T3/metrics_v2.py", "common/shape_metrics.py"}


def test_source_pin_rejects_shape_module_from_another_root(tmp_path, monkeypatch):
    import common.shape_metrics as shape

    monkeypatch.setattr(shape, "__file__", str(tmp_path / "shape_metrics.py"))
    with pytest.raises(SourcePinError, match="common.shape_metrics was imported from"):
        verify_source_pin()


@pytest.mark.parametrize("dependency", ["t1_metrics", "t2_metrics"])
def test_source_pin_rejects_reuse_dependencies_from_another_root(tmp_path, monkeypatch, dependency):
    import score_h5ad

    original_loader = score_h5ad._load_task_metrics
    def wrong_root(task):
        modules = original_loader(task)
        reuse = sys.modules["_reuse"]
        monkeypatch.setattr(getattr(reuse, dependency), "__file__", str(tmp_path / "metrics.py"))
        return modules
    monkeypatch.setattr(score_h5ad, "_load_task_metrics", wrong_root)
    with pytest.raises(SourcePinError, match=f"_reuse.{dependency} was imported from"):
        run_official(_bundle("T3", _inputs(tmp_path)))


@pytest.mark.parametrize("task", ["T2", "T3"])
def test_degenerate_spatial_cloud_has_clean_policy_error(tmp_path, task):
    bundle = _bundle(task, _inputs(tmp_path))
    coords = np.ones_like(bundle.target_coords)
    coords.setflags(write=False)
    with pytest.raises(InputPolicyError, match="target.*zero spatial extent"):
        run_official(replace(bundle, target_coords=coords))
