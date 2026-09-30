from __future__ import annotations

import json
import importlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from vec_scorelens.report import write_report


ROOT = Path(__file__).resolve().parents[1]


def _bundle(**provenance):
    return SimpleNamespace(provenance=provenance)


def _diagnostics(**overrides):
    result = {
        "genes": [{"gene": "GATA3", "change_truth": 0.3, "change_pred": -0.1, "absolute_error": 0.4, "squared_error": 0.16,
                   "variance_pred": 1, "variance_truth": 1, "truth_de_direction": "up", "pred_local_de_direction": "down", "official_slot_hit": False}],
        "gene_summary": {"reconciliation_abs_error": 0.0, "de_guard": "defined"},
        "composition": [{"celltype": "T cell", "target_fraction": 0.7, "pred_soft_fraction": 0.4, "proportion_error": -0.3, "jsd_term": 0.02}],
        "composition_summary": {"reconciliation_abs_error": 0.0},
        "pairs": [{"gene_i": "GATA3", "gene_j": "IL7R", "sampled_occurrences": 4, "variogram_squared_error": 0.12, "covariance_error": -0.2}],
        "pair_summary": {"n_sampled_pairs": 20_000, "n_unique_ordered_pairs": 15, "reconciliation_abs_error": 0.0},
        "distribution": {"mmd_sensitivity": {"before": 0.5, "after": 0.4, "change": -0.1, "affected_genes": ["GATA3"]}},
        "warnings": [],
    }
    result.update(overrides)
    return result


def _scored():
    return {"metrics": {"pb_rel_err": 0.11, "mmd": 0.23}, "provenance": {"commit": "pinned"}}


def test_write_report_creates_required_private_artifacts(tmp_path):
    report = write_report(_bundle(source_path=str(tmp_path / "input.h5ad")), _scored(), _diagnostics(), tmp_path / "report")

    assert report == tmp_path / "report"
    assert {path.name for path in report.iterdir()} == {
        "summary.md", "diagnostics.json", "gene_errors.csv", "composition_errors.csv",
        "covariance_errors.csv", "provenance.json", "report.html", "figures",
    }
    assert {path.name for path in (report / "figures").iterdir()} == {
        "gene_changes.png", "composition.png", "pair_discrepancies.png"}
    page = (report / "report.html").read_text()
    assert "<table>" in page and "<figure>" in page and "<pre>" not in page
    assert "src='figures/gene_changes.png'" in page
    assert "Within listed thresholds" in page
    summary = (report / "summary.md").read_text()
    assert "selected subset" in summary
    assert "different scales" in summary
    assert "false-positive means a prediction-only" in summary
    assert "target-aware" in summary
    assert "GATA3" in summary and "T cell" in summary and "IL7R" in summary
    covariance = (report / "covariance_errors.csv").read_text()
    assert "variogram_squared_error" in covariance and "covariance_error" in covariance
    provenance = json.loads((report / "provenance.json").read_text())
    assert provenance["scorer"]["commit"] == "pinned"
    assert json.loads((report / "diagnostics.json").read_text())["official_metrics"] == {"pb_rel_err": 0.11, "mmd": 0.23}


def test_write_report_replaces_nonfinite_with_null_and_warning(tmp_path):
    scored = {"metrics": {"mmd": float("nan")}, "provenance": {}}
    report = write_report(_bundle(), scored, _diagnostics(genes=[{"gene": "A", "squared_error": float("nan")}]), tmp_path / "report")

    raw = (report / "diagnostics.json").read_text()
    parsed = json.loads(raw)
    assert "NaN" not in raw
    assert parsed["genes"][0]["squared_error"] is None
    assert parsed["official_metrics"]["mmd"] is None
    assert any("non-finite value" in warning for warning in parsed["warnings"])


def test_write_report_escapes_adversarial_labels_for_all_renderings(tmp_path):
    bad_gene = '=HYPERLINK("https://example.invalid","click")'
    bad_class = "<img src=x onerror=alert(1)>"
    report = write_report(
        _bundle(), _scored(),
        _diagnostics(
            genes=[{"gene": bad_gene, "absolute_error": 1, "squared_error": 1}],
            composition=[{"celltype": bad_class, "proportion_error": 1}],
            pairs=[{"gene_i": bad_gene, "gene_j": bad_class, "variogram_squared_error": 1}],
        ),
        tmp_path / "report",
    )

    markdown = (report / "summary.md").read_text()
    page = (report / "report.html").read_text()
    csv = (report / "gene_errors.csv").read_text()
    assert "&lt;img" in markdown
    assert "&lt;img" in page
    assert "<img src=x" not in page
    assert "'=HYPERLINK" in csv
    assert "<script" not in page.lower()


def test_write_report_merges_and_deduplicates_provenance_sampling_warnings(tmp_path):
    omitted = "Uniform target sampling omitted valid full-file celltype classes: rare_type"
    report = write_report(
        _bundle(warnings=[omitted, omitted]), _scored(),
        _diagnostics(warnings=["diagnostic warning", omitted],
                     genes=[{"gene": "GATA3", "squared_error": .16}]), tmp_path / "report",
    )

    written = json.loads((report / "diagnostics.json").read_text())["warnings"]
    assert written[:2] == ["diagnostic warning", omitted]
    assert written.count(omitted) == 1
    assert len(written) == 3
    assert written[2].startswith("Measured figure evidence unavailable: figures/gene_changes.png.")
    assert "Uniform target sampling omitted valid full\\-file celltype classes: rare\\_type" in (report / "summary.md").read_text()
    assert omitted in (report / "report.html").read_text()
    assert "Measured figure evidence unavailable" in (report / "report.html").read_text()


def test_write_report_refuses_collision_and_source_overlap(tmp_path):
    existing = tmp_path / "existing"
    existing.mkdir()
    try:
        write_report(_bundle(), _scored(), _diagnostics(), existing)
    except FileExistsError:
        pass
    else:
        raise AssertionError("existing directory was overwritten")

    symlink = tmp_path / "symlink"
    symlink.symlink_to(existing, target_is_directory=True)
    try:
        write_report(_bundle(), _scored(), _diagnostics(), symlink)
    except FileExistsError:
        pass
    else:
        raise AssertionError("symlink destination was followed")

    real_parent = tmp_path / "real-parent"
    real_parent.mkdir()
    symlink_parent = tmp_path / "symlink-parent"
    symlink_parent.symlink_to(real_parent, target_is_directory=True)
    try:
        write_report(_bundle(), _scored(), _diagnostics(), symlink_parent / "report")
    except ValueError as exc:
        assert "symlink" in str(exc)
    else:
        raise AssertionError("symlink parent was followed")

    source = tmp_path / "input"
    source.mkdir()
    try:
        write_report(_bundle(input_path=str(source)), _scored(), _diagnostics(), source / "report")
    except ValueError as exc:
        assert "overlaps" in str(exc)
    else:
        raise AssertionError("report output was accepted inside a source directory")


def test_write_report_removes_its_staging_directory_on_failure(tmp_path, monkeypatch):
    import vec_scorelens.report as reports

    def fail(*args, **kwargs):
        raise OSError("injected write failure")

    monkeypatch.setattr(reports, "_write_text", fail)
    try:
        write_report(_bundle(), _scored(), _diagnostics(), tmp_path / "report")
    except OSError:
        pass
    else:
        raise AssertionError("injected failure did not propagate")
    assert not (tmp_path / "report").exists()
    assert not list(tmp_path.glob(".report.stage-*"))


def test_report_renders_actual_diagnostics_schema_without_null_field_mismatch(tmp_path, monkeypatch):
    """Exercise the real diagnose() output, rather than only a report-shaped mock."""
    score_h5ad = importlib.import_module("score_h5ad")
    core = importlib.import_module("common.core_metrics")
    legacy, metrics_v2 = score_h5ad._load_task_metrics("T1")
    from vec_scorelens.diagnostics import diagnose
    from vec_scorelens.io import load_bundle
    from vec_scorelens.synthetic import write_case

    inputs = write_case(tmp_path / "inputs", "direction", seed=0)
    bundle = load_bundle(inputs / "prediction.h5ad", inputs / "target.h5ad", inputs / "reference.h5ad", max_cells=64, seed=0)
    probe = legacy.train_frozen_probe(bundle.target, bundle.target_labels)
    metrics = metrics_v2.score_task1_v2(
        bundle.prediction, ["NA"] * len(bundle.prediction), bundle.target,
        bundle.target_labels, bundle.reference, probe=probe, seed=0,
    )
    scored = {"metrics": metrics, "probe": probe, "core": core, "legacy": legacy, "provenance": {"commit": "pinned"}}
    result = diagnose(bundle, scored, seed=0)
    report = write_report(bundle, scored, result, tmp_path / "report")

    summary = (report / "summary.md").read_text()
    assert "Target local DE calls:" in summary
    assert "Distribution variance-ratio evidence" in summary
    assert "Pseudobulk allocation reconciliation error: null" not in summary
    assert "Sampled pair count:" in summary
    assert "Possible ordered nonself pairs:" in summary
    assert "The listed delta is after minus before." in summary
    assert "synthetic_gene_" in (report / "gene_errors.csv").read_text()
    assert "variogram_contribution" in (report / "covariance_errors.csv").read_text()
    assert json.loads((report / "diagnostics.json").read_text())["distribution"]["mmd_sensitivity"]["change"] is not None


@pytest.mark.parametrize("task", ["T2", "T3"])
def test_spatial_report_figures_and_task_tables_use_supplied_evidence(tmp_path, task):
    coords = np.array([[0, 0, 0], [1, 0, 0], [0, 2, 0], [0, 0, 3]], dtype=float)
    bundle = SimpleNamespace(task=task, provenance={}, prediction_coords=coords, target_coords=coords)
    diagnostics = _diagnostics(
        scope={"selected_counts": {"prediction": 4, "target": 4}, "original_counts": {"prediction": 400, "target": 500}},
        spatial={"shape": {"d2": {"definition": "Sampled D2", "status": "defined", "official_value": .02, "replayed_value": .02001, "rows": []}},
                 "growth": {"scale": {"replayed_value": 0, "prediction_rms": 1, "target_rms": 1},
                            "count": {"prediction_count": 4, "target_count": 4, "original_prediction_count": 400, "original_target_count": 500, "original_log_ratio": -.22314}},
                 "local": {"neighborhood_mmd": {"anchor_rows": [
                     {"side": side, "row_position": 0, "coordinates": [0, 0, 0], "allocation": term, "obs_name": "<script>private</script>"}
                     for side, term in (("prediction", -.001), ("truth", .002))], "anchor_allocation": {"sum": .001, "passed": True}},
                           "moran": {"official_panel_includes_metric": task == "T2", "replayed_value": .9, "rows": []}}},
        response={"genes": [{"gene": "g", "truth_delta": .5, "pred_delta": .8, "truth_de_eligible": True,
                            "magnitude_ratio": 1.6, "magnitude_status": "over", "direct_error": .3}],
                  "summary": {"status": "accepted", "official_logbeta": .47, "beta": 1.6,
                              "official_r2": 1, "n_sign_reversals": 1}})
    scored = {"metrics": {"count_log_ratio": 0, "mmd_u": -.01}, "provenance": {"task": task, "seed": 1}}
    report = write_report(bundle, scored, diagnostics, tmp_path / task)
    page = (report / "report.html").read_text()
    assert "<pre>" not in page and "<script>" not in page
    assert "&lt;script&gt;private&lt;/script&gt;" in page
    assert "Signed anchor kernel terms" in page and "not causal per-cell errors" in page
    assert "source-frame caveat" in page and "Laterality is unassessed" in page
    assert (report / "spatial_errors.csv").is_file()
    assert (report / "neighborhood_anchors.csv").is_file()
    manifest = json.loads((report / "diagnostics.json").read_text())["figures"]
    assert {item["path"] for item in manifest} >= {"figures/spatial_clouds.png", "figures/shape_distances.png", "figures/spatial_size_counts.png", "figures/local_neighborhood.png"}
    for item in manifest:
        assert (report / item["path"]).read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
        assert item["status"] == "defined"
    if task == "T3":
        assert (report / "response_errors.csv").is_file()
        assert (report / "figures" / "response.png").is_file()
        assert "constant nonzero uniform offset" in page
        assert "n_sign_reversals" in page
        assert "including genes without a target DE call" in page
        assert "including genes without a target DE call" in (report / "summary.md").read_text()
        assert "not part of the official T3 panel (extra diagnostic only)" in page
    else:
        assert not (report / "response_errors.csv").exists()
    assert json.loads((report / "diagnostics.json").read_text())["official_metrics"] == scored["metrics"]


def test_figure_failure_rolls_back_report_and_staging(tmp_path, monkeypatch):
    import vec_scorelens.figures as figures
    def fail(*args, **kwargs):
        raise RuntimeError("injected mandatory figure failure")
    monkeypatch.setattr(figures, "write_figures", fail)
    with pytest.raises(RuntimeError, match="mandatory figure failure"):
        write_report(_bundle(), _scored(), _diagnostics(), tmp_path / "report")
    assert not (tmp_path / "report").exists()
    assert not list(tmp_path.glob(".report.stage-*"))


def test_matplotlib_labels_are_literal_untrusted_text(tmp_path):
    diagnostics = _diagnostics()
    diagnostics["genes"][0]["gene"] = r"$\untrusted_command$"
    diagnostics["composition"][0]["celltype"] = r"$\not_a_math_command$"
    report = write_report(_bundle(), _scored(), diagnostics, tmp_path / "report")
    assert (report / "figures" / "gene_changes.png").is_file()
    assert (report / "figures" / "composition.png").is_file()
