"""Independent synthetic acceptance through public ScoreLens APIs.

Truth is supplied by the fixture intervention, never inferred from the product.
Recall, precision and healthy false alerts are reported as measured evidence.
"""

from functools import lru_cache
import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from scipy.spatial import cKDTree

from vec_scorelens.diagnostics import diagnose
from vec_scorelens.io import Bundle, load_task_bundle
from vec_scorelens.scorer import run_official
from vec_scorelens.synthetic_multitask import (
    DEFAULT_CASES, MODULES, TASK_CASES, create_task_case, write_task_case,
)


def _bundle(case):
    arrays = []
    coords = []
    for name in ("prediction", "target", "reference"):
        matrix = np.asarray(case[name].X, dtype=np.float32)
        matrix.flags.writeable = False
        arrays.append(matrix)
        coordinate = (np.asarray(case[name].obsm["X_spatial"], dtype=np.float64)
                      if "X_spatial" in case[name].obsm else None)
        if coordinate is not None:
            coordinate.flags.writeable = False
        coords.append(coordinate)
    return Bundle(prediction=arrays[0], target=arrays[1], reference=arrays[2],
                  genes=case["target"].var_names.tolist(),
                  target_labels=case["target"].obs["celltype"].astype(str).to_numpy(),
                  task=case["planted"]["task"], prediction_coords=coords[0],
                  target_coords=coords[1], reference_coords=coords[2],
                  provenance={"synthetic_only": True, "source_shapes": case["planted"]["shapes"],
                              "inputs": {name: {"original_shape": shape} for name, shape in case["planted"]["shapes"].items()},
                              "assessed_scope": "full synthetic fixture"})


def _record(case, scored, diagnostic, **acceptance):
    destination = os.environ.get("VEC_SCORELENS_MULTITASK_EVIDENCE")
    if not destination:
        return
    record = {"task": case["planted"]["task"], "case": case["planted"]["case"],
              "seed": case["planted"]["seed"], "synthetic_only": True,
              "matrix_sha256": case["planted"]["matrix_sha256"],
              "official_metrics": scored["metrics"], "acceptance": acceptance}
    if "summary" in diagnostic.get("response", {}):
        record["response_summary"] = diagnostic["response"]["summary"]
    if "spatial" in diagnostic:
        record["spatial"] = diagnostic["spatial"]
    record["priorities"] = diagnostic.get("priorities", diagnostic.get("priority", {}))
    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(record, allow_nan=False) + "\n")


@lru_cache(maxsize=72)
def _evaluate(task, name, seed):
    case = create_task_case(task, name, seed)
    bundle = _bundle(case)
    scored = run_official(bundle, seed=seed)
    diagnostic = diagnose(bundle, scored, seed=seed)
    return case, scored, diagnostic


def _recovery(rows, planted, key="direct_error"):
    planted = set(planted)
    order = sorted(rows, key=lambda row: (-abs(row[key]), row["gene"]))
    selected = {row["gene"] for row in order[:len(planted)]}
    hits = len(selected & planted)
    return {"precision_at_k": hits / len(selected), "recall_at_k": hits / len(planted),
            "k": len(planted), "top_genes": [row["gene"] for row in order[:len(planted)]]}


@pytest.mark.parametrize("task", ["T2", "T3"])
def test_default_case_and_atomic_file_contract(task, tmp_path):
    sample = create_task_case(task)
    assert sample["planted"]["case"] == DEFAULT_CASES[task]
    out = write_task_case(tmp_path / task, task=task)
    assert set(path.name for path in out.iterdir()) == {
        "prediction.h5ad", "target.h5ad", "reference.h5ad", "truth.json"}
    before = {path.name: path.read_bytes() for path in out.iterdir()}
    with pytest.raises(FileExistsError):
        write_task_case(out, task=task)
    assert before == {path.name: path.read_bytes() for path in out.iterdir()}
    for name in ("prediction", "target", "reference"):
        loaded = ad.read_h5ad(out / f"{name}.h5ad")
        assert loaded.shape == (256, 500)
        assert loaded.obsm["X_spatial"].shape == (256, 3)
        assert ("celltype" in loaded.obs) == (name == "target")
        np.testing.assert_array_equal(loaded.X, sample[name].X)


@pytest.mark.parametrize("task", ["T2", "T3"])
def test_all_fixtures_are_finite_nonnegative_full_width_and_deterministic(task):
    for name in TASK_CASES[task]:
        first, second = create_task_case(task, name, 2), create_task_case(task, name, 2)
        assert first["planted"] == second["planted"]
        for population in ("prediction", "target", "reference"):
            assert first[population].shape == (256, 500)
            assert np.isfinite(first[population].X).all()
            assert np.min(first[population].X) >= 0
            coordinate = first[population].obsm["X_spatial"]
            assert np.isfinite(coordinate).all()
            assert np.linalg.matrix_rank(coordinate - coordinate.mean(axis=0)) == 3
            np.testing.assert_array_equal(first[population].X, second[population].X)
            np.testing.assert_array_equal(coordinate, second[population].obsm["X_spatial"])


@pytest.mark.parametrize("task", ["T1", "T2", "T3"])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_independent_controls_have_ordinary_sampling_noise(task, seed):
    case = create_task_case(task, "independent_healthy", seed)
    assert not case["planted"]["generation"]["within_state_centered"]
    assert case["planted"]["generation"]["random_mixture_counts"]
    assert case["planted"]["causes"] == []
    assert not np.array_equal(case["prediction"].X, case["target"].X)
    # Means must be independently sampled, rather than engineered to match.
    assert np.max(abs(case["prediction"].X.mean(axis=0) - case["target"].X.mean(axis=0))) > 0.005
    labels = np.asarray(case["planted"]["privileged_test_only"]["prediction_state_ids"])
    target_counts = case["target"].obs["celltype"].value_counts().sort_index().to_numpy()
    pred_counts = np.array([np.sum(labels == f"state_{i}") for i in range(4)])
    assert not np.array_equal(pred_counts, target_counts)


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("name", ["wrong_sign", "half_magnitude", "double_magnitude",
                                  "missed_response", "false_positive", "no_response", "mixed"])
def test_t3_gene_interventions_recover_specific_response_evidence(name, seed):
    case, scored, diagnostic = _evaluate("T3", name, seed)
    response = diagnostic["response"]
    rows = {row["gene"]: row for row in response["genes"]}
    truth = case["planted"]
    localized = truth["planted_genes"]
    if name == "mixed":
        localized = [f"synthetic_gene_{i:03d}" for i in truth["sign_gene_indices"] + truth["false_positive_gene_indices"]]
    recovery = _recovery(response["genes"], localized)
    assert recovery["precision_at_k"] >= 0.95
    assert recovery["recall_at_k"] >= 0.95
    if name in ("wrong_sign", "mixed"):
        for index in truth["sign_gene_indices"]:
            row = rows[f"synthetic_gene_{index:03d}"]
            assert row["literal_sign_status"] == "opposite"
            assert row["magnitude_ratio"] is None
            assert abs(row["pred_delta"] / row["truth_delta"]) == pytest.approx(1.0, abs=0.002)
        assert response["summary"]["n_sign_reversals"] >= len(truth["sign_gene_indices"])
    if name in ("half_magnitude", "double_magnitude"):
        expected = 0.5 if name == "half_magnitude" else 2.0
        assert response["summary"]["status"] == "accepted"
        assert response["summary"]["magnitude_median_ratio"] == pytest.approx(expected, abs=0.002)
        assert scored["metrics"]["severity_slope"] == pytest.approx(np.log(expected), abs=0.001)
        assert "response_magnitude" in {finding["rule_id"] for finding in diagnostic["priorities"]["findings"]}
        for index in truth["magnitude_gene_indices"]:
            row = rows[f"synthetic_gene_{index:03d}"]
            assert row["magnitude_ratio"] == pytest.approx(expected, abs=0.002)
            assert row["literal_sign_status"] == "agree"
            assert row["magnitude_status"] == ("under" if name == "half_magnitude" else "over")
    if name in ("missed_response", "no_response"):
        for index in truth["missed_gene_indices"]:
            row = rows[f"synthetic_gene_{index:03d}"]
            assert row["local_de_miss"]
            assert abs(row["pred_delta"]) < 0.002
        assert response["summary"]["n_local_misses"] >= len(truth["missed_gene_indices"])
    if name in ("false_positive", "mixed"):
        for index in truth["false_positive_gene_indices"]:
            assert rows[f"synthetic_gene_{index:03d}"]["local_de_false_positive"]
        assert response["summary"]["n_local_false_positives"] >= len(truth["false_positive_gene_indices"])
    _record(case, scored, diagnostic, **recovery)


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("name", ["wrong_distribution", "mixed"])
def test_t3_distribution_intervention_recovers_prespecified_covariance_modules(seed, name):
    case, scored, diagnostic = _evaluate("T3", name, seed)
    if name == "wrong_distribution":
        assert max(abs(row["direct_error"]) for row in diagnostic["response"]["genes"]) < 1e-4
    # Independent prespecified module statistic, evaluated against reported pairs.
    errors = []
    for module in MODULES[:4]:
        genes = {f"synthetic_gene_{i:03d}" for i in module}
        pairs = [row for row in diagnostic["pairs"] if row["gene_i"] in genes and row["gene_j"] in genes]
        assert pairs, "full-width diagnostic pair selection omitted an entire planted module"
        errors.append(float(np.mean([abs(row["covariance_error"]) for row in pairs])))
    assert sum(value > 0.08 for value in errors) == 4
    _record(case, scored, diagnostic, prespecified_module_covariance_error=errors)


@pytest.mark.parametrize("task", ["T2", "T3"])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_real_h5ad_loader_keeps_all_planted_fixture_rows(task, seed, tmp_path):
    out = write_task_case(tmp_path / f"{task}_{seed}", task=task, seed=seed)
    bundle = load_task_bundle(task=task, prediction=out / "prediction.h5ad",
                              target=out / "target.h5ad", seed=seed,
                              **({"wt": out / "reference.h5ad"} if task == "T3" else {"reference": out / "reference.h5ad"}))
    assert bundle.prediction.shape == (256, 500)
    assert bundle.task == task
    assert bundle.prediction_coords.shape == (256, 3)
    case = create_task_case(task, seed=seed)
    np.testing.assert_allclose(bundle.prediction, case["prediction"].X, rtol=1e-6)
    np.testing.assert_allclose(bundle.prediction_coords, case["prediction"].obsm["X_spatial"], rtol=1e-6)


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("name", ["translation", "rotation", "expansion", "shrinkage"])
def test_t2_rigid_and_uniform_scale_invariants(name, seed):
    _, _, healthy = _evaluate("T2", "healthy", seed)
    case, scored, changed = _evaluate("T2", name, seed)
    before, after = healthy["spatial"], changed["spatial"]
    # D2 uses the same independently sampled pair identities per run. Rigid
    # transforms and uniform scale preserve those normalized distances.
    assert after["shape"]["d2"]["replayed_value"] == pytest.approx(before["shape"]["d2"]["replayed_value"], abs=1e-10)
    assert after["local"]["neighborhood_mmd"]["replayed_value"] == pytest.approx(before["local"]["neighborhood_mmd"]["replayed_value"], abs=1e-8)
    assert after["local"]["moran"]["replayed_value"] == pytest.approx(before["local"]["moran"]["replayed_value"], abs=1e-10)
    expected = {"translation": 0.0, "rotation": 0.0,
                "expansion": np.log(1.7), "shrinkage": np.log(0.55)}[name]
    assert after["growth"]["scale"]["replayed_value"] == pytest.approx(expected, abs=1e-10)
    rule_ids = {row["rule_id"] for row in changed["priorities"]["findings"]}
    assert rule_ids == ({"spatial_size"} if name in ("expansion", "shrinkage") else set())
    # SW/Dice invariance is deliberately NOT asserted: exact healthy controls
    # can expose the independently documented pinned-source PCA frame defect.
    for key in ("sliced_wasserstein", "occupancy"):
        assert after["shape"][key]["status"] == "provisional_frame"
    _record(case, scored, changed, invariants_preserved=["D2", "neighborhood_MMD", "Moran"],
            expected_log_size=expected)


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("name", ["anisotropic", "wrong_geometry", "random_cube", "local_scramble",
                                  "wrong_spatial_assignment", "neighborhood_corruption", "wrong_expression", "mixed"])
def test_t2_planted_geometry_expression_and_assignment_evidence(name, seed):
    _, _, healthy = _evaluate("T2", "healthy", seed)
    case, scored, diagnostic = _evaluate("T2", name, seed)
    spatial = diagnostic["spatial"]
    hspatial = healthy["spatial"]
    acceptance = {}
    if name in ("anisotropic", "wrong_geometry", "random_cube", "mixed"):
        gap = spatial["shape"]["d2"]["replayed_value"] - hspatial["shape"]["d2"]["replayed_value"]
        assert gap > 0.05
        acceptance["D2_increase"] = gap
    if name in ("local_scramble", "wrong_spatial_assignment", "neighborhood_corruption", "random_cube", "mixed"):
        gap = spatial["local"]["neighborhood_mmd"]["replayed_value"] - hspatial["local"]["neighborhood_mmd"]["replayed_value"]
        if name in ("wrong_spatial_assignment", "random_cube", "mixed"):
            assert gap > 0.10
            assert spatial["local"]["moran"]["replayed_value"] < 0.90
        else:
            # Partial interventions need not have the effect size of a global
            # scramble. Independently verify changed neighborhood identities
            # and an actual reported effect rather than sharing its cutoff.
            pc, tc = case["prediction"].obsm["spatial_3D"], case["target"].obsm["spatial_3D"]
            pred_neighbors = cKDTree(pc).query(pc, k=15)[1]
            truth_neighbors = cKDTree(tc).query(tc, k=15)[1]
            changed_fraction = float(np.mean([set(p) != set(t) for p, t in zip(pred_neighbors, truth_neighbors)]))
            assert changed_fraction > 0.25
            assert gap > 1e-5
            assert abs(spatial["local"]["moran"]["replayed_value"] - 1.0) > 1e-5
            acceptance["independent_changed_neighborhood_fraction"] = changed_fraction
        acceptance["neighborhood_MMD_increase"] = gap
        acceptance["Moran_agreement"] = spatial["local"]["moran"]["replayed_value"]
    if name in ("wrong_expression", "mixed"):
        recovery = _recovery(diagnostic["genes"], case["planted"]["planted_genes"], key="absolute_error")
        assert recovery["precision_at_k"] == 1.0
        assert recovery["recall_at_k"] == 1.0
        acceptance.update(recovery)
        assert "mean_residual" in {row["rule_id"] for row in diagnostic["priorities"]["findings"]}
    if name in ("anisotropic", "wrong_geometry", "mixed"):
        assert "spatial_size" in {row["rule_id"] for row in diagnostic["priorities"]["findings"]}
    if name in ("local_scramble", "wrong_spatial_assignment", "neighborhood_corruption"):
        np.testing.assert_array_equal(case["prediction"].X, case["target"].X)
        if name != "neighborhood_corruption":
            assert sorted(map(tuple, case["prediction"].obsm["spatial_3D"])) == sorted(map(tuple, case["target"].obsm["spatial_3D"]))
        assert not {row["rule_id"] for row in diagnostic["priorities"]["findings"]}.intersection({"mean_residual", "opposite_sign", "response_magnitude", "composition_gap", "variance_ratio"})
        acceptance["no_local_priority_threshold_adopted"] = True
    _record(case, scored, diagnostic, **acceptance)


@pytest.mark.parametrize("task", ["T1", "T2", "T3"])
@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("name", ["healthy", "independent_healthy"])
def test_multitask_healthy_controls_have_no_default_followup_alerts(task, seed, name):
    case, scored, diagnostic = _evaluate(task, name, seed)
    priority = diagnostic["priorities"]
    _record(case, scored, diagnostic, false_followup_alerts=priority["n_triggered_rules"],
            tested_default_profile=priority["profile"]["id"],
            max_gene_mean_residual=max(row["absolute_error"] for row in diagnostic["genes"]))
    assert priority["findings"] == []
    assert priority["n_triggered_rules"] == 0
    assert priority["coverage"] == "listed rules evaluated"
    if task in ("T2", "T3"):
        assert abs(diagnostic["spatial"]["growth"]["scale"]["replayed_value"]) < np.log(1.25)
        assert diagnostic["spatial"]["local"]["moran"]["replayed_value"] > 0.95
    if task == "T3":
        assert diagnostic["response"]["summary"]["status"] == "accepted"
        assert diagnostic["response"]["summary"]["magnitude_median_ratio"] == pytest.approx(1.0, abs=0.08)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_t3_named_response_evidence_is_equivariant_to_aligned_gene_order(seed):
    case, _, original = _evaluate("T3", "wrong_sign", seed)
    order = np.random.default_rng(1000 + seed).permutation(500)
    reordered = {name: case[name][:, order].copy() for name in ("prediction", "target", "reference")}
    reordered["planted"] = case["planted"]
    bundle = _bundle(reordered)
    scored = run_official(bundle, seed=seed)
    changed = diagnose(bundle, scored, seed=seed)
    before = {row["gene"]: row for row in original["response"]["genes"]}
    after = {row["gene"]: row for row in changed["response"]["genes"]}
    assert before.keys() == after.keys()
    for gene in before:
        for key in ("truth_delta", "pred_delta", "direct_error"):
            assert after[gene][key] == pytest.approx(before[gene][key], abs=1e-4)
        if before[gene]["truth_de_eligible"]:
            assert after[gene]["literal_sign_status"] == before[gene]["literal_sign_status"]
    assert changed["response"]["summary"]["beta"] == pytest.approx(original["response"]["summary"]["beta"], abs=1e-4)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_t2_joint_expression_coordinate_row_permutation_keeps_local_evidence(seed):
    case, _, original = _evaluate("T2", "healthy", seed)
    rows = np.random.default_rng(2500 + seed).permutation(256)
    reordered = {name: case[name].copy() for name in ("prediction", "target", "reference")}
    reordered["prediction"] = case["prediction"][rows].copy()
    reordered["planted"] = case["planted"]
    bundle = _bundle(reordered)
    scored = run_official(bundle, seed=seed)
    changed = diagnose(bundle, scored, seed=seed)
    assert max(row["absolute_error"] for row in changed["genes"]) < 1e-4
    for key in ("neighborhood_mmd", "moran"):
        assert changed["spatial"]["local"][key]["replayed_value"] == pytest.approx(original["spatial"]["local"][key]["replayed_value"], abs=1e-6)
    assert changed["priorities"]["findings"] == []
    # D2's independent sampled realization and provisional SW/Dice frame can
    # change with row order. No exact-score invariance claim is made for these.
