"""Spatial source conformance and scientific invariance regressions."""
from types import SimpleNamespace
import json

import numpy as np
import pytest

from vec_scorelens.spatial import diagnose_spatial, _shape_grid, _d2


@pytest.fixture(scope="module")
def modules():
    import score_h5ad
    from common import core_metrics, shape_metrics
    legacy, _ = score_h5ad._load_task_metrics("T2")
    return {"core": core_metrics, "shape": shape_metrics, "spatial_legacy": legacy}


def fixture():
    rng = np.random.default_rng(71)
    c = rng.normal(size=(96, 3))*[3, 1, .4]
    x = (3+np.tanh(c[:, 0, None])*np.linspace(.2, 1, 36)+rng.normal(scale=.2, size=(96, 36))).astype(np.float32)
    return SimpleNamespace(prediction=x.copy(), target=x, reference=x.copy(),
        prediction_coords=c.copy(), target_coords=c, reference_coords=c.copy(),
        genes=[f"g{i}" for i in range(36)],
        provenance={"inputs": {"prediction": {"original_shape": [500, 36]}, "target": {"original_shape": [1000, 36]}}})


def scored(bundle, modules, seed=0):
    p, t = bundle.prediction_coords, bundle.target_coords
    s, l = modules["shape"], modules["spatial_legacy"]
    return {**modules, "metrics": {
        "d2_shape": round(s.d2_distance(p, t, seed=seed), 5),
        "sliced_wasserstein": round(s.sliced_wasserstein(p, t, seed=seed)[0], 5),
        "occupancy_dice": round(s.occupancy_dice(p, t, seed=seed)[0], 4),
        "scale_log_ratio": round(s.scale_log_ratio(p, t), 4),
        "count_log_ratio": round(s.count_log_ratio(len(p), len(t)), 4),
        "neighborhood_mmd": round(l.neighborhood_mmd(bundle.prediction, p, bundle.target, t), 5),
        "morans_I_agreement": round(l.morans_I_agreement(bundle.prediction, p, bundle.target, t), 4),
    }}


def test_spatial_allocations_reconcile_actual_pinned_source_and_preserve_inputs(modules):
    b = fixture()
    originals = [a.copy() for a in (b.prediction, b.target, b.prediction_coords, b.target_coords)]
    d = diagnose_spatial(b, scored(b, modules, seed=5), seed=5)
    json.dumps(d, allow_nan=False)
    for family in (d["shape"]["d2"], d["shape"]["sliced_wasserstein"], d["shape"]["occupancy"], d["local"]["moran"]):
        assert family["reconciliation_abs_error"] < 1e-12
        assert family["rounded_official_abs_error"] <= 5e-5
    assert sum(r["contribution"] for r in d["shape"]["d2"]["rows"]) == pytest.approx(d["shape"]["d2"]["replayed_value"], abs=1e-12)
    assert sum(r["loss_contribution"] for r in d["shape"]["occupancy"]["rows"]) == pytest.approx(1-d["shape"]["occupancy"]["replayed_value"])
    assert sum(r["agreement_contribution"] for r in d["local"]["moran"]["rows"]) == pytest.approx(1)
    for old, current in zip(originals, (b.prediction, b.target, b.prediction_coords, b.target_coords)):
        np.testing.assert_array_equal(old, current)


def test_count_scope_must_not_hide_original_counts(modules):
    b = fixture()
    d = diagnose_spatial(b, {**modules, "metrics": {}})
    counts = d["growth"]["count"]
    assert counts["replayed_value"] == 0
    assert counts["original_prediction_count"] == 500
    assert counts["original_target_count"] == 1000
    assert counts["original_log_ratio"] == pytest.approx(np.log(.5))
    assert counts["original_relative_count_difference"] == -.5


def test_local_scramble_sensitive_and_negative_identity_is_retained(modules):
    b = fixture()
    d = diagnose_spatial(b, {**modules, "metrics": {}}, seed=71)
    local = d["local"]["neighborhood_mmd"]
    assert local["effective_seed"] == 0
    assert local["replayed_value"] < 0
    assert local["sensitivity"]["after"] > local["replayed_value"]+.1
    assert local["sensitivity"]["evidence_type"] == "noncausal_sensitivity"
    assert local["includes_self_for_unique_coordinates"] is True


def test_rigid_and_scale_invariance_where_source_actually_has_it(modules):
    b = fixture()
    s, l = modules["shape"], modules["spatial_legacy"]
    rng = np.random.default_rng(333)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    q[:, 0] *= np.linalg.det(q)
    c = b.target_coords
    baseline_d2 = s.d2_distance(c, c)
    baseline_local = l.neighborhood_mmd(b.prediction, c, b.target, c)
    for transform, expected_scale in ((c+[2, -1, 7], 0), (c@q, 0), (c*2, np.log(2))):
        assert s.d2_distance(transform, c) == pytest.approx(baseline_d2, abs=1e-12)
        assert s.scale_log_ratio(transform, c) == pytest.approx(expected_scale, abs=1e-12)
        assert l.neighborhood_mmd(b.prediction, transform, b.target, c) == pytest.approx(baseline_local, abs=1e-7)
        assert l.morans_I_agreement(b.prediction, transform, b.target, c) == pytest.approx(1)


def test_known_identical_cloud_and_rotation_frame_defect_is_exposed(modules):
    rng = np.random.default_rng(219)
    c = np.column_stack([rng.exponential(2, 384), rng.exponential(.8, 384), rng.exponential(.3, 384)])
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    identity_sw, identity_dice = _shape_grid(c, c, modules["shape"], {"metrics": {}}, 0)
    rotated_sw, rotated_dice = _shape_grid(c@q, c, modules["shape"], {"metrics": {}}, 0)
    # This is a regression for the PINNED implementation and numerical runtime;
    # changing the accepted upstream requires explicitly revisiting this test.
    assert identity_sw["replayed_value"] > .01
    assert identity_dice["replayed_value"] < .9
    assert rotated_sw["replayed_value"] < 1e-12
    assert rotated_dice["replayed_value"] == 1
    assert identity_sw["status"] == identity_dice["status"] == "provisional_frame"
    assert "handedness defect" in identity_sw["limitations"][0]


def test_d2_degenerate_branch_replays_rng_and_empty_interval_case(modules):
    b = fixture()
    zeros = np.zeros_like(b.target_coords)
    for p, t in ((zeros, zeros), (zeros, b.target_coords), (b.target_coords, zeros)):
        d = _d2(p, t, modules["shape"], {"metrics": {}}, 3)
        assert d["reconciliation_abs_error"] < 1e-12
    assert _d2(zeros, zeros, modules["shape"], {"metrics": {}}, 0)["rows"] == []


def test_moran_constant_profiles_report_undefined(modules):
    b = fixture()
    b.prediction = np.full_like(b.prediction, 3)
    from vec_scorelens.spatial import _moran
    with np.errstate(invalid="ignore", divide="ignore"):
        d = _moran(b, modules["spatial_legacy"], {"metrics": {}})
    assert d["status"] == "undefined"
    assert not np.isfinite(d["replayed_value"])


def test_resource_guard_before_upstream_allocation(modules):
    b = fixture()
    b.prediction = np.zeros((96, 2001), dtype=np.float32)
    from vec_scorelens.io import InputPolicyError
    with pytest.raises(InputPolicyError, match="2000 genes"):
        diagnose_spatial(b, {**modules, "metrics": {}})


def test_anchor_terms_match_independent_dense_kernel_oracle_and_actual_coordinates(modules):
    from scipy.spatial.distance import cdist
    from sklearn.decomposition import PCA
    from vec_scorelens.spatial import _mmd_anchors

    b = fixture()
    # Unequal counts exercise the two distinct denominators. The truth fit is
    # unchanged and includes every truth neighborhood before kernel sampling.
    b.prediction = b.prediction[:67].copy()
    b.prediction[:, :8] += .3
    b.prediction_coords = b.prediction_coords[:67].copy()
    b.provenance["selection"] = {"rows": {
        "prediction": [{"position": 3*i, "obs_name": f"cell-{i}"} for i in range(67)],
        "target": [{"position": 2*i, "obs_name": f"truth-{i}"} for i in range(96)]}}
    legacy, core = modules["spatial_legacy"], modules["core"]
    hp = legacy._knn_neighborhood_pb(b.prediction, b.prediction_coords)
    ht = legacy._knn_neighborhood_pb(b.target, b.target_coords)
    originals = (hp.copy(), ht.copy(), b.prediction_coords.copy(), b.target_coords.copy())
    for array in (hp, ht, b.prediction_coords, b.target_coords):
        array.setflags(write=False)
    official = core.mmd_unbiased(hp, ht, seed=0)
    rows, info = _mmd_anchors(hp, ht, b, official)
    assert info["passed"]
    assert info["source_replay_abs_error"] < 1e-12
    assert sum(row["allocation"] for row in rows) == pytest.approx(official, abs=info["tolerance"])
    assert len(rows) == 67+96

    rng = np.random.default_rng(0)
    pi, ti = rng.choice(67, 67, replace=False), rng.choice(96, 96, replace=False)
    pca = PCA(n_components=30, random_state=0).fit(ht)
    a, t = pca.transform(hp[pi]), pca.transform(ht[ti])
    truth_distances = ((t[:, None]-t[None, :])**2).sum(-1)
    gamma = 1/(np.median(truth_distances[truth_distances > 0])+1e-9)
    # cdist + explicit Gaussian is independent of sklearn's rbf_kernel path.
    pp, tt, pt = cdist(a, a, "sqeuclidean"), cdist(t, t, "sqeuclidean"), cdist(a, t, "sqeuclidean")
    expected_p, expected_t = np.zeros(67), np.zeros(96)
    for scale in (.25, .5, 1, 2, 4):
        kpp, ktt, kpt = np.exp(-gamma*scale*pp), np.exp(-gamma*scale*tt), np.exp(-gamma*scale*pt)
        np.fill_diagonal(kpp, 0)
        np.fill_diagonal(ktt, 0)
        expected_p += (kpp.sum(1)/(67*66)-kpt.sum(1)/(67*96))/5
        expected_t += (ktt.sum(1)/(96*95)-kpt.sum(0)/(67*96))/5
    np.testing.assert_allclose([r["allocation"] for r in rows[:67]], expected_p, atol=2e-8, rtol=0)
    np.testing.assert_allclose([r["allocation"] for r in rows[67:]], expected_t, atol=2e-8, rtol=0)
    for row in rows:
        coords = b.prediction_coords if row["side"] == "prediction" else b.target_coords
        multiplier = 3 if row["side"] == "prediction" else 2
        np.testing.assert_array_equal(row["coordinates"], coords[row["row_position"]])
        assert row["original_row_position"] == multiplier*row["row_position"]
    for old, current in zip(originals, (hp, ht, b.prediction_coords, b.target_coords)):
        np.testing.assert_array_equal(old, current)


def test_anchor_identity_keeps_negative_terms_and_reconciles(modules):
    b = fixture()
    result = diagnose_spatial(b, {**modules, "metrics": {}})["local"]["neighborhood_mmd"]
    assert result["anchor_allocation"]["passed"]
    assert result["anchor_allocation"]["sum"] < 0
    assert any(r["allocation"] < 0 for r in result["anchor_rows"])
    assert {r["side"] for r in result["anchor_rows"]} == {"prediction", "truth"}
    assert "nonadditive" in result["reconciliation_status"]


def test_copermuting_expression_and_coords_preserves_graph_and_anchor_terms(modules):
    from vec_scorelens.spatial import _mmd_anchors
    b = fixture()
    b.prediction[:, 0] += .25
    legacy, core = modules["spatial_legacy"], modules["core"]
    hp = legacy._knn_neighborhood_pb(b.prediction, b.prediction_coords)
    ht = legacy._knn_neighborhood_pb(b.target, b.target_coords)
    rows, info = _mmd_anchors(hp, ht, b, core.mmd_unbiased(hp, ht))
    perm = np.random.default_rng(44).permutation(len(b.prediction))
    b.prediction = b.prediction[perm]
    b.prediction_coords = b.prediction_coords[perm]
    permuted_hp = legacy._knn_neighborhood_pb(b.prediction, b.prediction_coords)
    np.testing.assert_allclose(permuted_hp, hp[perm], atol=1e-7, rtol=0)
    permuted_rows, permuted_info = _mmd_anchors(permuted_hp, ht, b, core.mmd_unbiased(permuted_hp, ht))
    assert info["passed"] and permuted_info["passed"]
    first = {(r["side"], r["row_position"]): r["allocation"] for r in rows}
    for row in permuted_rows:
        original_position = int(perm[row["row_position"]]) if row["side"] == "prediction" else row["row_position"]
        assert row["allocation"] == pytest.approx(first[row["side"], original_position], abs=2e-8)
