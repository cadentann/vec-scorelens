"""Independent identities and fidelity checks for the diagnostic evidence."""

from types import SimpleNamespace

import json
import numpy as np
import pytest

from vec_scorelens.diagnostics import diagnose


@pytest.fixture(scope="module")
def upstream():
    import score_h5ad
    from common import core_metrics

    legacy, v2 = score_h5ad._load_task_metrics("T1")
    return core_metrics, legacy, v2


def make_bundle(*, mode="errors", seed=71):
    rng = np.random.default_rng(seed)
    n, G = 64, 36
    baseline = 3.0 + np.linspace(0.0, 1.0, G)
    R = (baseline + rng.normal(0, 0.025, (n, G))).astype(np.float32)
    delta = np.zeros(G, np.float32)
    delta[[1, 7, 15, 23]] = 0.6
    delta[[4, 10, 18, 29]] = -0.6
    T = R + delta
    P = T.copy()
    if mode == "errors":
        P[:, 5] += 0.7
        P[:, 10] += 0.9
        P[:, 27] -= 0.35
    elif mode == "negative_slots":
        shifted = np.linspace(-0.5, -0.9, G).astype(np.float32)
        shifted[[1, 7, 15, 23]] = [-0.10, -0.11, -0.12, -0.13]
        shifted[[4, 10, 18, 29]] = [-1.5, -1.6, -1.7, -1.8]
        P = R + shifted
    elif mode == "no_change":
        P = R.copy()
    elif mode == "no_truth_de":
        T = R.copy()
        P = R.copy()
    elif mode == "collapse":
        P = np.broadcast_to(T.mean(0), T.shape).copy()
    for X in (P, T, R):
        X.setflags(write=False)
    return SimpleNamespace(
        prediction=P, target=T, reference=R,
        genes=[f"g{i:02}" for i in range(G)],
        target_labels=np.array(["a"] * 60 + ["b"] * 4),
        provenance={"scope": "synthetic selected subset"},
    )


class FixedProbe:
    classes_ = np.array(["a", "b"])

    def predict_proba(self, X):
        return np.tile(np.array([0.05, 0.95]), (len(X), 1))


def score_fixture(bundle, upstream, seed=0, probe=None):
    core, legacy, v2 = upstream
    probe = probe if probe is not None else legacy.train_frozen_probe(bundle.target, bundle.target_labels)
    metrics = v2.score_task1_v2(
        bundle.prediction, np.full(len(bundle.prediction), "NA"),
        bundle.target, bundle.target_labels, bundle.reference,
        probe=probe, seed=seed,
    )
    return {"core": core, "legacy": legacy, "probe": probe, "metrics": metrics}


def test_gene_allocation_reconstructs_squared_norm_and_official_de(upstream):
    bundle = make_bundle()
    scored = score_fixture(bundle, upstream)
    diagnostic = diagnose(bundle, scored)
    rows, summary = diagnostic["genes"], diagnostic["gene_summary"]
    mu_p = bundle.prediction.mean(0)
    mu_t = bundle.target.mean(0)
    expected_squared_relative_norm = float(np.linalg.norm(mu_p - mu_t) / np.linalg.norm(mu_t)) ** 2
    assert sum(r["pb_squared_allocation"] for r in rows) == pytest.approx(expected_squared_relative_norm, abs=1e-8)
    assert summary["de_guard"] == "none"
    assert sum(r["de_score_allocation"] for r in rows) == pytest.approx(summary["official_de_score"], abs=1e-12)
    assert sum(r["de_raw_allocation"] for r in rows) == pytest.approx(summary["official_de_raw"], abs=1e-12)
    assert summary["pb_reconciliation_abs_error"] < 1e-8
    # A gene's large magnitude mismatch need not belong to the truth DE set.
    row = next(r for r in rows if r["gene"] == "g05")
    assert row["truth_de_direction"] == "none"
    assert row["pred_local_de_direction"] == "up"
    assert row["local_de_false_positive"] is True
    assert row["absolute_error"] > 0.65
    assert row["de_score_allocation"] == 0


def test_rank_slot_hit_is_not_literal_sign_agreement(upstream):
    bundle = make_bundle(mode="negative_slots")
    result = diagnose(bundle, score_fixture(bundle, upstream))
    true_up = [r for r in result["genes"] if r["truth_de_direction"] == "up"]
    assert len(true_up) == 4
    assert all(r["official_slot_hit"] for r in true_up)
    assert all(r["literal_pred_direction"] == "down" for r in true_up)
    assert all(r["slot_hit_wrong_literal_sign"] for r in true_up)
    assert all(not r["local_de_hit"] for r in true_up)
    assert result["gene_summary"]["official_de_raw"] == 1.0


@pytest.mark.parametrize("mode,guard", [
    ("no_change", "no_change_std_guard"),
    ("no_truth_de", "undefined_no_truth_de"),
])
def test_zero_and_undefined_de_branches_are_explicit_and_json_safe(upstream, mode, guard):
    bundle = make_bundle(mode=mode)
    result = diagnose(bundle, score_fixture(bundle, upstream))
    summary = result["gene_summary"]
    assert summary["de_guard"] == guard
    if mode == "no_change":
        assert summary["official_de_score"] == 0
        assert all(r["de_score_allocation"] == 0 for r in result["genes"])
        assert all(r["de_raw_allocation"] == 0 for r in result["genes"])
    else:
        assert summary["official_de_score"] is None
        assert all(r["de_score_allocation"] is None for r in result["genes"])
        assert result["distribution"]["mmd_sensitivity"]["affected_genes"] == []
        assert result["distribution"]["mmd_sensitivity"]["change"] == 0
    json.dumps(result, allow_nan=False)


def test_composition_terms_use_count_vs_fraction_smoothing_exactly(upstream):
    bundle = make_bundle()
    result = diagnose(bundle, score_fixture(bundle, upstream, probe=FixedProbe()))
    row = result["composition"][0]
    expected_count_smoothing = (60 + 1e-6) / (64 + 2e-6)
    wrong_fraction_smoothing = (60 / 64 + 1e-6) / (1 + 2e-6)
    assert row["smoothed_target"] == pytest.approx(expected_count_smoothing, abs=1e-15)
    assert abs(row["smoothed_target"] - wrong_fraction_smoothing) > 1e-7
    assert row["smoothed_pred"] == pytest.approx((0.05 + 1e-6) / (1 + 2e-6), abs=1e-15)
    summary = result["composition_summary"]
    assert sum(r["jsd_term"] for r in result["composition"]) == pytest.approx(summary["jsd"], abs=1e-15)
    assert summary["reconciliation_abs_error"] < 1e-15


@pytest.mark.parametrize("seed", [0, 3])
def test_chunked_variogram_preserves_rng_and_duplicate_pair_weights(upstream, seed):
    core, _, _ = upstream
    bundle = make_bundle()
    result = diagnose(bundle, score_fixture(bundle, upstream, seed=seed), seed=seed)
    summary = result["pair_summary"]
    unrounded = core.variogram_score(bundle.prediction, bundle.target, seed=seed)
    assert summary["variogram"] == pytest.approx(unrounded, abs=2e-7)
    assert sum(r["variogram_contribution"] for r in result["pairs"]) == pytest.approx(unrounded, abs=2e-7)
    assert sum(r["sampled_occurrences"] for r in result["pairs"]) == summary["n_sampled_pairs"]
    assert any(r["sampled_occurrences"] > 1 for r in result["pairs"])
    assert all(r["gene_i"] != r["gene_j"] for r in result["pairs"])
    # Independently check a descriptive covariance column against a 2x2
    # covariance calculation rather than against the variogram formula.
    row = result["pairs"][0]
    g, h = bundle.genes.index(row["gene_i"]), bundle.genes.index(row["gene_j"])
    selected = bundle.prediction[summary["sampled_prediction_rows"]]
    expected_covariance = np.cov(selected[:, [g, h]].astype(float), rowvar=False, ddof=0)[0, 1]
    assert row["covariance_pred"] == pytest.approx(expected_covariance, abs=1e-12)


def test_mmd_sensitivity_is_exact_paired_replay_and_preserves_inputs(upstream):
    core, _, _ = upstream
    bundle = make_bundle()
    originals = tuple(X.tobytes() for X in (bundle.prediction, bundle.target, bundle.reference))
    result = diagnose(bundle, score_fixture(bundle, upstream, seed=2), seed=2)
    sensitivity = result["distribution"]["mmd_sensitivity"]
    cf = bundle.prediction.copy()
    selected = sensitivity["affected_genes"]
    assert 0 < len(selected) <= 5
    for row in sensitivity["edit_rows"]:
        g = bundle.genes.index(row["gene"])
        shift = bundle.target.mean(0)[g] - bundle.prediction.mean(0)[g]
        cf[:, g] = np.maximum(cf[:, g] + shift, 0)
    expected_before = core.mmd_unbiased(bundle.prediction, bundle.target, seed=2)
    expected_after = core.mmd_unbiased(cf, bundle.target, seed=2)
    assert sensitivity["before"] == pytest.approx(expected_before, abs=1e-12)
    assert sensitivity["after"] == pytest.approx(expected_after, abs=1e-12)
    assert sensitivity["change"] == pytest.approx(expected_after - expected_before, abs=1e-12)
    assert all(row["after_mean_error"] == pytest.approx(0, abs=3e-6) for row in sensitivity["edit_rows"])
    assert np.isfinite(cf).all() and (cf >= 0).all()
    assert originals == tuple(X.tobytes() for X in (bundle.prediction, bundle.target, bundle.reference))
    assert all(not X.flags.writeable for X in (bundle.prediction, bundle.target, bundle.reference))


def test_collapsed_population_has_good_mean_but_zero_variance(upstream):
    bundle = make_bundle(mode="collapse")
    result = diagnose(bundle, score_fixture(bundle, upstream))
    assert result["gene_summary"]["pb_rel_err"] < 3e-6
    assert result["distribution"]["variance_ratio"] < 1e-8
    assert all(r["variance_pred"] < 1e-8 for r in result["genes"])


def test_nonnegative_clipping_exposes_incomplete_mean_alignment(upstream):
    bundle = make_bundle()
    P = bundle.prediction.copy()
    P[:, 30] = 0
    P[:8, 30] = 40
    P.setflags(write=False)
    bundle.prediction = P
    result = diagnose(bundle, score_fixture(bundle, upstream))
    row = next(r for r in result["distribution"]["mmd_sensitivity"]["edit_rows"] if r["gene"] == "g30")
    assert row["requested_shift"] < 0
    assert row["nonnegative_clipped_cells"] == 56
    assert abs(row["achieved_mean_shift"]) < abs(row["requested_shift"])
    assert row["after_mean_error"] > 0.5
