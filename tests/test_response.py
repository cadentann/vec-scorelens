"""Source-grounded T3 severity gates and descriptive matched-WT response evidence."""
from types import SimpleNamespace

import json
import numpy as np
import pytest

from vec_scorelens.response import diagnose_response


@pytest.fixture(scope="module")
def core():
    from common import core_metrics
    return core_metrics


def make_bundle(mode="correct", n_de=12):
    rng = np.random.default_rng(9030)
    n, G = 64, 40
    W = 20 + rng.normal(0, .025, (n, G))
    delta = np.zeros(G)
    delta[:6] = np.linspace(.4, .9, 6)
    delta[6:12] = -np.linspace(.4, .9, 6)
    if n_de != 12:
        delta[:] = 0
        delta[:n_de] = np.linspace(.4, .9, n_de)
    prediction_delta = delta.copy()
    if mode == "identity":
        prediction_delta[:] = 0
    elif mode == "half":
        prediction_delta *= .5
    elif mode == "double":
        prediction_delta *= 2
    elif mode == "inverted":
        prediction_delta *= -1
    elif mode == "unrelated":
        v = np.tile([1., -1.], 6)
        v -= v.mean()
        v -= delta[:12] * (v @ delta[:12]) / (delta[:12] @ delta[:12])
        v /= np.linalg.norm(v)
        prediction_delta[:12] = .1 * delta[:12] + 5 * v
    elif mode == "tinybeta":
        prediction_delta[:12] *= .0001
        prediction_delta[12:] = np.linspace(-2, 2, G - 12)
    elif mode == "below_guard":
        prediction_delta *= .009
    elif mode == "above_guard":
        prediction_delta *= .011
    elif mode == "local_calls":
        prediction_delta[0] = 0
        prediction_delta[6] *= -1
        prediction_delta[20] = .7
    elif mode == "no_truth_response":
        W = np.broadcast_to(np.linspace(4, 4.1, n)[:, None], (n, G)).copy()
        delta[:] = .5
        prediction_delta = delta.copy()
    elif mode == "constant_target":
        W = np.full((n, G), 4.)
        delta[:] = 0
        prediction_delta[:] = 0
    elif mode == "constant_shift":
        # Exact repeated columns preserve constant pseudobulk shift in float32.
        W = np.broadcast_to(np.linspace(4, 4.1, n)[:, None], (n, G)).copy()
        prediction_delta[:] = 1
    P, T = W + prediction_delta, W + delta
    P, T, W = (np.asarray(X, np.float32) for X in (P, T, W))
    for X in (P, T, W):
        X.setflags(write=False)
    return SimpleNamespace(task="T3", prediction=P, target=T, reference=W,
                           genes=[f"g{g:02}" for g in range(G)])


def explain(bundle, core):
    official, r2 = core.severity_slope(bundle.prediction, bundle.target, bundle.reference)
    scored = {"core": core, "metrics": {
        "severity_slope": round(official, 4) if np.isfinite(official) else None,
        "_slope_r2": r2,
    }}
    return diagnose_response(bundle, scored)


@pytest.mark.parametrize("mode,status", [
    ("correct", "accepted"), ("half", "accepted"), ("double", "accepted"),
    ("identity", "no_predicted_response"), ("inverted", "nonpositive_beta"),
    ("unrelated", "low_r2"), ("below_guard", "no_predicted_response"),
    ("above_guard", "accepted"), ("tinybeta", "accepted"),
    ("no_truth_response", "undefined_no_truth_response"),
    ("constant_target", "undefined_too_few_truth_de"),
    ("constant_shift", "no_predicted_response"),
])
def test_ordered_severity_branches_reconcile_unchanged_core(core, mode, status):
    bundle = make_bundle(mode)
    originals = [X.tobytes() for X in (bundle.prediction, bundle.target, bundle.reference)]
    result = explain(bundle, core)
    summary = result["summary"]
    assert summary["status"] == status
    assert summary["reconciliation_matches"] is True
    assert summary["rounded_official_reconciliation_matches"] is True
    assert summary["stored_r2_reconciliation_matches"] is True
    assert originals == [X.tobytes() for X in (bundle.prediction, bundle.target, bundle.reference)]
    assert all(not X.flags.writeable for X in (bundle.prediction, bundle.target, bundle.reference))
    json.dumps(result, allow_nan=False)
    evaluated = [gate["evaluated"] for gate in summary["ordered_gates"]]
    assert evaluated == sorted(evaluated, reverse=True)
    for component in summary["component_reconciliation"].values():
        assert component["matches"] is True


@pytest.mark.parametrize("n_de", [0, 1, 4])
def test_less_than_five_truth_de_is_undefined_before_prediction_guard(core, n_de):
    result = explain(make_bundle("identity", n_de), core)
    summary = result["summary"]
    assert summary["status"] == "undefined_too_few_truth_de"
    assert summary["n_truth_de"] == n_de
    assert summary["official_logbeta"] is None
    assert summary["official_r2"] is None
    assert summary["beta"] is None
    assert summary["reported_logbeta_origin"] == "undefined_board"
    assert all(row["slope_residual"] is None for row in result["genes"])


def test_five_truth_de_permits_fit(core):
    summary = explain(make_bundle(n_de=5), core)["summary"]
    assert summary["status"] == "accepted"
    assert summary["n_truth_de"] == 5
    assert summary["beta"] == 1


@pytest.mark.parametrize("mode,beta,logbeta", [
    ("correct", 1, 0), ("half", .5, np.log(.5)), ("double", 2, np.log(2)),
])
def test_interpretable_severity_and_magnitude_scale(core, mode, beta, logbeta):
    result = explain(make_bundle(mode), core)
    summary = result["summary"]
    assert summary["beta"] == pytest.approx(beta, abs=3e-5)
    assert summary["official_logbeta"] == pytest.approx(logbeta, abs=3e-5)
    assert summary["official_r2"] == 1
    assert summary["magnitude_median_ratio"] == pytest.approx(beta, abs=3e-5)
    assert summary["n_magnitude_eligible"] == 12
    assert sum(row["beta_allocation"] or 0 for row in result["genes"]) == pytest.approx(summary["beta"], abs=1e-5)
    assert not any("logbeta_allocation" in row for row in result["genes"])


def test_inverted_fit_can_have_perfect_r2_but_fails_beta(core):
    summary = explain(make_bundle("inverted"), core)["summary"]
    assert summary["beta"] == pytest.approx(-1, abs=3e-5)
    assert summary["official_r2"] == 1
    assert summary["official_logbeta"] == summary["failure_sentinel"]
    assert summary["reported_logbeta_origin"] == "failure_sentinel"
    assert summary["n_sign_reversals"] == 12
    assert summary["magnitude_status"] == "not_evaluated"


def test_tiny_accepted_beta_is_not_clipped_to_failure_sentinel(core):
    summary = explain(make_bundle("tinybeta"), core)["summary"]
    assert 0 < summary["beta"] < .001
    assert summary["official_logbeta"] < summary["failure_sentinel"]
    assert summary["reported_logbeta_origin"] == "accepted_log_beta"
    assert summary["r2_unrounded"] > .01


def test_local_calls_signs_and_undefined_ratios_are_descriptive(core):
    result = explain(make_bundle("local_calls"), core)
    rows = {row["gene"]: row for row in result["genes"]}
    assert rows["g00"]["local_de_miss"] is True
    assert rows["g00"]["literal_sign_status"] == "pred_zero"
    assert rows["g00"]["magnitude_ratio"] is None
    assert rows["g06"]["local_de_miss"] is True
    assert rows["g06"]["local_de_false_positive"] is True
    assert rows["g06"]["literal_sign_status"] == "opposite"
    assert rows["g06"]["magnitude_ratio"] is None
    assert rows["g20"]["truth_de_eligible"] is False
    assert rows["g20"]["local_de_false_positive"] is True
    assert rows["g20"]["magnitude_ratio"] is None
    assert rows["g20"]["slope_residual"] is None


def test_no_response_fit_fields_are_null_not_estimated_beta_point001(core):
    result = explain(make_bundle("identity"), core)
    summary = result["summary"]
    assert summary["official_logbeta"] == np.log(.001)
    assert summary["official_r2"] == 0
    assert summary["beta"] is None
    assert summary["r2_unrounded"] is None
    assert summary["fit_reached"] is False
    assert summary["component_reconciliation"] == {}
    assert all(row["beta_allocation"] is None for row in result["genes"])


def test_coordinates_and_predicted_labels_have_no_role_in_response_fit(core):
    bundle = make_bundle()
    before = explain(bundle, core)
    bundle.prediction_coords = np.zeros((64, 3))
    bundle.target_coords = np.ones((64, 3))
    bundle.reference_coords = np.full((64, 3), 100.)
    bundle.prediction_labels = np.array(["untrusted"] * 64)
    assert explain(bundle, core) == before


def test_response_row_permutations_preserve_results_with_float_tolerance(core):
    bundle = make_bundle("double")
    original = explain(bundle, core)
    shuffled = SimpleNamespace(**vars(bundle))
    rng = np.random.default_rng(118)
    for name in ("prediction", "target", "reference"):
        array = getattr(shuffled, name)
        setattr(shuffled, name, array[rng.permutation(len(array))])
    summary = explain(shuffled, core)["summary"]
    assert summary["status"] == original["summary"]["status"]
    assert summary["beta"] == pytest.approx(original["summary"]["beta"], abs=1e-5)


def test_official_metric_context_mismatch_is_visible(core):
    result = diagnose_response(make_bundle(), {"core": core, "metrics": {"severity_slope": 3.}})
    assert result["summary"]["rounded_official_reconciliation_matches"] is False
    assert any("does not reconcile" in warning for warning in result["warnings"])


def test_official_r2_context_mismatch_is_visible(core):
    result = diagnose_response(make_bundle(), {"core": core, "metrics": {"_slope_r2": -.3}})
    assert result["summary"]["stored_r2_reconciliation_matches"] is False
    assert any("does not reconcile" in warning for warning in result["warnings"])


def test_non_t3_context_is_rejected(core):
    bundle = make_bundle()
    bundle.task = "T1"
    with pytest.raises(ValueError, match="T3"):
        explain(bundle, core)
