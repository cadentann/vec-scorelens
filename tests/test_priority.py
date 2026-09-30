"""Priority conventions: explicit numeric boundaries and incomplete coverage."""
from copy import deepcopy
import math

import pytest

from vec_scorelens.priority import DEFAULT_PROFILE, prioritize


def evidence():
    return {"scope": {"selected_counts": {"prediction": 256, "target": 256}},
            "genes": [{"gene": "g", "change_pred": .5, "change_truth": .5, "absolute_error": 0,
                       "variance_pred": 1, "variance_truth": 1, "truth_de_direction": "up"}],
            "composition": [{"celltype": "a", "target_fraction": .5, "pred_soft_fraction": .5}],
            "distribution": {"variance_ratio": 1}}


def triggered(result):
    return [finding["rule_id"] for finding in result["findings"]]


def test_matching_selected_evidence_within_listed_thresholds_and_exclusions():
    result = prioritize(evidence())
    assert result["assessment"] == "within_listed_thresholds"
    assert result["coverage"] == "listed rules evaluated"
    assert result["incomplete_rules"] == []
    assert {row["rule_id"] for row in result["not_evaluated"]} == {"pair_priority", "sensitivity_priority"}
    assert "not statistically calibrated" in result["profile"]["calibration"]
    assert result["profile"]["rules"]["mean_residual"]["floor"] == .25
    assert result["profile"]["rules"]["composition_gap"]["floor"] == .075


def test_mean_residual_floor_and_cell_noise_guard_are_visible():
    data = evidence()
    data["genes"][0]["absolute_error"] = .3
    result = prioritize(data)
    row = next(item for item in result["evaluations"] if item["rule_id"] == "mean_residual")["rows"][0]
    assert row["threshold"] == pytest.approx(math.sqrt(1 / 8))
    assert "mean_residual" not in triggered(result)
    data["genes"][0]["absolute_error"] = .4
    assert "mean_residual" in triggered(prioritize(data))
    data["genes"][0].update(variance_pred=0, variance_truth=0, absolute_error=.25)
    assert "mean_residual" not in triggered(prioritize(data))
    data["genes"][0]["absolute_error"] = .25001
    assert "mean_residual" in triggered(prioritize(data))


@pytest.mark.parametrize("delta,eligible,expected", [(-.250001, "up", True), (-.25, "up", False), (-.8, "none", False), (.8, "up", False)])
def test_literal_reversal_requires_target_de_and_meaningful_prediction(delta, eligible, expected):
    data = evidence()
    data["genes"][0].update(change_pred=delta, truth_de_direction=eligible)
    assert ("opposite_sign" in triggered(prioritize(data))) is expected


def test_composition_rule_exposes_counts_and_conservative_threshold():
    data = evidence()
    data["composition"][0]["pred_soft_fraction"] = .7
    result = prioritize(data)
    row = next(item for item in result["evaluations"] if item["rule_id"] == "composition_gap")["rows"][0]
    assert row["threshold"] == pytest.approx(math.sqrt(1 / 32))
    assert "composition_gap" in triggered(result)
    data["scope"]["selected_counts"] = {"prediction": 32, "target": 32}
    assert "composition_gap" not in triggered(prioritize(data))


@pytest.mark.parametrize("ratio,expected", [(0, True), (.49, True), (.5, False), (2, False), (2.01, True)])
def test_variance_bounds_are_inclusive(ratio, expected):
    data = evidence()
    data["distribution"]["variance_ratio"] = ratio
    result = prioritize(data)
    assert ("variance_ratio" in triggered(result)) is expected
    if expected:
        assert result["findings"][0]["evidence_class"] == "D"


def test_missing_counts_or_null_variance_prevents_all_clear():
    data = evidence()
    data["scope"] = {}
    data["distribution"]["variance_ratio"] = None
    result = prioritize(data)
    assert result["assessment"] == "incomplete"
    assert {"mean_residual", "composition_gap", "variance_ratio"} <= set(result["incomplete_rules"])


def test_sensitivity_pair_and_raw_metric_values_do_not_trigger():
    data = evidence()
    data["official_metrics"] = {"mmd_u": 1e99, "de_score": -1e99}
    data["pairs"] = [{"variogram_squared_error": 1e99, "covariance_error": 1e99}]
    data["distribution"]["mmd_sensitivity"] = {"change": -1e99}
    assert triggered(prioritize(data)) == []


def spatial_evidence():
    data = evidence()
    data["spatial"] = {"growth": {"scale": {"replayed_value": 0},
        "count": {"original_prediction_count": 1000, "original_target_count": 1000}}}
    return data


def test_spatial_size_and_original_count_do_not_rewrite_official_selected_count():
    data = spatial_evidence()
    data["official_metrics"] = {"count_log_ratio": 0}
    data["spatial"]["growth"]["scale"]["replayed_value"] = math.log(1.25001)
    data["spatial"]["growth"]["count"]["original_prediction_count"] = 1251
    before = deepcopy(data)
    result = prioritize(data, task="T2")
    assert triggered(result) == ["spatial_size", "original_count"]
    assert data == before
    assert result["findings"][1]["scope"] == "original file-count metadata only"
    assert {"shape_and_local_priority", "laterality"} <= {row["rule_id"] for row in result["not_evaluated"]}


@pytest.mark.parametrize("status,dp,ratio,expected", [("accepted", .9, 1.8, True), ("accepted", .5, 1, False), ("accepted", -.9, 1.8, False), ("no_predicted_response", .9, 1.8, False)])
def test_response_magnitude_requires_accepted_gate_and_matching_sign(status, dp, ratio, expected):
    data = spatial_evidence()
    data["response"] = {"summary": {"status": status}, "genes": [{"gene": "g", "truth_de_eligible": True,
        "truth_delta": .5, "pred_delta": dp, "magnitude_ratio": ratio, "magnitude_status": "over"}]}
    result = prioritize(data, task="T3")
    assert ("response_magnitude" in triggered(result)) is expected
    if status != "accepted" or dp < 0:
        assert "response_magnitude" in result["incomplete_rules"]


def test_domain_order_and_six_finding_cap_retain_all_evaluations():
    data = spatial_evidence()
    data["genes"][0].update(change_pred=-2, absolute_error=2)
    data["genes"].append({**data["genes"][0], "gene": "matching", "change_pred": 1})
    data["composition"][0]["pred_soft_fraction"] = 1
    data["distribution"]["variance_ratio"] = 3
    data["response"] = {"summary": {"status": "accepted"}, "genes": [{"gene": "matching", "truth_de_eligible": True,
        "truth_delta": .5, "pred_delta": 1, "magnitude_ratio": 2, "magnitude_status": "over"}]}
    data["spatial"]["growth"]["scale"]["replayed_value"] = 1
    data["spatial"]["growth"]["count"]["original_prediction_count"] = 2000
    result = prioritize(data, task="T3")
    assert triggered(result) == ["opposite_sign", "mean_residual", "response_magnitude", "composition_gap", "variance_ratio", "spatial_size"]
    assert result["n_triggered_rules"] == 7 and result["n_findings_omitted"] == 1
    assert any(item["rule_id"] == "original_count" for item in result["evaluations"])
    assert prioritize(data, task="T3") == result


def test_bad_profiles_rejected_and_default_profile_not_mutated():
    profile = deepcopy(DEFAULT_PROFILE)
    profile["rules"]["variance_ratio"]["lower"] = 3
    with pytest.raises(ValueError, match="lower < upper"):
        prioritize(evidence(), profile=profile)
    result = prioritize(evidence())
    result["profile"]["rules"]["variance_ratio"]["lower"] = 999
    assert DEFAULT_PROFILE["rules"]["variance_ratio"]["lower"] == .5
