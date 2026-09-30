"""T3 matched-WT response evidence for the unchanged pinned severity metric.

Authoritative values come from upstream. The replay explains its ordered gates
and through-origin fit, while literal signs and magnitude categories are
descriptive. Nothing here allocates log severity or implies causal importance.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from .diagnostics import _finite_json


SIGN_ABS_TOL = 1e-6
SIGN_REL_TOL = 1e-5
MAGNITUDE_RATIO_TOL = 1e-3
RECONCILIATION_ATOL = 1e-5
RECONCILIATION_RTOL = 1e-5


def _agreement(left: float, right: float) -> tuple[bool, float | None]:
    if not np.isfinite(left) or not np.isfinite(right):
        return bool(np.isnan(left) and np.isnan(right)), None
    return bool(np.isclose(left, right, atol=RECONCILIATION_ATOL,
                           rtol=RECONCILIATION_RTOL)), abs(float(left - right))


def _literal(value: float, tolerance: float = 0.0) -> str:
    return "up" if value > tolerance else "down" if value < -tolerance else "zero"


def diagnose_response(bundle: Any, scored: dict) -> dict:
    """Explain T3 response against the supplied WT on the selected arrays.

    The official module is supplied by the verified scorer. Coordinates have no
    role in this expression-only fit; their validation belongs to the adapter.
    """
    if getattr(bundle, "task", "T3") != "T3":
        raise ValueError("diagnose_response requires a T3 bundle with matched WT reference")
    core = scored["core"]
    P, T, W = bundle.prediction, bundle.target, bundle.reference
    mu_p, mu_t, mu_w = (core.pseudobulk(X) for X in (P, T, W))
    dp_all, dt_all = mu_p - mu_w, mu_t - mu_w
    up_t, dn_t, _ = core.de_genes(T, W)
    up_p, dn_p, _ = core.de_genes(P, W)
    idx = np.concatenate([up_t, dn_t]).astype(int)
    official_logbeta, official_r2 = core.severity_slope(P, T, W)
    sentinel = float(core.SEVERITY_LOG_WORST)
    std_p, std_t = float(np.std(dp_all)), float(np.std(dt_all))
    no_response_threshold = .01 * std_t
    # Diagnostics have their own explicit tolerance; upstream guards never use it.
    sign_tolerance = max(SIGN_ABS_TOL, SIGN_REL_TOL * std_t)

    gate_names = ["too_few_truth_de", "no_truth_response", "no_predicted_response",
                  "degenerate_truth_de_response", "nonpositive_beta_or_low_r2"]
    gates = [{"condition": name, "evaluated": False, "triggered": None}
             for name in gate_names]

    def gate(number: int, condition: bool) -> bool:
        gates[number].update(evaluated=True, triggered=bool(condition))
        return bool(condition)

    nan = float("nan")
    ss = numerator = beta = sse = centered_ss = r2_denominator = r2 = nan
    logbeta_candidate = replay_logbeta = replay_r2 = nan
    fit_reached = False
    residual = None
    if gate(0, len(idx) < 5):
        status = "undefined_too_few_truth_de"
    elif gate(1, std_t < 1e-12):
        status = "undefined_no_truth_response"
    elif gate(2, std_p < no_response_threshold):
        status = "no_predicted_response"
        replay_logbeta, replay_r2 = sentinel, 0.0
    else:
        dp, dt = dp_all[idx], dt_all[idx]
        ss = float(dt @ dt)
        if gate(3, ss < 1e-12):
            status = "undefined_degenerate_truth_de_response"
        else:
            fit_reached = True
            numerator = float(dt @ dp)
            beta = numerator / ss
            # Retain the upstream dtype/operation order before scalar conversion.
            residual = dp - beta * dt
            sse = float(residual @ residual)
            centered_ss = float(((dp - dp.mean()) ** 2).sum())
            r2_denominator = max(centered_ss, 1e-12)
            r2 = 1.0 - sse / r2_denominator
            replay_r2 = round(r2, 4)
            if beta > 0:
                logbeta_candidate = float(np.log(beta))
            if gate(4, beta <= 0 or r2 <= .01):
                if beta <= 0 and r2 <= .01:
                    status = "nonpositive_beta_and_low_r2"
                elif beta <= 0:
                    status = "nonpositive_beta"
                else:
                    status = "low_r2"
                replay_logbeta = sentinel
            else:
                status = "accepted"
                replay_logbeta = logbeta_candidate

    origin = ("accepted_log_beta" if status == "accepted" else
              "undefined_board" if status.startswith("undefined_") else
              "failure_sentinel")
    log_match, log_error = _agreement(replay_logbeta, float(official_logbeta))
    r2_match, r2_error = _agreement(replay_r2, float(official_r2))
    metrics = scored.get("metrics", {})
    rounded_value = metrics.get("severity_slope")
    if "severity_slope" in metrics:
        rounded_match = (rounded_value is None and not np.isfinite(official_logbeta)) or (
            rounded_value is not None and np.isfinite(official_logbeta)
            and float(rounded_value) == round(float(official_logbeta), 4))
        rounded_error = (abs(float(rounded_value) - float(official_logbeta))
                         if rounded_value is not None and np.isfinite(official_logbeta) else None)
    else:
        rounded_match, rounded_error = None, None
    if "_slope_r2" in metrics:
        reported_r2 = metrics["_slope_r2"]
        stored_r2 = nan if reported_r2 is None else float(reported_r2)
        stored_r2_match, stored_r2_error = _agreement(stored_r2, float(official_r2))
    else:
        stored_r2_match, stored_r2_error = None, None

    truth_de = np.full(len(bundle.genes), "none", dtype="<U4")
    local_de = truth_de.copy()
    truth_de[up_t], truth_de[dn_t] = "up", "down"
    local_de[up_p], local_de[dn_p] = "up", "down"
    fit_position = {int(g): k for k, g in enumerate(idx)}
    direct_error = mu_p - mu_t
    rows = []
    for g, gene in enumerate(bundle.genes):
        truth_sign = _literal(float(dt_all[g]), sign_tolerance)
        pred_sign = _literal(float(dp_all[g]), sign_tolerance)
        if truth_sign == "zero":
            sign_status = "truth_zero"
        elif pred_sign == "zero":
            sign_status = "pred_zero"
        elif truth_sign == pred_sign:
            sign_status = "agree"
        else:
            sign_status = "opposite"
        eligible = truth_de[g] != "none"
        ratio = None
        if not eligible:
            magnitude_status = "not_interpretable_non_de"
        elif sign_status != "agree":
            magnitude_status = "not_interpretable_" + sign_status
        else:
            ratio = abs(float(dp_all[g]) / float(dt_all[g]))
            magnitude_status = ("matched" if abs(ratio - 1.0) <= MAGNITUDE_RATIO_TOL
                                else "under" if ratio < 1.0 else "over")
        row = {
            "gene": gene, "mean_pred": mu_p[g], "mean_truth": mu_t[g], "mean_wt": mu_w[g],
            "pred_delta": dp_all[g], "truth_delta": dt_all[g], "direct_error": direct_error[g],
            "absolute_error": abs(float(direct_error[g])),
            "truth_de_direction": truth_de[g], "truth_de_eligible": eligible,
            "pred_local_de_direction": local_de[g],
            "local_de_hit": bool(eligible and local_de[g] == truth_de[g]),
            "local_de_miss": bool(eligible and local_de[g] != truth_de[g]),
            "local_de_false_positive": bool(local_de[g] != "none" and local_de[g] != truth_de[g]),
            "literal_pred_direction": _literal(float(dp_all[g])),
            "literal_truth_direction": _literal(float(dt_all[g])),
            "tolerance_pred_direction": pred_sign, "tolerance_truth_direction": truth_sign,
            "literal_sign_status": sign_status, "sign_tolerance": sign_tolerance,
            "magnitude_ratio": ratio, "magnitude_status": magnitude_status,
            "slope_eligible": eligible, "slope_fit_reached": bool(eligible and fit_reached),
            "truth_delta_squared": float(dt_all[g]) ** 2 if eligible else None,
            "slope_numerator_term": None, "beta_allocation": None,
            "fitted_pred_delta": None, "slope_residual": None,
            "slope_residual_squared": None, "r2_deficit_allocation": None,
        }
        if fit_reached and eligible:
            k = fit_position[g]
            term = float(dt_all[g]) * float(dp_all[g])
            residual_g = float(residual[k])
            row.update(slope_numerator_term=term, beta_allocation=term / ss,
                       fitted_pred_delta=float(beta * dt_all[g]), slope_residual=residual_g,
                       slope_residual_squared=residual_g ** 2,
                       r2_deficit_allocation=residual_g ** 2 / r2_denominator)
        rows.append(row)

    component_reconciliation = {}
    if fit_reached:
        for key, expected in [("slope_numerator_term", numerator), ("beta_allocation", beta),
                              ("slope_residual_squared", sse), ("r2_deficit_allocation", 1 - r2)]:
            actual = sum(row[key] for row in rows if row[key] is not None)
            matches, error = _agreement(actual, expected)
            component_reconciliation[key] = {"sum": actual, "expected": expected,
                                              "abs_error": error, "matches": matches}
    magnitude_ratios = [row["magnitude_ratio"] for row in rows if row["magnitude_ratio"] is not None]
    summary = {
        "definition": "Predicted and observed mean log-expression shifts relative to supplied matched WT; ordered official severity gates and through-origin response fit.",
        "evidence_type": "algebraic_fit_components_and_descriptive_response_mismatch",
        "limitations": [
            "Mean log-expression shifts are not fold changes of aggregated raw counts.",
            "WT matching is user-supplied and is not independently verified by expression similarity.",
            "The fit is over truth DE genes; its no-response guard uses standard deviation across all genes, so even a large uniform nonzero offset can trigger it. This source branch is not proof of absent predicted change.",
            "R2 uses centered predicted response in the denominator and can be negative; an inverted response can have R2=1 but still fail beta.",
            "The failure sentinel is not an estimated beta or a lower bound on accepted log beta.",
            "Fit components reconstruct beta and residual fit, not additive log-severity contributions or causal gene importance.",
            "Local prediction-vs-WT significant calls and tolerance-based signs do not define official DE rank slots or partial rank correlation.",
            "Magnitude categories are descriptive and require truth DE eligibility and matching nonzero tolerance-based signs.",
            "The sign-reversal count includes all genes with opposite tolerance-based signs, including non-DE genes with tiny changes; it is not a count of significant or biologically meaningful response reversals.",
            "Cell-level significance is conditional on supplied cells, not independent embryo replication.",
            "All evidence applies to the selected supplied arrays; coordinates are evaluated separately and no cross-embryo cell pairing is assumed.",
        ],
        "reference_role": "matched_wild_type", "matching_status": "user_supplied_not_independently_verified",
        "status": status, "ordered_gates": gates, "fit_reached": fit_reached,
        "n_truth_de": len(idx), "n_truth_up": len(up_t), "n_truth_down": len(dn_t),
        "std_pred_all": std_p, "std_truth_all": std_t,
        "relative_no_response_threshold": no_response_threshold,
        "de_truth_ss": ss, "beta_numerator": numerator, "beta": beta,
        "residual_ss": sse, "pred_centered_ss": centered_ss, "r2_denominator": r2_denominator,
        "r2_unrounded": r2, "official_r2": official_r2, "official_logbeta": official_logbeta,
        "failure_sentinel": sentinel, "logbeta_candidate": logbeta_candidate,
        "reported_logbeta_origin": origin, "replayed_logbeta": replay_logbeta,
        "replayed_r2": replay_r2, "reconciliation_matches": log_match and r2_match,
        "reconciliation_abs_error": log_error, "r2_reconciliation_abs_error": r2_error,
        "rounded_official_reconciliation_matches": rounded_match,
        "rounded_official_reconciliation_abs_error": rounded_error,
        "stored_r2_reconciliation_matches": stored_r2_match,
        "stored_r2_reconciliation_abs_error": stored_r2_error,
        "component_reconciliation": component_reconciliation,
        "reconciliation_tolerance": {"atol": RECONCILIATION_ATOL, "rtol": RECONCILIATION_RTOL},
        "sign_tolerance": sign_tolerance,
        "sign_tolerance_definition": "max(1e-6, 1e-5 * std(truth delta)); descriptive only, never used in official gates",
        "magnitude_ratio_tolerance": MAGNITUDE_RATIO_TOL,
        "magnitude_status": "evaluated" if magnitude_ratios else "not_evaluated",
        "n_magnitude_eligible": len(magnitude_ratios),
        "magnitude_median_ratio": float(np.median(magnitude_ratios)) if magnitude_ratios else None,
        "n_local_hits": sum(row["local_de_hit"] for row in rows),
        "n_local_misses": sum(row["local_de_miss"] for row in rows),
        "n_local_false_positives": sum(row["local_de_false_positive"] for row in rows),
        "n_sign_reversals": sum(row["literal_sign_status"] == "opposite" for row in rows),
    }
    warnings = []
    if status == "no_predicted_response":
        warnings.append("Source no_predicted_response means low variation of the predicted change across genes, not necessarily zero change: a large uniform nonzero offset can also trigger this guard.")
    if status != "accepted":
        warnings.append(f"Severity explanation branch: {status}; reported value origin: {origin}.")
    if not log_match or not r2_match or rounded_match is False or stored_r2_match is False:
        warnings.append("Severity explanation does not reconcile with the unchanged official value; inspect source, dtype and task context.")
    if any(not item["matches"] for item in component_reconciliation.values()):
        warnings.append("Severity component allocation exceeds the declared float32 reconciliation tolerance.")
    return _finite_json({"genes": rows, "summary": summary, "warnings": warnings})
