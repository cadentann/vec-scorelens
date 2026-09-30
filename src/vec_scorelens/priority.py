"""Conservative, explicit diagnostic follow-up conventions, never official scores.

These rules are heuristics on supplied selected cells, not calibrated tests,
false-discovery rates, biological-health labels, or model-repair instructions.
"""
from __future__ import annotations

import math
from collections.abc import Mapping
from copy import deepcopy
from typing import Any


DEFAULT_PROFILE = {
    "id": "scorelens-diagnostic-conventions-v1",
    "version": 1,
    "origin": "diagnostic_convention",
    "calibration": "heuristic; not statistically calibrated; no guaranteed false-alarm rate or FDR",
    "scope": "selected supplied cells and supplied gene panel; original counts are metadata only",
    "max_findings": 6,
    "rules": {
        "mean_residual": {"floor": 0.25, "se_multiplier": 4.0, "comparator": ">",
                          "units": "supplied transformed expression", "formula": "max(floor, se_multiplier*sqrt(variance_pred/n_pred + variance_truth/n_target))"},
        "opposite_sign": {"pred_delta_floor": 0.25, "comparator": ">",
                          "eligibility": "target operational DE call and literal opposite signs", "units": "supplied transformed expression"},
        "composition_gap": {"floor": 0.075, "se_multiplier": 4.0, "comparator": ">",
                            "units": "proportion", "formula": "max(floor, se_multiplier*sqrt(target_fraction*(1-target_fraction)*(1/n_pred + 1/n_target)))"},
        "variance_ratio": {"lower": 0.5, "upper": 2.0, "comparator": "outside_inclusive_range", "units": "ratio"},
        "response_magnitude": {"lower": 0.75, "upper": 1.25, "comparator": "outside_inclusive_range",
                               "eligibility": "target DE, reliable magnitude denominator, matching literal sign", "units": "median absolute response ratio"},
        "spatial_size": {"threshold": math.log(1.25), "comparator": ">", "units": "absolute selected log RMS-radius ratio"},
        "original_count": {"threshold": 0.25, "comparator": ">", "units": "absolute original prediction/target count ratio minus one"},
    },
    "ordering": ["opposite_sign", "mean_residual", "response_magnitude", "composition_gap", "variance_ratio", "spatial_size", "original_count"],
    "excluded": ["raw official metric ranking", "SW/Dice shape priority", "sampled covariance priority",
                 "local-neighborhood priority", "laterality", "biological health", "sensitivity-based repair advice"],
}


def _mapping(value: Any) -> Mapping:
    return value if isinstance(value, Mapping) else {}


def _number(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        converted = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return converted if math.isfinite(converted) else None


def _rows(value: Any) -> list[Mapping]:
    return [row for row in value if isinstance(row, Mapping)] if isinstance(value, list) else []


def _counts(diagnostics: Mapping) -> tuple[float | None, float | None]:
    scope = _mapping(diagnostics.get("scope"))
    counts = _mapping(scope.get("selected_counts"))
    pred, target = _number(counts.get("prediction")), _number(counts.get("target"))
    if pred is None or target is None or pred <= 0 or target <= 0:
        return None, None
    return pred, target


def _profile(value: Mapping | None) -> dict:
    profile = deepcopy(DEFAULT_PROFILE if value is None else dict(value))
    rules = _mapping(profile.get("rules"))
    if not all(rule in rules for rule in DEFAULT_PROFILE["rules"]):
        raise ValueError("priority profile must explicitly define all seven diagnostic conventions")
    for rule, spec in rules.items():
        if rule not in DEFAULT_PROFILE["rules"]:
            raise ValueError(f"unsupported priority rule: {rule}")
        for field in DEFAULT_PROFILE["rules"][rule]:
            if field not in spec:
                raise ValueError(f"priority profile rule {rule} is missing {field}")
            if isinstance(DEFAULT_PROFILE["rules"][rule][field], (int, float)):
                numeric = _number(spec[field])
                if numeric is None or numeric <= 0:
                    raise ValueError(f"priority profile {rule}.{field} must be finite and positive")
            elif spec[field] != DEFAULT_PROFILE["rules"][rule][field]:
                raise ValueError(f"priority profile {rule}.{field} must match the implemented convention")
        if "lower" in spec and float(spec["lower"]) >= float(spec["upper"]):
            raise ValueError(f"priority profile {rule} requires lower < upper")
    if profile.get("ordering") != DEFAULT_PROFILE["ordering"]:
        raise ValueError("priority profile must retain the documented deterministic domain ordering")
    limit = _number(profile.get("max_findings"))
    if limit is None or int(limit) != limit or not 1 <= int(limit) <= 6:
        raise ValueError("priority profile max_findings must be an integer from 1 to 6")
    return profile


def prioritize(diagnostics: Mapping, *, task: str = "T1", profile: Mapping | None = None) -> dict:
    """Return auditable follow-up findings using only defined diagnostic evidence.

    At most six findings are shown in fixed domain order. All trigger rows and
    threshold calculations are retained in ``evaluations`` independently of
    that display limit. Missing required evidence makes coverage incomplete.
    """
    if task not in ("T1", "T2", "T3"):
        raise ValueError("priority task must be T1, T2, or T3")
    configured = _profile(profile)
    rules = configured["rules"]
    n_pred, n_target = _counts(diagnostics)
    evaluations: list[dict] = []
    findings: list[dict] = []
    unavailable: list[dict] = []
    incomplete: list[str] = []

    def missing(rule: str, reason: str, *, required: bool = True) -> None:
        unavailable.append({"rule_id": rule, "reason": reason})
        if required:
            incomplete.append(rule)

    def add(rule: str, domain: str, title: str, rows: list[dict], check: str,
            *, evidence_class: str = "A", scope: str = "selected supplied cells") -> None:
        evaluations.append({"rule_id": rule, "status": "evaluated", "n_evaluated": len(rows),
                            "rows": rows, "definition": rules[rule]})
        triggered = [row for row in rows if row["triggered"]]
        if not triggered:
            return
        # Only comparable rows within a rule are ranked; never rank across domains.
        triggered.sort(key=lambda row: (-row["exceedance"], str(row["entity"])))
        findings.append({
            "id": f"{task.lower()}:{rule}", "rule_id": rule, "diagnostic_domain": domain,
            "title": title, "n_triggered": len(triggered), "evidence_class": evidence_class,
            "scope": scope, "evidence_ids": [row["evidence_id"] for row in triggered],
            "evidence": triggered[:5], "rationale": "Measured evidence exceeds this explicit diagnostic convention.",
            "suggested_check": check,
            "limitations": [configured["calibration"], "Follow-up order is a convention; it is not a comparison of unlike official metrics or a biological-health conclusion."],
        })

    genes = _rows(diagnostics.get("genes"))
    wrong: list[dict] = []
    residuals: list[dict] = []
    skipped_mean = 0
    skipped_sign = 0
    for row in genes:
        label = str(row.get("gene", "unnamed gene"))
        dp, dt = _number(row.get("change_pred")), _number(row.get("change_truth"))
        if dp is not None and dt is not None:
            if row.get("truth_de_direction") in ("up", "down"):
                floor = float(rules["opposite_sign"]["pred_delta_floor"])
                opposite = dp * dt < 0 and abs(dp) > floor
                wrong.append({"entity": label, "evidence_id": f"genes:{label}:literal_sign",
                              "measured_value": dp, "target_delta": dt, "threshold": floor,
                              "comparator": "opposite sign and absolute predicted delta >", "units": rules["opposite_sign"]["units"],
                              "triggered": opposite, "exceedance": abs(dp) / floor if opposite else 0.0})
        else:
            skipped_sign += 1
        error = _number(row.get("absolute_error"))
        vp, vt = _number(row.get("variance_pred")), _number(row.get("variance_truth"))
        if error is None or vp is None or vt is None or vp < 0 or vt < 0 or n_pred is None:
            skipped_mean += 1
            continue
        se = math.sqrt(vp / n_pred + vt / n_target)
        threshold = max(float(rules["mean_residual"]["floor"]), float(rules["mean_residual"]["se_multiplier"]) * se)
        residuals.append({"entity": label, "evidence_id": f"genes:{label}:absolute_error", "measured_value": error,
                          "threshold": threshold, "estimated_cell_se": se, "comparator": ">",
                          "units": rules["mean_residual"]["units"], "triggered": error > threshold,
                          "exceedance": error / threshold})
    if genes and not skipped_sign:
        add("opposite_sign", "response", "Meaningful literal response reversals", wrong,
            "Inspect named gene contrasts, input alignment, and normalization; operational target DE is not confirmed biological regulation.")
    else:
        missing("opposite_sign", "Gene contrasts are missing or undefined.")
    if residuals:
        add("mean_residual", "expression", "Large selected-sample mean residuals", residuals,
            "Inspect named gene means, input alignment, normalization, and the supplied target/reference contrast.")
    if not residuals or skipped_mean:
        missing("mean_residual", f"Counts/variances/residuals missing or undefined for {skipped_mean or len(genes)} gene row(s).")

    if task == "T3":
        response = _mapping(diagnostics.get("response"))
        response_rows = _rows(response.get("genes"))
        response_accepted = _mapping(response.get("summary")).get("status") == "accepted"
        eligible: list[tuple[str, float]] = []
        for row in response_rows:
            ratio = _number(row.get("magnitude_ratio"))
            dp, dt = _number(row.get("pred_delta")), _number(row.get("truth_delta"))
            # Source response module must mark the denominator reliable. Both
            # signs must agree; amplitude is not interpretable for a reversal.
            status = str(row.get("magnitude_status", ""))
            reliable = status in ("defined", "reliable", "eligible", "defined_reliable", "under", "over", "matched")
            if response_accepted and row.get("truth_de_eligible") is True and reliable and ratio is not None and ratio >= 0 and dp is not None and dt is not None and dp * dt > 0:
                eligible.append((str(row.get("gene", "unnamed gene")), ratio))
        if eligible:
            values = sorted(value for _, value in eligible)
            middle = len(values) // 2
            median = values[middle] if len(values) % 2 else (values[middle - 1] + values[middle]) / 2
            low, high = float(rules["response_magnitude"]["lower"]), float(rules["response_magnitude"]["upper"])
            add("response_magnitude", "response", "Eligible WT-relative response magnitude mismatch",
                [{"entity": "reliable same-sign target-DE genes", "evidence_id": "response:magnitude_ratio_median",
                  "measured_value": median, "threshold": [low, high], "comparator": "outside inclusive range",
                  "units": rules["response_magnitude"]["units"], "triggered": median < low or median > high,
                  "exceedance": low / median if 0 < median < low else median / high,
                  "n_eligible": len(eligible), "eligible_genes": sorted(label for label, _ in eligible)}],
                "Inspect WT matching and eligible response deltas; this ratio is separate from the guarded official severity slope.")
        else:
            missing("response_magnitude", "The accepted severity reliability gate and reliable same-sign target-DE ratios are required; no amplitude all-clear is possible.")

    composition: list[dict] = []
    comp_rows = _rows(diagnostics.get("composition"))
    skipped_comp = 0
    for row in comp_rows:
        target, pred = _number(row.get("target_fraction")), _number(row.get("pred_soft_fraction"))
        if n_pred is None or target is None or pred is None or not 0 <= target <= 1 or not 0 <= pred <= 1:
            skipped_comp += 1
            continue
        se = math.sqrt(target * (1 - target) * (1 / n_pred + 1 / n_target))
        threshold = max(float(rules["composition_gap"]["floor"]), float(rules["composition_gap"]["se_multiplier"]) * se)
        label = str(row.get("celltype", "unnamed class"))
        gap = abs(pred - target)
        composition.append({"entity": label, "evidence_id": f"composition:{label}:proportion_error",
                            "measured_value": gap, "threshold": threshold, "estimated_cell_se": se,
                            "comparator": ">", "units": "proportion", "triggered": gap > threshold, "exceedance": gap / threshold})
    if composition:
        add("composition_gap", "composition", "Frozen-probe composition discrepancies", composition,
            "Inspect target labels and the fitted probe domain. Low predicted soft proportions do not establish biological absence.")
    if not composition or skipped_comp:
        missing("composition_gap", "Target/probe fractions or selected counts are missing or undefined.")

    ratio = _number(_mapping(diagnostics.get("distribution")).get("variance_ratio"))
    if ratio is None or ratio < 0:
        missing("variance_ratio", "Selected global variance ratio is undefined or missing.")
    else:
        low, high = float(rules["variance_ratio"]["lower"]), float(rules["variance_ratio"]["upper"])
        add("variance_ratio", "distribution", "Selected variability discrepancy",
            [{"entity": "global selected variance", "evidence_id": "distribution:variance_ratio", "measured_value": ratio,
              "threshold": [low, high], "comparator": "outside inclusive range", "units": "ratio",
              "triggered": ratio < low or ratio > high, "exceedance": low / max(ratio, 1e-300) if ratio < low else ratio / high}],
            "Inspect selected cell heterogeneity, normalization, and sampling; this descriptive ratio does not identify a biological mechanism.", evidence_class="D")

    if task in ("T2", "T3"):
        spatial = _mapping(diagnostics.get("spatial"))
        growth = _mapping(spatial.get("growth"))
        scale, count = _mapping(growth.get("scale")), _mapping(growth.get("count"))
        value = _number(scale.get("replayed_value"))
        if value is None:
            missing("spatial_size", "Selected RMS-radius ratio is missing or undefined.")
        else:
            threshold = float(rules["spatial_size"]["threshold"])
            add("spatial_size", "spatial", "Selected coordinate size discrepancy",
                [{"entity": "selected RMS radius", "evidence_id": "spatial:growth:scale", "measured_value": abs(value),
                  "signed_value": value, "threshold": threshold, "comparator": ">", "units": rules["spatial_size"]["units"],
                  "triggered": abs(value) > threshold, "exceedance": abs(value) / threshold}],
                "Inspect coordinate units and selected-cloud size. This comparison does not establish biological growth.", evidence_class="D")
        np_orig, nt_orig = _number(count.get("original_prediction_count")), _number(count.get("original_target_count"))
        if np_orig is None or nt_orig is None or np_orig <= 0 or nt_orig <= 0:
            missing("original_count", "Original file-count metadata is missing or undefined.")
        else:
            value = abs(np_orig / nt_orig - 1)
            threshold = float(rules["original_count"]["threshold"])
            add("original_count", "spatial", "Original file-count discrepancy",
                [{"entity": "original file cell count", "evidence_id": "spatial:growth:original_count", "measured_value": value,
                  "threshold": threshold, "comparator": ">", "units": rules["original_count"]["units"],
                  "triggered": value > threshold, "exceedance": value / threshold,
                  "original_prediction_count": np_orig, "original_target_count": nt_orig}],
                "Inspect input cell counts and sampling. Original metadata counts do not establish proliferation or a full-file expression match.",
                evidence_class="D", scope="original file-count metadata only")
        missing("shape_and_local_priority", "Shape and local organization are measured separately but have no adopted priority thresholds; SW/Dice have a source-frame caveat.", required=False)
        missing("laterality", "No calibrated anatomical laterality diagnostic is computed.", required=False)
    missing("pair_priority", "Sampled variogram and descriptive covariance remain evidence tables without an adopted priority threshold.", required=False)
    missing("sensitivity_priority", "Target-aware sensitivity is noncausal and does not determine follow-up priority.", required=False)
    order = {name: index for index, name in enumerate(configured["ordering"])}
    findings.sort(key=lambda finding: (order[finding["rule_id"]], finding["id"]))
    for index, finding in enumerate(findings, 1):
        finding["priority_order"] = index
    all_findings = len(findings)
    findings = findings[:int(configured["max_findings"])]
    assessment = "thresholds_exceeded" if all_findings else "incomplete" if incomplete else "within_listed_thresholds"
    return {"schema_version": "vec-scorelens-priorities-v1", "profile": configured, "assessment": assessment,
            "coverage": "incomplete" if incomplete else "listed rules evaluated", "scope": "selected supplied cells and explicit metadata counts only",
            "findings": findings, "n_triggered_rules": all_findings, "n_findings_omitted": all_findings - len(findings),
            "evaluations": evaluations, "not_evaluated": unavailable, "incomplete_rules": list(dict.fromkeys(incomplete)),
            "interpretation": "Within listed thresholds means only that evaluated heuristic rules did not trigger. It is not a statement of biological health, full-file equivalence, or expected model improvement."}
