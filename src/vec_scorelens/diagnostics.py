"""T1 explanations for the bounded, pinned upstream evaluation context.

Authoritative metric values always come from the supplied upstream modules.
Gene allocations, pair discrepancies, and matrix-edit sensitivity have distinct
meanings; none identifies a causal gene or a model-internal training defect.
"""

from __future__ import annotations

from typing import Any

import numpy as np


PAIR_BLOCK = 128
VARIOGRAM_PAIRS = 20_000
VARIOGRAM_CELL_CAP = 1_500
MMD_EDIT_GENE_CAP = 5


def _finite_json(value: Any) -> Any:
    """Convert numerical scalars to strict JSON values, preserving null status."""
    if isinstance(value, dict):
        return {str(k): _finite_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_finite_json(v) for v in value]
    if isinstance(value, np.ndarray):
        return _finite_json(value.tolist())
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    return value


def _directions(values: np.ndarray) -> np.ndarray:
    return np.where(values > 0, "up", np.where(values < 0, "down", "zero"))


def _slots(values: np.ndarray, n_up: int, n_down: int) -> tuple[np.ndarray, np.ndarray]:
    # The upstream scorer sorts float64 conversions and does not enforce signs.
    order = np.argsort(-np.asarray(values, float))
    slots = np.full(len(values), "none", dtype="<U4")
    slots[order[:n_up]] = "up"
    if n_down:
        slots[order[len(order) - n_down :]] = "down"
    ranks = np.empty(len(order), dtype=int)
    ranks[order] = np.arange(1, len(order) + 1)
    return slots, ranks


def _change_ledger(bundle: Any, core: Any) -> tuple[list[dict], dict, list[str]]:
    P, T, R = bundle.prediction, bundle.target, bundle.reference
    mu_p, mu_t, mu_r = (core.pseudobulk(X) for X in (P, T, R))
    up_t, dn_t, dt = core.de_genes(T, R)
    up_p, dn_p, dp = core.de_genes(P, R)
    de = core.de_score(P, T, R)
    direction = core.de_direction(P, T, R)
    # Direct cancellation matches pb_rel_err; subtracting two reference-relative
    # float32 differences need not produce the same last bit.
    error = mu_p - mu_t
    var_p, var_t = P.var(0), T.var(0)
    true_de = np.full(len(bundle.genes), "none", dtype="<U4")
    pred_de = true_de.copy()
    true_de[up_t], true_de[dn_t] = "up", "down"
    pred_de[up_p], pred_de[dn_p] = "up", "down"
    literal = _directions(dp)
    true_literal = _directions(dt)
    slots, ranks = _slots(dp, len(up_t), len(dn_t))
    hit = (true_de != "none") & (slots == true_de)
    local_hit = (true_de != "none") & (pred_de == true_de)
    local_miss = (true_de != "none") & ~local_hit
    local_fp = (pred_de != "none") & (pred_de != true_de)
    K = len(up_t) + len(dn_t)
    warnings: list[str] = []

    if not K:
        guard = "undefined_no_truth_de"
    elif np.std(dt) < 1e-12:
        guard = "undefined_degenerate_truth_change"
    elif np.std(dp) < 1e-2 * np.std(dt):
        guard = "no_change_std_guard"
    elif not np.isfinite(de["score"]):
        guard = "undefined_saturated_expression_null"
    else:
        guard = "none"
    if guard != "none":
        warnings.append(f"DE explanation branch: {guard}; rank-slot counts are descriptive in this branch.")

    if np.std(dt) < 1e-12:
        direction_guard = "degenerate_truth_change_returns_zero"
    elif np.std(dp) < 1e-2 * np.std(dt):
        direction_guard = "no_change_std_guard_returns_zero"
    else:
        rank_correlation = np.corrcoef(np.vstack([
            core._ranks(dp), core._ranks(dt), core._ranks(mu_r),
        ]))
        if not np.isfinite(rank_correlation).all():
            direction_guard = "nonfinite_rank_correlation_returns_zero"
        elif 1.0 - abs(rank_correlation[0, 2]) < 1e-6:
            direction_guard = "reference_rank_collinearity_returns_zero"
        else:
            direction_guard = "none"

    raw_allocation = np.zeros(len(bundle.genes), dtype=float)
    de_allocation = np.full(len(bundle.genes), np.nan)
    null_orientation: str | None = None
    c = float(de["chance"])
    if K:
        plus, _ = core._signed_overlap(mu_r, up_t, dn_t)
        minus, _ = core._signed_overlap(-mu_r, up_t, dn_t)
        null_orientation = "+reference_mean" if plus >= minus else "-reference_mean"
        null_slots, _ = _slots(mu_r if plus >= minus else -mu_r, len(up_t), len(dn_t))
        null_hit = (true_de != "none") & (null_slots == true_de)
        if np.isfinite(de["raw"]) and guard != "no_change_std_guard":
            raw_allocation = hit.astype(float) / K
        if guard == "no_change_std_guard":
            de_allocation[:] = 0.0
        elif guard == "none":
            de_allocation = (hit.astype(float) - null_hit.astype(float)) / (K * (1.0 - c))
    # Undefined branches carry null allocations. The guarded raw score is zero,
    # even though the independently shown rank-slot selection can have hits.
    raw_defined = np.isfinite(de["raw"])
    if not raw_defined:
        raw_allocation[:] = np.nan

    norm_t = float(np.linalg.norm(mu_t))
    pb = core.pb_rel_err(P, T)
    squared = np.asarray(error, dtype=float) ** 2
    pb_alloc = squared / (norm_t * norm_t) if norm_t > 0 else np.full(len(error), np.nan)
    if norm_t <= 0:
        warnings.append("Pseudobulk relative error is undefined because the truth mean has zero norm.")
    rows = []
    for g, gene in enumerate(bundle.genes):
        rows.append({
            "gene": gene,
            "mean_pred": mu_p[g], "mean_truth": mu_t[g], "mean_reference": mu_r[g],
            "change_pred": dp[g], "change_truth": dt[g], "change_error": error[g],
            "absolute_error": abs(float(error[g])), "squared_error": squared[g],
            "pb_squared_allocation": pb_alloc[g],
            "variance_pred": var_p[g], "variance_truth": var_t[g],
            "variance_ratio": float(var_p[g] / var_t[g]) if var_t[g] > 0 else None,
            "truth_de_direction": true_de[g], "pred_local_de_direction": pred_de[g],
            "local_de_hit": local_hit[g], "local_de_miss": local_miss[g],
            "local_de_false_positive": local_fp[g],
            "official_rank": ranks[g], "official_slot": slots[g],
            "official_slot_hit": hit[g], "literal_pred_direction": literal[g],
            "literal_sign_agreement": bool(literal[g] == true_literal[g]),
            "slot_hit_wrong_literal_sign": bool(hit[g] and literal[g] != true_de[g]),
            "de_raw_allocation": raw_allocation[g], "de_score_allocation": de_allocation[g],
        })
    pb_sum = float(pb_alloc.sum())
    de_sum = float(de_allocation.sum())
    summary = {
        "definition": "Mean log-expression shift and error; local significant DE calls are separate from upstream target-sized rank slots.",
        "evidence_type": "algebraic_allocation_and_diagnostic_mismatch",
        "limitations": [
            "Mean log-expression shifts are not log fold-changes of aggregated raw counts.",
            "DE rank slots do not require literal positive/negative shifts; literal signs are shown separately.",
            "DE allocations reconstruct only the guarded upstream score, not causal gene importance.",
            "Prediction-vs-reference significant DE sets are diagnostic and do not define upstream de_score.",
            "Local false-positive means a prediction-only significant signed DE call; truth nonsignificance is not no biological change or an estimate of false discovery rate.",
            "Cell-level significance is conditional on supplied cells, not independent embryo replication.",
            "Per-gene variance ratios are undefined when truth variance is zero.",
        ],
        "pb_rel_err": pb, "pb_squared_allocation_sum": pb_sum,
        "pb_reconciliation_abs_error": abs(pb_sum - pb * pb),
        "official_de_score": de["score"], "official_de_raw": de["raw"],
        "official_de_direction": direction, "de_guard": guard,
        "de_direction_guard": direction_guard,
        "n_truth_up": len(up_t), "n_truth_down": len(dn_t),
        "n_pred_local_up": len(up_p), "n_pred_local_down": len(dn_p),
        "n_slot_hits": int(hit.sum()), "n_local_hits": int(local_hit.sum()),
        "n_local_misses": int(local_miss.sum()), "n_local_false_positives": int(local_fp.sum()),
        "local_de_status_definition": "Hit: matching significant signed calls. Miss: truth signed call without matching prediction call. False-positive: prediction signed call without matching truth call; a direction reversal can be both miss and prediction-only signed call.",
        "n_slot_hits_wrong_literal_sign": sum(bool(row["slot_hit_wrong_literal_sign"]) for row in rows),
        "de_allocation_sum": de_sum,
        "de_reconciliation_abs_error": abs(de_sum - de["score"]),
        "null_orientation": null_orientation, "expression_null": c,
        "legacy_uniform_null": de["chance_uniform"],
        "signed_uniform_expectation": (len(up_t) ** 2 + len(dn_t) ** 2) / (len(bundle.genes) * K) if K else None,
        "partial_direction_definition": "Official rank-partial correlation controlling reference mean rank; literal sign counts do not reconstruct it.",
    }
    return rows, summary, warnings


def _composition(bundle: Any, scored: dict) -> tuple[list[dict], dict]:
    probe = scored["probe"]
    classes = probe.classes_
    soft = probe.predict_proba(bundle.prediction).mean(0)
    labels, counts = np.unique(bundle.target_labels, return_counts=True)
    count_by_label = dict(zip(labels, counts))
    observed_counts = np.array([count_by_label.get(c, 0) for c in classes], float)
    # The source smooths TARGET COUNTS and PREDICTED FRACTIONS before separate
    # normalization. Smoothing two fraction vectors would not match its JSD.
    p = observed_counts + 1e-6
    q = soft + 1e-6
    p, q = p / p.sum(), q / q.sum()
    mixture = 0.5 * (p + q)
    terms = 0.5 * (p * np.log2(p / mixture) + q * np.log2(q / mixture))
    official = scored["legacy"].composition_jsd(bundle.prediction, probe, bundle.target_labels)
    rows = [{
        "celltype": str(c), "target_count": int(observed_counts[i]),
        "target_fraction": observed_counts[i] / len(bundle.target_labels),
        "pred_soft_fraction": soft[i],
        "proportion_error": soft[i] - observed_counts[i] / len(bundle.target_labels),
        "smoothed_target": p[i], "smoothed_pred": q[i], "jsd_term": terms[i],
    } for i, c in enumerate(classes)]
    return rows, {
        "definition": "Frozen official target probe soft proportions; class terms sum to upstream base-2 JSD using count-vs-fraction smoothing.",
        "evidence_type": "algebraic_allocation",
        "limitations": [
            "Classifier-defined proportions depend on target labels and probe domain/calibration.",
            "Prediction labels are ignored; low soft proportion is not proof of biological absence.",
            "The probe is fitted on supplied target cells, not an independent validation population.",
        ],
        "jsd": official, "sum_terms": float(terms.sum()),
        "reconciliation_abs_error": abs(float(terms.sum()) - official),
        "smoothing_epsilon": 1e-6, "target_smoothing_input": "counts",
        "prediction_smoothing_input": "soft_fractions",
        "target_class_coverage": float(observed_counts.sum() / len(bundle.target_labels)),
    }


def _joint_ledger(bundle: Any, core: Any, seed: int) -> tuple[list[dict], dict]:
    rng = np.random.default_rng(seed)
    P, T = bundle.prediction, bundle.target
    pi = rng.choice(len(P), min(VARIOGRAM_CELL_CAP, len(P)), replace=False)
    ti = rng.choice(len(T), min(VARIOGRAM_CELL_CAP, len(T)), replace=False)
    G = P.shape[1]
    i = rng.integers(0, G, VARIOGRAM_PAIRS)
    j = rng.integers(0, G, VARIOGRAM_PAIRS)
    keep = i != j
    i, j = i[keep], j[keep]
    pair_ids = i.astype(np.int64) * G + j
    unique, inverse, multiplicity = np.unique(pair_ids, return_inverse=True, return_counts=True)
    ui, uj = unique // G, unique % G
    n_pairs = len(i)
    squared_unique = np.empty(len(unique), dtype=np.float32)
    rows: list[dict] = []
    for start in range(0, len(unique), PAIR_BLOCK):
        stop = min(start + PAIR_BLOCK, len(unique))
        left, right = ui[start:stop], uj[start:stop]
        a_left, a_right = P[np.ix_(pi, left)], P[np.ix_(pi, right)]
        b_left, b_right = T[np.ix_(ti, left)], T[np.ix_(ti, right)]
        va = (np.abs(a_left - a_right) ** 0.5).mean(0)
        vb = (np.abs(b_left - b_right) ** 0.5).mean(0)
        error = va - vb
        squared_unique[start:stop] = error ** 2
        # Descriptive population covariance on the SAME sampled cells, in
        # float64 for cancellation stability. It is not a variogram term.
        al, ar = a_left.astype(float), a_right.astype(float)
        bl, br = b_left.astype(float), b_right.astype(float)
        ca = ((al - al.mean(0)) * (ar - ar.mean(0))).mean(0)
        cb = ((bl - bl.mean(0)) * (br - br.mean(0))).mean(0)
        for k in range(stop - start):
            idx = start + k
            sq = float(squared_unique[idx])
            rows.append({
                "gene_i": bundle.genes[int(left[k])], "gene_j": bundle.genes[int(right[k])],
                "sampled_occurrences": int(multiplicity[idx]),
                "variogram_pred": va[k], "variogram_truth": vb[k],
                "variogram_error": error[k], "variogram_squared_error": sq,
                "variogram_contribution": sq * int(multiplicity[idx]) / n_pairs,
                "covariance_pred": ca[k], "covariance_truth": cb[k],
                "covariance_error": ca[k] - cb[k],
            })
    # Re-expand only the small pair-value vector in upstream occurrence order.
    # This retains the exact float32 averaging semantics without an unbounded
    # cells x 20,000 advanced-indexing array.
    replay = float(squared_unique[inverse].mean()) if n_pairs else float("nan")
    weighted = sum(row["variogram_contribution"] for row in rows)
    return rows, {
        "definition": "Sampled E|X_i-X_j|^0.5 mismatch; repeated ordered pairs retain occurrence weights. Covariance columns are separate descriptive population covariances.",
        "evidence_type": "sampled_pair_allocation_and_diagnostic_mismatch",
        "limitations": [
            "Fractional variogram moments depend on means, marginals and dependence; they are not pure covariance.",
            "Pairs are sampled with replacement; unsampled pairs have no evidence here.",
            "Covariance can reflect composition changes and does not establish gene regulation.",
            "The companion uses a bounded supplied-cell subset, not the full original files.",
        ],
        "variogram": replay, "weighted_pair_mean": weighted,
        "reconciliation_abs_error": abs(weighted - replay),
        "n_sampled_pairs": n_pairs, "n_unique_ordered_pairs": len(unique),
        "total_possible_ordered_pairs": G * (G - 1),
        "unique_ordered_pair_coverage": len(unique) / (G * (G - 1)),
        "n_genes_observed_in_pairs": int(len(np.unique(np.concatenate([ui, uj])))),
        "observed_gene_coverage": len(np.unique(np.concatenate([ui, uj]))) / G,
        "coverage_interpretation": "Coverage is of ordered nonself pairs, not modules. Unsampled pairs have no diagnostic evidence; a small gene module may have no sampled within-module pair even when every gene appears somewhere.",
        "n_pair_draws_before_self_filter": VARIOGRAM_PAIRS,
        "pair_exponent": 0.5, "covariance_ddof": 0,
        "sampled_prediction_rows": pi.tolist(), "sampled_target_rows": ti.tolist(),
        "seed": seed, "pair_block_size": PAIR_BLOCK,
        "rounded_official_reconciliation_abs_error": None,
    }


def _mmd_sensitivity(bundle: Any, core: Any, seed: int) -> dict:
    P, T = bundle.prediction, bundle.target
    mu_p, mu_t = core.pseudobulk(P), core.pseudobulk(T)
    error = mu_p - mu_t
    order = np.argsort(-np.abs(error), kind="stable")
    selected = [int(g) for g in order if error[g] != 0][:MMD_EDIT_GENE_CAP]
    counterfactual = P.copy()
    edit_rows = []
    for g in selected:
        requested = float(mu_t[g] - mu_p[g])
        proposed = P[:, g] + np.float32(requested)
        counterfactual[:, g] = np.maximum(proposed, np.float32(0.0))
        edit_rows.append({
            "gene": bundle.genes[g], "requested_shift": requested,
            "applied_shift_before_nonnegative_clip": requested,
            "nonnegative_clipped_cells": int((proposed < 0).sum()),
            "max_abs_cell_edit": float(np.max(np.abs(counterfactual[:, g] - P[:, g]))),
            "before_mean_error": float(error[g]),
        })
    mu_cf = core.pseudobulk(counterfactual)
    for g, row in zip(selected, edit_rows):
        row["achieved_mean_shift"] = float(mu_cf[g] - mu_p[g])
        row["after_mean_error"] = float(mu_cf[g] - mu_t[g])
    before = core.mmd_unbiased(P, T, seed=seed)
    after = core.mmd_unbiased(counterfactual, T, seed=seed) if selected else before
    changed = counterfactual[:, selected].astype(float) - P[:, selected].astype(float)
    rms = float(np.sqrt(np.mean(np.sum(changed ** 2, axis=1))))
    return {
        "definition": "One target-informed, nonnegative-clipped target-mean shift of at most five largest pseudobulk-error genes, followed by unchanged official MMD replay; the gene count is bounded.",
        "evidence_type": "noncausal_sensitivity",
        "limitations": [
            "This uses the target and is an audit experiment, not expected training improvement or held-out generalization.",
            "MMD change is non-additive and is not a per-gene contribution or causal effect.",
            "The edit can worsen MMD and does not restore covariance or biological states.",
            "The log-space edit can change implied library size; no independent amplitude bound is imposed.",
            "Finite-sample unbiased MMD can be negative; values are not clipped.",
        ],
        "before": before, "after": after, "change": after - before,
        "lower_is_better": True, "affected_genes": [bundle.genes[g] for g in selected],
        "bound": "at most five genes; target-mean displacement; nonnegative floor", "gene_cap": MMD_EDIT_GENE_CAP,
        "selection": "largest absolute pseudobulk mean error; stable input-panel tie order",
        "edit_rows": edit_rows, "rms_row_displacement": rms,
        "library_size_ratio_before": core.library_size_ratio(P, T),
        "library_size_ratio_after": core.library_size_ratio(counterfactual, T),
        "same_target_and_seed": True, "same_prediction_cell_count": True,
        "seed": seed, "input_mutated": False,
    }


def diagnose(bundle: Any, scored: dict, seed: int = 0) -> dict:
    """Explain the pinned task on supplied bounded arrays without mutation.

    ``scored`` supplies the authoritative core, legacy probe helpers, and fitted
    official probe. Output is strict JSON-safe and explicitly labels undefined
    values and the three kinds of evidence.
    """
    core = scored["core"]
    genes, gene_summary, warnings = _change_ledger(bundle, core)
    composition, composition_summary = _composition(bundle, scored)
    pairs, pair_summary = _joint_ledger(bundle, core, seed)
    rounded = scored["metrics"].get("variogram")
    if rounded is not None:
        pair_summary["rounded_official_reconciliation_abs_error"] = abs(pair_summary["variogram"] - rounded)
    sensitivity = _mmd_sensitivity(bundle, core, seed)
    warnings.extend([
        "All scores and diagnostics apply to the selected supplied-cell subset; they do not certify unscanned original rows.",
        "Biological states are assessed only through global gene errors and frozen-probe proportions; no state-conditional expression attribution is computed.",
        "The target-informed matrix edit is noncausal sensitivity and may alter normalization.",
    ])
    task = getattr(bundle, "task", "T1")
    if task not in ("T1", "T2", "T3"):
        raise ValueError(f"Unsupported diagnostic task: {task}")
    inputs = bundle.provenance.get("inputs", {})
    role = "wt" if task == "T3" else "reference"
    result = {
        "schema_version": "vec-scorelens-diagnostics-v2",
        "task": task,
        "scope": {
            "assessment": "sampled_subset_only", "gene_panel": "supplied_ordered_panel",
            "contrast_role": role,
            "selected_counts": {"prediction": len(bundle.prediction), "target": len(bundle.target), role: len(bundle.reference)},
            "original_counts": {name: entry["original_shape"][0] for name, entry in inputs.items() if "original_shape" in entry},
        },
        "genes": genes, "gene_summary": gene_summary,
        "composition": composition, "composition_summary": composition_summary,
        "pairs": pairs, "pair_summary": pair_summary,
        "distribution": {
            "definition": "Official distribution guard and one explicitly bounded MMD sensitivity experiment.",
            "evidence_type": "diagnostic_mismatch_and_noncausal_sensitivity",
            "limitations": ["No energy allocation or biological causal attribution is included."],
            "variance_ratio": core.variance_ratio(bundle.prediction, bundle.target),
            "mmd_sensitivity": sensitivity,
        },
        "warnings": warnings,
        "response": {"status": "not_applicable"},
        "spatial": {"status": "not_applicable"},
        "evidence": [
            {"id": "gene_means", "class": "A — EXACT / DIRECT", "source": "genes", "scope": "selected cells", "definition": "Direct mean-expression differences; squared PB allocations reconcile numerically, not causal gene importance."},
            {"id": "composition", "class": "A — EXACT / DIRECT", "source": "composition", "scope": "selected cells and fitted probe", "definition": "Per-class JSD terms using the pinned count-versus-fraction smoothing."},
            {"id": "variogram", "class": "B — SAMPLED", "source": "pairs", "scope": "recorded selected cells and sampled ordered pairs", "definition": "Occurrence-weighted fractional moment discrepancies; unsampled pairs unassessed."},
            {"id": "covariance", "class": "D — DESCRIPTIVE / ASSOCIATIONAL", "source": "pairs", "scope": "recorded pair/cell sample", "definition": "Population covariance residuals are not variogram allocations or regulation."},
            {"id": "mmd_edit", "class": "C — COUNTERFACTUAL / SENSITIVITY", "source": "distribution.mmd_sensitivity", "scope": "same selected target", "definition": "Target-informed mean-shift edit; no causal or expected training improvement claim."},
        ],
    }
    if task == "T3":
        from .response import diagnose_response
        result["response"] = diagnose_response(bundle, scored)
        result["warnings"].extend(result["response"].get("warnings", []))
        result["evidence"].append({"id": "response", "class": "A — EXACT / DIRECT", "source": "response", "scope": "selected cells and supplied WT", "definition": "WT-relative transformed-expression deltas and guarded severity regression; no paired-cell or biological-response truth."})
    if task in ("T2", "T3"):
        from .spatial import diagnose_spatial
        result["spatial"] = diagnose_spatial(bundle, scored, seed=seed)
        result["warnings"].extend(result["spatial"].get("warnings", []))
        result["evidence"].append({"id": "spatial", "class": "MIXED — see per-field classification", "source": "spatial", "scope": "selected local coordinate frames", "definition": "Direct geometry, sampled shape and descriptive spatial-expression diagnostics; never a registered developmental trajectory."})
    from .priority import prioritize
    safe_result = _finite_json(result)
    safe_result["priorities"] = prioritize(safe_result, task=task)
    return _finite_json(safe_result)
