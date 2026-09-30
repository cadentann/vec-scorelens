"""Pinned-source spatial explanations; allocations are not biological causes.

These diagnostics apply to the same selected cells as the official panel.
Independent spatial samples have no assumed cellular correspondence.
"""
from __future__ import annotations

from typing import Any
import numpy as np


MAX_SPATIAL_GENES = 2_000
MAX_SPATIAL_CELLS = 1_024
D2_PAIRS = 200_000
PROJECTIONS = 200
SPATIAL_K = 15
FRAME_WARNING = (
    "Pinned sliced Wasserstein and occupancy Dice have a PCA-frame handedness defect: "
    "proper rotation or independent row permutation can change their scores, even for identical clouds. "
    "Treat these as provisional frame-dependent diagnostics, not standalone failure priorities or laterality tests."
)


def _json(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_json(v) for v in value]
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    return value


def _section(key, value, scored, definition, evidence, rows, **extra):
    official = scored.get("metrics", {}).get(key)
    finite = bool(np.isfinite(value))
    return {
        "official_metric_key": key, "official_value": official,
        "replayed_value": value, "definition": definition,
        "evidence_type": evidence, "scope": "selected_supplied_cells",
        "status": "finite" if finite else "undefined",
        "rounded_official_abs_error": abs(value - official) if official is not None and finite else None,
        "rows": rows, **extra,
    }


def _d2(P, T, shape, scored, seed):
    rng = np.random.default_rng(seed)

    def distances(C):
        X = np.asarray(C, float) - np.asarray(C, float).mean(0)
        rms = float(np.sqrt((X ** 2).sum(1).mean()))
        if rms < 1e-12:
            return np.zeros(1)
        i = rng.integers(0, len(X), D2_PAIRS)
        j = rng.integers(0, len(X), D2_PAIRS)
        keep = i != j
        return np.linalg.norm(X[i[keep]] - X[j[keep]], axis=1) / rms

    a, b = distances(P), distances(T)
    sa, sb = np.sort(a), np.sort(b)
    edges = np.unique(np.concatenate([sa, sb]))
    left, right = edges[:-1], edges[1:]
    delta = np.searchsorted(sa, left, side="right") / len(sa) - np.searchsorted(sb, left, side="right") / len(sb)
    terms = np.abs(delta) * (right - left)
    # Aggregate exact CDF intervals, never approximate distance distributions
    # with histogram-bin centers. The rows still sum to the exact W1 integral.
    rows = []
    groups = np.array_split(np.arange(len(left)), min(32, len(left))) if len(left) else []
    for indices in groups:
        lo, hi = float(left[indices[0]]), float(right[indices[-1]])
        widths = right[indices] - left[indices]
        rows.append({
            "distance_lo": lo, "distance_hi": hi,
            "contribution": float(terms[indices].sum()),
            "cdf_delta_mean": float(np.average(delta[indices], weights=widths)),
            "pred_pairs": int(((a >= lo) & (a < hi)).sum()),
            "truth_pairs": int(((b >= lo) & (b < hi)).sum()),
        })
    if rows:
        # Last interval is closed on its right endpoint for descriptive counts.
        rows[-1]["pred_pairs"] += int((a == edges[-1]).sum())
        rows[-1]["truth_pairs"] += int((b == edges[-1]).sum())
    value = shape.d2_distance(P, T, seed=seed)
    total = float(terms.sum())
    return _section("d2_shape", value, scored,
        "Wasserstein distance between independently sampled RMS-normalized within-cloud pair distances.",
        "sampled_distribution_allocation", rows,
        effective_seed=seed, prediction_pairs=len(a), truth_pairs=len(b),
        draws_per_side=D2_PAIRS, allocation_sum=total,
        reconciliation_abs_error=abs(total-value),
        limitations=[
            "Independent sampled pairs give a nonzero identical-cloud Monte Carlo floor; row-order changes can change the sampled realization.",
            "Distance-band contributions are not anatomical regions, homologous cell pairs, or causal attributions.",
            "D2 is reflection-blind; finite distance distributions do not uniquely identify a shape.",
        ])


def _frame(C):
    X = np.asarray(C, float) - np.asarray(C, float).mean(0)
    _, s, vt = np.linalg.svd(X-X.mean(0), full_matrices=False)
    eigen = s*s / len(X)
    return {"determinant": float(np.linalg.det(vt)), "eigenvalues": eigen,
            "relative_eigenvalue_gaps": np.diff(-eigen) / max(float(eigen[0]), 1e-12)}


def _shape_grid(P, T, shape, scored, seed):
    A, B, n = shape._match_n(np.asarray(P, float), np.asarray(T, float), seed)
    Za, _ = shape._canonicalise(A)
    Zb, _ = shape._canonicalise(B)
    rng = np.random.default_rng(seed)
    directions = rng.normal(size=(PROJECTIONS, 3))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    truth_proj = np.sort(Zb @ directions.T, axis=0)
    flip_values, projection_values = [], []
    for f in shape._PROPER_FLIPS:
        residual = np.abs(np.sort((Za*np.asarray(f)) @ directions.T, axis=0)-truth_proj)
        flip_values.append(float(residual.mean()))
        projection_values.append(residual.mean(0))
    best = int(np.argmin(flip_values))
    selected = projection_values[best]
    sw, spread = shape.sliced_wasserstein(P, T, seed=seed)
    frames = {"prediction": _frame(A), "truth": _frame(B)}
    sw_rows = [{"projection": i, "direction": directions[i],
                "distance": selected[i], "contribution": selected[i]/PROJECTIONS}
               for i in range(PROJECTIONS)]
    sw_report = _section("sliced_wasserstein", sw, scored,
        "Mean projected sorted-rank distance at the one globally selected upstream sign flip.",
        "algebraic_allocation", sw_rows, status="provisional_frame",
        effective_seed=seed, matched_cells=n, frames=frames,
        selected_flip=shape._PROPER_FLIPS[best], flip_values=flip_values, flip_spread=spread,
        allocation_sum=float(selected.mean()), reconciliation_abs_error=abs(float(selected.mean())-sw),
        limitations=[FRAME_WARNING, "Projection rank matching is not biological cell correspondence.",
                     "Flip spread is not a calibrated uncertainty interval or sufficient test for eigenvalue degeneracy."])

    def occupancy(Z):
        vox = np.floor((np.clip(Z, -3, 3-1e-9)+3)/6*16).astype(int)
        unique, counts = np.unique(vox, axis=0, return_counts=True)
        return {tuple(v): int(c) for v, c in zip(unique, counts)}

    tb = occupancy(Zb)
    alternatives, scores = [], []
    for f in shape._PROPER_FLIPS:
        pa = occupancy(Za*np.asarray(f))
        alternatives.append(pa)
        scores.append(2*len(set(pa) & set(tb))/max(len(pa)+len(tb), 1))
    best_dice = int(np.argmax(scores))
    pa = alternatives[best_dice]
    denom = len(pa)+len(tb)
    voxel_rows = []
    for voxel in sorted(set(pa) | set(tb)):
        inp, intruth = voxel in pa, voxel in tb
        voxel_rows.append({"voxel": voxel, "lower_bound": np.array(voxel)*.375-3,
            "upper_bound": (np.array(voxel)+1)*.375-3,
            "status": "shared" if inp and intruth else "prediction_only" if inp else "truth_only",
            "prediction_points": pa.get(voxel, 0), "truth_points": tb.get(voxel, 0),
            "loss_contribution": 0.0 if inp and intruth else 1.0/denom})
    dice, resolution = shape.occupancy_dice(P, T, seed=seed)
    loss_sum = sum(r["loss_contribution"] for r in voxel_rows)
    selected_Za = Za*np.asarray(shape._PROPER_FLIPS[best_dice])
    dice_report = _section("occupancy_dice", dice, scored,
        "Binary occupied-voxel Dice on the pinned canonical grid; unmatched voxels allocate 1−Dice.",
        "algebraic_allocation", voxel_rows, status="provisional_frame",
        effective_seed=seed, matched_cells=n, selected_flip=shape._PROPER_FLIPS[best_dice],
        flip_values=scores, frames=frames, voxel_edge_over_truth_nn=resolution,
        resolution_warning=bool(resolution <= 10),
        prediction_clipped_fraction=float(np.any((selected_Za < -3)|(selected_Za > 3-1e-9), axis=1).mean()),
        truth_clipped_fraction=float(np.any((Zb < -3)|(Zb > 3-1e-9), axis=1).mean()),
        allocation_target="1 - occupancy_dice", allocation_sum=loss_sum,
        reconciliation_abs_error=abs(loss_sum-(1-dice)),
        limitations=[FRAME_WARNING, "Voxels are metric-grid regions, not registered anatomical regions.",
                     "Voxel/NN ratio below roughly 10 risks measuring sampling density; clipped outliers occupy boundary voxels.",
                     "Point counts are descriptive, not Dice weights."])
    return sw_report, dice_report


def _growth(bundle, shape, scored):
    P, T = bundle.prediction_coords, bundle.target_coords
    radial = []
    radii = []
    for side, C in (("prediction", P), ("truth", T)):
        X = np.asarray(C, float)-np.asarray(C, float).mean(0)
        squared = (X*X).sum(1)
        rms = float(np.sqrt(squared.mean()))
        radii.append(rms)
        shell = np.minimum(np.floor(np.sqrt(squared)/max(rms, 1e-12)).astype(int), 3)
        for r in range(4):
            mask = shell == r
            radial.append({"side": side, "rms_shell": r, "cells": int(mask.sum()),
                           "squared_radius_contribution": float(squared[mask].sum()/len(C))})
    scale = shape.scale_log_ratio(P, T)
    scale_report = _section("scale_log_ratio", scale, scored,
        "Signed log RMS-radius ratio; shell contributions sum to each side's squared RMS, not to the log ratio.",
        "descriptive_mismatch", radial, prediction_rms=radii[0], target_rms=radii[1],
        reconciliation_abs_error=max(abs(sum(r["squared_radius_contribution"] for r in radial if r["side"] == side)-radius**2)
            for side, radius in zip(("prediction", "truth"), radii)),
        limitations=["RMS extent is not occupied volume, biological growth, or anatomical region error.",
                     "Coordinate units must match; all coordinates here are selected-scope."])
    npred, ntruth = len(bundle.prediction), len(bundle.target)
    inputs = bundle.provenance.get("inputs", {})
    def original(name):
        return inputs.get(name, {}).get("original_shape", [None])[0]
    op, ot = original("prediction"), original("target")
    count = shape.count_log_ratio(npred, ntruth)
    count_report = _section("count_log_ratio", count, scored,
        "Signed log ratio of supplied selected expression row counts; original-file counts are separate metadata diagnostics.",
        "descriptive_mismatch", [], prediction_count=npred, target_count=ntruth,
        original_prediction_count=op, original_target_count=ot,
        original_log_ratio=shape.count_log_ratio(op, ot) if op is not None and ot is not None else None,
        original_relative_count_difference=op/ot-1 if op is not None and ot else None,
        reconciliation_abs_error=0.0,
        limitations=["Equal cell caps can hide original-file count differences in the selected-scope official metric.",
                     "Observed cell counts do not by themselves establish proliferation or cell loss."])
    return {"scale": scale_report, "count": count_report}


def _moran(bundle, legacy, scored):
    P, T = bundle.prediction, bundle.target
    genes = np.argsort(-T.var(0))[:200]
    ip = legacy.morans_I_profile(P, bundle.prediction_coords, genes=genes)
    it = legacy.morans_I_profile(T, bundle.target_coords, genes=genes)
    up, ut = ip-ip.mean(), it-it.mean()
    denom = float(np.sqrt((up*up).sum()*(ut*ut).sum()))
    contributions = up*ut/denom if denom > 0 else np.full(len(genes), np.nan)
    value = legacy.morans_I_agreement(P, bundle.prediction_coords, T, bundle.target_coords)
    rows = [{"gene": bundle.genes[g], "gene_index": int(g),
        "prediction_I": ip[i], "truth_I": it[i], "difference": ip[i]-it[i],
        "agreement_contribution": contributions[i]} for i, g in enumerate(genes)]
    return _section("morans_I_agreement", value, scored,
        "Pearson agreement of Moran's I on genes selected once by truth variance; centered gene products allocate correlation.",
        "algebraic_allocation_and_descriptive_mismatch", rows,
        allocation_sum=float(contributions.sum()),
        reconciliation_abs_error=abs(float(contributions.sum())-value),
        selected_genes=len(genes), small_profile_warning=bool(len(genes)<10),
        official_panel_includes_metric="morans_I_agreement" in scored.get("metrics", {}),
        limitations=["Autocorrelation strength agreement does not establish matching expression locations or laterality.",
                     "Constant I profiles make correlation undefined; a small gene profile can give uninformative extreme correlations.",
                     "Moran drops the first of 16 queried neighbors; duplicate-coordinate ties may affect self-exclusion."])


def _mmd_anchors(hp, ht, bundle, official_value):
    """Signed two-sided kernel allocation in the exact upstream fitted space.

    Kernels are evaluated sequentially. Bandwidth distances use the upstream
    subtraction/square/sum expression in row blocks, retaining its dtype and
    within-distance reduction, without an N x N x PC allocation.
    """
    from sklearn.decomposition import PCA
    from sklearn.metrics.pairwise import rbf_kernel

    scales = (.25, .5, 1., 2., 4.)
    rng = np.random.default_rng(0)
    pca = PCA(n_components=min(30, ht.shape[1]), random_state=0).fit(ht)
    pi = rng.choice(len(hp), min(2000, len(hp)), replace=False)
    ti = rng.choice(len(ht), min(2000, len(ht)), replace=False)
    A, B = pca.transform(hp[pi]), pca.transform(ht[ti])
    na, nb = len(A), len(B)
    d2 = np.empty((nb, nb), dtype=B.dtype)
    for start in range(0, nb, 64):
        d2[start:start+64] = np.sum((B[start:start+64, None]-B[None, :])**2, -1)
    positive = d2[d2 > 0]
    gamma0 = 1.0/(np.median(positive)+1e-9) if len(positive) else float("nan")
    del d2, positive
    ap, at = np.zeros(na, dtype=float), np.zeros(nb, dtype=float)
    source_total = 0.0
    if np.isfinite(gamma0):
        for scale in scales:
            kernel = rbf_kernel(A, A, gamma0*scale)
            np.fill_diagonal(kernel, 0.)
            source_pp = kernel.sum()/(na*(na-1))
            ap += kernel.sum(axis=1, dtype=np.float64)/(na*(na-1))
            del kernel
            kernel = rbf_kernel(B, B, gamma0*scale)
            np.fill_diagonal(kernel, 0.)
            source_tt = kernel.sum()/(nb*(nb-1))
            at += kernel.sum(axis=1, dtype=np.float64)/(nb*(nb-1))
            del kernel
            kernel = rbf_kernel(A, B, gamma0*scale)
            source_cross = kernel.mean()
            ap -= kernel.sum(axis=1, dtype=np.float64)/(na*nb)
            at -= kernel.sum(axis=0, dtype=np.float64)/(na*nb)
            source_total += source_pp+source_tt-2*source_cross
            del kernel
        ap /= len(scales)
        at /= len(scales)
        source_value = float(source_total/len(scales))
    else:
        ap[:] = np.nan
        at[:] = np.nan
        source_value = float("nan")

    anchor_rows = []
    selected_rows = bundle.provenance.get("selection", {}).get("rows", {})
    for side, role, indices, allocations in (("prediction", "prediction", pi, ap), ("truth", "target", ti, at)):
        coords = getattr(bundle, role+"_coords")
        provenance_rows = selected_rows.get(role, [])
        for index, allocation in zip(indices, allocations):
            row = {"side": side, "row_position": int(index),
                   "coordinates": coords[index].tolist(), "allocation": float(allocation)}
            if len(provenance_rows) > index:
                row["original_row_position"] = provenance_rows[index].get("position")
                row["obs_name"] = provenance_rows[index].get("obs_name")
            anchor_rows.append(row)
    total = float(ap.sum()+at.sum())
    residual = abs(total-official_value)
    # Native float32 global reductions in upstream and float64 row reductions
    # here are different summation orders, not different kernel estimators.
    tolerance = max(1e-12, 32*float(np.finfo(A.dtype).eps))
    return anchor_rows, {
        "definition": "Both sides' within-kernel means minus cross-kernel means, each normalized by its own anchor count; signed terms sum to MMD-u.",
        "evidence_type": "algebraic_allocation", "sum": total,
        "reconciliation_abs_error": residual, "tolerance": tolerance,
        "passed": bool(np.isfinite(residual) and residual <= tolerance),
        "status": "finite" if np.isfinite(total) else "undefined_truth_kernel",
        "source_replayed_value": source_value,
        "source_replay_abs_error": abs(source_value-official_value),
        "effective_seed": 0, "kernel_scales": scales, "gamma0": gamma0,
        "n_prediction_anchors": na, "n_truth_anchors": nb,
        "coordinate_frame": "each_clouds_own_supplied_frame",
        "limitations": [
            "Signed allocations can be negative; magnitude is a kernel term, not nonnegative error, causal importance, or a failure priority.",
            "An anchor identifies a neighborhood center in its own cloud, not a matched prediction-to-truth cell or anatomical region.",
            "Both prediction and truth terms are required for reconciliation; gene residual rows remain nonadditive.",
            "Float64 anchor reductions reconcile to native-dtype upstream reductions within the reported tolerance.",
        ],
    }


def _local(bundle, core, legacy, scored, seed):
    P, T = bundle.prediction, bundle.target
    pc, tc = bundle.prediction_coords, bundle.target_coords
    hp = legacy._knn_neighborhood_pb(P, pc)
    ht = legacy._knn_neighborhood_pb(T, tc)
    value = core.mmd_unbiased(hp, ht, seed=0)
    anchor_rows, anchor_allocation = _mmd_anchors(hp, ht, bundle, value)
    mp, mt = hp.mean(0), ht.mean(0)
    rows = [{"gene": gene, "prediction_neighborhood_mean": mp[g],
             "truth_neighborhood_mean": mt[g], "mean_residual": mp[g]-mt[g]}
            for g, gene in enumerate(bundle.genes)]
    # One bounded, deterministic reassignment preserves coordinates as a set and
    # every expression value while testing dependence on their supplied pairing.
    permutation = np.random.default_rng(seed).permutation(len(pc))
    cf_h = legacy._knn_neighborhood_pb(P, pc[permutation])
    cf = core.mmd_unbiased(cf_h, ht, seed=0)
    sensitivity = {"evidence_type": "noncausal_sensitivity", "operation": "prediction_coordinate_assignment_scramble",
        "permutation": permutation, "seed": seed, "effective_mmd_seed": 0,
        "before": value, "after": cf, "change": cf-value,
        "same_expression_values": True, "same_coordinate_set": True, "input_mutated": False,
        "limitations": ["This counterfactual breaks expression-location assignment; it is not a predicted intervention or expected model improvement.",
                        "No pointwise correspondence between prediction and truth is assumed; a scramble need not worsen every input."]}
    report = _section("neighborhood_mmd", value, scored,
        "Truth-fitted PCA multi-kernel MMD-u of spatial-neighborhood pseudobulks; gene means are descriptive residuals.",
        "algebraic_anchor_allocation_descriptive_gene_mismatch_and_noncausal_sensitivity", rows,
        effective_seed=0, neighbors=SPATIAL_K, includes_self_for_unique_coordinates=True,
        anchor_rows=anchor_rows, anchor_allocation=anchor_allocation,
        sensitivity=sensitivity, reconciliation_abs_error=anchor_allocation["reconciliation_abs_error"],
        reconciliation_status="signed_two_sided_anchors_sum_to_scalar_gene_means_are_nonadditive",
        limitations=["Gene residuals are not additive MMD contributions; no gene-wise reconciliation is claimed.",
                     "MMD-u can be negative, including for identical arrays; no clipping to zero.",
                     "Neighborhoods are constructed after cell selection and may differ from full-file neighborhoods.",
                     "Constant neighborhood truth can make the truth-fitted kernel undefined even when original expression varies."])
    return {"neighborhood_mmd": report, "moran": _moran(bundle, legacy, scored)}


def diagnose_spatial(bundle: Any, scored: dict, seed: int = 0) -> dict:
    """Explain T2/T3 spatial scores without editing input arrays or scorer code."""
    from .io import InputPolicyError

    if max(len(bundle.prediction), len(bundle.target)) > MAX_SPATIAL_CELLS or bundle.prediction.shape[1] > MAX_SPATIAL_GENES:
        raise InputPolicyError("spatial diagnostics support at most 1024 selected cells per side and 2000 genes")
    for name in ("prediction", "target"):
        C = np.asarray(getattr(bundle, name+"_coords"))
        if C.shape != (len(getattr(bundle, name)), 3) or not np.isfinite(C).all():
            raise InputPolicyError(f"{name} coordinates must be finite selected cells × 3")
        if len(C) <= SPATIAL_K:
            raise InputPolicyError("spatial diagnostics require at least 16 selected cells per side")
    shape, core, legacy = scored["shape"], scored["core"], scored["spatial_legacy"]
    P, T = bundle.prediction_coords, bundle.target_coords
    sw, dice = _shape_grid(P, T, shape, scored, seed)
    return _json({"shape": {"d2": _d2(P, T, shape, scored, seed), "sliced_wasserstein": sw, "occupancy": dice},
        "growth": _growth(bundle, shape, scored), "local": _local(bundle, core, legacy, scored, seed),
        "warnings": [FRAME_WARNING,
            "Spatial scores and neighborhoods apply only to selected supplied cells; full-file count is a separate descriptive diagnostic.",
            "D2 uses independent random pair draws and can have a nonzero identical-cloud floor.",
            "No anatomical correspondence, regional causation, or calibrated laterality is assessed."]})
