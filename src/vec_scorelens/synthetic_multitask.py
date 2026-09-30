"""Synthetic-only T1/T2/T3 fixtures with privileged, independently usable truth.

These expression-like floats and spatial gradients describe no biological assay.
Planted fixtures use exact population means to isolate each intervention. The
independent healthy controls use ordinary uncentered noise and random counts.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import tempfile

import anndata as ad
import numpy as np
import pandas as pd

from .synthetic import CASES as T1_CASES, create_case


TASK_CASES = {
    "T1": (*T1_CASES, "independent_healthy"),
    "T2": ("translation", "rotation", "expansion", "shrinkage", "anisotropic",
           "local_scramble", "random_cube", "wrong_expression", "wrong_geometry",
           "wrong_spatial_assignment", "neighborhood_corruption", "mixed", "healthy",
           "independent_healthy"),
    "T3": ("wrong_sign", "half_magnitude", "double_magnitude", "missed_response",
           "false_positive", "wrong_distribution", "no_response", "mixed", "healthy",
           "independent_healthy"),
}
DEFAULT_CASES = {"T1": "direction", "T2": "mixed", "T3": "mixed"}
COUNTS = np.array([96, 72, 56, 32])
RESPONSE = np.arange(80, 144)
SELECTED = np.r_[80:88, 112:120]
FALSE_POSITIVE = np.arange(144, 160)
MODULES = [np.arange(160 + 16 * i, 176 + 16 * i) for i in range(8)]


def _hash(array):
    return hashlib.sha256(np.asarray(array, dtype="<f8").tobytes(order="C")).hexdigest()


def _adata(matrix, labels, coords, prefix, *, expose_labels=False):
    obs = pd.DataFrame(index=[f"{prefix}_{i:04d}" for i in range(len(matrix))])
    if expose_labels:
        obs["celltype"] = pd.Categorical(labels)
    var = pd.DataFrame(index=[f"synthetic_gene_{i:03d}" for i in range(matrix.shape[1])])
    result = ad.AnnData(np.asarray(matrix, dtype=np.float64), obs=obs, var=var)
    if coords is not None:
        result.obsm["spatial_3D"] = np.asarray(coords, dtype=np.float64)
        result.obsm["X_spatial"] = result.obsm["spatial_3D"].copy()
    result.uns["synthetic_only"] = True
    return result


def _coordinates(counts, rng):
    """Asymmetric curved clouds; each state has a spatial and expression niche."""
    centers = np.array([[-2.4, -0.8, 0.2], [-0.7, 1.2, 0.7],
                        [1.1, -0.4, 1.5], [2.8, 0.9, -0.8]])
    blocks = []
    for state, count in enumerate(counts):
        latent = rng.normal(size=(int(count), 3))
        latent *= np.array([0.45 + state * 0.05, 0.24, 0.16])
        latent[:, 2] += 0.12 * latent[:, 0] ** 2
        blocks.append(centers[state] + latent)
    return np.vstack(blocks)


def _draw(counts, coords, baseline, delta, rng, *, centered):
    labels = np.repeat([f"state_{s}" for s in range(4)], counts)
    X = np.broadcast_to(baseline + delta, (len(labels), 500)).copy()
    for state in range(4):
        rows = np.flatnonzero(labels == f"state_{state}")
        X[np.ix_(rows, np.arange(20 * state, 20 * state + 20))] += 2.5
        noise = rng.normal(0.0, 0.08, (len(rows), 500))
        for module in MODULES:
            noise[:, module] = rng.normal(0.0, 0.4, (len(rows), 1)) + rng.normal(
                0.0, 0.04, (len(rows), len(module)))
        if centered:
            noise -= noise.mean(axis=0)
        X[rows] += noise
    # Prespecified gradients retained within states as well as between states.
    for axis, genes in enumerate((np.arange(288, 320), np.arange(320, 352), np.arange(352, 384))):
        X[:, genes] += 0.35 * coords[:, axis, None]
    if np.min(X) < 0:
        raise RuntimeError("synthetic draw unexpectedly became negative")
    return X, labels


def _independent_t1(seed):
    from .synthetic import COUNTS as old_counts, _template
    sample = create_case("healthy", seed)
    rng = np.random.default_rng(np.random.SeedSequence([seed, 171]))
    means, delta = _template(seed)
    for name, response in (("reference", 0), ("target", 1), ("prediction", 1)):
        counts = rng.multinomial(int(old_counts.sum()), old_counts / old_counts.sum())
        labels = np.repeat([f"state_{s}" for s in range(4)], counts)
        noise = rng.normal(0, 0.08, (len(labels), 128))
        for m in range(4):
            genes = np.arange(96 + 8 * m, 104 + 8 * m)
            noise[:, genes] = rng.normal(0, 0.55, (len(labels), 1)) + rng.normal(0, 0.04, (len(labels), 8))
        X = means[np.repeat(np.arange(4), counts)] + response * delta + noise
        sample[name] = _adata(X, labels, None, name, expose_labels=name == "target")
        sample["planted"]["matrix_sha256"][name] = _hash(X)
        if name == "target":
            sample["planted"]["target_counts"] = counts.tolist()
        if name == "prediction":
            sample["planted"]["privileged_test_only"]["prediction_state_ids"] = labels.tolist()
    sample["planted"].update(task="T1", case="independent_healthy", causes=[])
    sample["planted"]["generation"].update(within_state_centered=False, random_mixture_counts=True)
    return sample


def create_task_case(task, case=None, seed=0):
    """Return AnnData inputs and privileged truth; only target exposes labels."""
    task = str(task).upper()
    if task not in TASK_CASES:
        raise ValueError("task must be T1, T2, or T3")
    case = DEFAULT_CASES[task] if case is None else case
    if case not in TASK_CASES[task]:
        raise ValueError(f"Unknown synthetic {task} case {case!r}; choose from {', '.join(TASK_CASES[task])}")
    if not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    seed = int(seed)
    if task == "T1":
        if case == "independent_healthy":
            return _independent_t1(seed)
        sample = create_case(case, seed)
        sample["planted"]["task"] = "T1"
        return sample
    rng = np.random.default_rng(np.random.SeedSequence([seed, 293]))
    baseline = rng.uniform(5.8, 6.2, 500)
    delta = np.zeros(500)
    if task == "T3":
        delta[80:112] = np.linspace(0.9, 1.5, 32)
        delta[112:144] = -np.linspace(0.9, 1.5, 32)
    independent = case == "independent_healthy"
    target_counts = rng.multinomial(256, COUNTS / 256) if independent else COUNTS.copy()
    target_coords = _coordinates(target_counts, rng)
    reference_counts = rng.multinomial(256, COUNTS / 256) if independent else COUNTS.copy()
    reference_coords = _coordinates(reference_counts, rng) if independent else target_coords.copy()
    reference, reference_labels = _draw(reference_counts, reference_coords, baseline, np.zeros(500), rng, centered=not independent)
    target, target_labels = _draw(target_counts, target_coords, baseline, delta, rng, centered=not independent)
    prediction = target.copy()
    prediction_labels = target_labels.copy()
    prediction_coords = target_coords.copy()
    if independent:
        counts = rng.multinomial(256, COUNTS / 256)
        prediction_coords = _coordinates(counts, rng)
        prediction, prediction_labels = _draw(counts, prediction_coords, baseline, delta, rng, centered=False)
    planted = []
    causes = []
    changed_coordinate_rows = []
    if task == "T3":
        factors = {"wrong_sign": -1.0, "half_magnitude": 0.5, "double_magnitude": 2.0,
                   "missed_response": 0.0, "mixed": -1.0}
        if case in factors:
            affected = RESPONSE if case in ("half_magnitude", "double_magnitude") else SELECTED
            prediction[:, affected] += (factors[case] - 1.0) * delta[affected]
            planted.extend(affected.tolist())
            causes.append(case if case != "mixed" else "wrong_sign")
        if case in ("false_positive", "mixed"):
            prediction[:, FALSE_POSITIVE] += 1.1
            planted.extend(FALSE_POSITIVE.tolist())
            causes.append("false_positive")
        if case == "no_response":
            prediction = reference.copy()
            prediction_labels = reference_labels.copy()
            prediction_coords = reference_coords.copy()
            planted.extend(RESPONSE.tolist())
            causes.append("no_response")
        if case in ("wrong_distribution", "mixed"):
            for state in range(4):
                rows = np.flatnonzero(prediction_labels == f"state_{state}")
                # Permute each gene separately: means and marginal distributions
                # are preserved while planted module covariances are destroyed.
                for gene in np.concatenate(MODULES[:4]):
                    prediction[rows, gene] = prediction[rng.permutation(rows), gene]
            planted.extend(np.concatenate(MODULES[:4]).tolist())
            causes.append("covariance_scrambling")
    else:
        if case == "translation":
            prediction_coords += np.array([11.0, -8.0, 3.0])
        elif case == "rotation":
            theta = 0.73
            rotation = np.array([[np.cos(theta), -np.sin(theta), 0],
                                 [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
            prediction_coords = prediction_coords @ rotation.T
        elif case in ("expansion", "shrinkage"):
            factor = 1.7 if case == "expansion" else 0.55
            prediction_coords *= factor
            causes.append(case)
        elif case in ("anisotropic", "wrong_geometry", "mixed"):
            prediction_coords *= np.array([1.8, 0.55, 1.15])
            causes.append("anisotropic_geometry")
        elif case == "random_cube":
            lo, hi = target_coords.min(axis=0), target_coords.max(axis=0)
            prediction_coords = rng.uniform(lo, hi, target_coords.shape)
            causes.append("random_cube_geometry")
        if case in ("wrong_spatial_assignment", "mixed"):
            # Coordinates and expression are individually unchanged for pure
            # assignment corruption; their correspondence is the intervention.
            order = rng.permutation(len(prediction_coords))
            prediction_coords = prediction_coords[order]
            changed_coordinate_rows = np.flatnonzero(order != np.arange(len(order))).tolist()
            causes.append("coordinate_expression_assignment_scrambling")
        elif case == "local_scramble":
            # A prespecified left half region loses its local assignment;
            # coordinates outside this region remain exactly in place.
            rows = np.flatnonzero(target_coords[:, 0] < np.median(target_coords[:, 0]))
            prediction_coords[rows] = prediction_coords[rng.permutation(rows)]
            changed_coordinate_rows = rows.tolist()
            causes.append("left_half_coordinate_assignment_scrambling")
        elif case == "neighborhood_corruption":
            # Geometric corruption in one state niche is distinct from an
            # assignment permutation: coordinate values themselves change.
            rows = np.flatnonzero(target_labels == "state_1")
            prediction_coords[rows] += rng.normal(0, [1.0, 0.7, 0.6], (len(rows), 3))
            changed_coordinate_rows = rows.tolist()
            causes.append("state_1_local_coordinate_jitter")
        if case in ("wrong_expression", "mixed"):
            prediction[:, SELECTED] += 1.3
            planted.extend(SELECTED.tolist())
            causes.append("expression_shift")
    gene_names = [f"synthetic_gene_{i:03d}" for i in range(500)]
    planted = sorted(set(planted))
    metadata = {
        "synthetic_only": True, "task": task, "case": case, "seed": seed, "causes": causes,
        "planted_gene_indices": planted, "planted_genes": [gene_names[i] for i in planted],
        "response_gene_indices": RESPONSE.tolist() if task == "T3" else [],
        "sign_gene_indices": SELECTED.tolist() if task == "T3" and case in ("wrong_sign", "mixed") else [],
        "magnitude_gene_indices": RESPONSE.tolist() if task == "T3" and case in ("half_magnitude", "double_magnitude") else [],
        "missed_gene_indices": (RESPONSE if case == "no_response" else SELECTED).tolist() if case in ("missed_response", "no_response") else [],
        "false_positive_gene_indices": FALSE_POSITIVE.tolist() if task == "T3" and case in ("false_positive", "mixed") else [],
        "covariance_modules": [m.tolist() for m in MODULES], "target_response": delta.tolist(),
        "target_counts": target_counts.tolist(), "planted_states": [], "precision": "float64",
        "generation": {"within_state_centered": not independent, "random_mixture_counts": independent,
                       "common_random_numbers_across_T2_T3": True,
                       "independent_noise_sd": 0.08, "module_factor_sd": 0.4,
                       "spatial_gradient_gene_indices": list(range(288, 384)),
                       "spatial_dimensions": 3},
        "shapes": {"prediction": list(prediction.shape), "target": list(target.shape), "reference": list(reference.shape)},
        "matrix_sha256": {"prediction": _hash(prediction), "target": _hash(target), "reference": _hash(reference)},
        "coordinate_sha256": {"prediction": _hash(prediction_coords), "target": _hash(target_coords), "reference": _hash(reference_coords)},
        "privileged_test_only": {"prediction_state_ids": prediction_labels.tolist(),
                                 "changed_coordinate_rows": changed_coordinate_rows},
        "case_alias_of": "anisotropic" if task == "T2" and case == "wrong_geometry" else None,
        "limitations": ["Expression-like synthetic floats, not a biological assay simulator",
                        "No enforced 10,000 implied library-size normalization",
                        "T2 rigid controls do not assume official occupancy/frame alignment invariance"],
    }
    return {"prediction": _adata(prediction, prediction_labels, prediction_coords, "prediction"),
            "target": _adata(target, target_labels, target_coords, "target", expose_labels=True),
            "reference": _adata(reference, reference_labels, reference_coords, "reference"), "planted": metadata}


def write_task_case(out, *, task, case=None, seed=0):
    """Atomically create prediction/target/reference H5AD and truth in a NEW directory."""
    out = Path(out).expanduser().absolute()
    if out.exists() or out.is_symlink():
        raise FileExistsError(f"Synthetic output already exists: {out}")
    sample = create_task_case(task, case, seed)
    out.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{out.name}.synthetic-", dir=out.parent))
    try:
        for name in ("prediction", "target", "reference"):
            sample[name].write_h5ad(staging / f"{name}.h5ad")
        (staging / "truth.json").write_text(json.dumps(sample["planted"], indent=2, allow_nan=False) + "\n")
        if out.exists() or out.is_symlink():
            raise FileExistsError(f"Synthetic output appeared during generation: {out}")
        staging.rename(out)
    except BaseException:
        shutil.rmtree(staging)
        raise
    return out
