"""Deterministic, synthetic-only planted failures; no biological data are used."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tempfile

import anndata as ad
import numpy as np
import pandas as pd


CASES = (
    "direction", "covariance", "missing", "overrepresented", "within_state",
    "scale", "magnitude", "response_magnitude", "collapse", "combined",
    "healthy", "row_permutation",
)
COUNTS = np.array([384, 288, 192, 96])
RESPONSE = np.arange(64, 96)
SELECTED = np.array([64, 65, 66, 67, 80, 81, 82, 83])
COVARIANCE = np.arange(96, 128)


def _template(seed):
    rng = np.random.default_rng(seed)
    baseline = rng.uniform(3.8, 4.2, 128)
    means = np.broadcast_to(baseline, (4, 128)).copy()
    for state in range(4):
        means[state, state * 16:(state + 1) * 16] += 2.0
    delta = np.zeros(128)
    delta[64:80] = np.linspace(0.8, 1.4, 16)
    delta[80:96] = -np.linspace(0.8, 1.4, 16)
    return means, delta


def _draw(means, counts, rng):
    """Center noise within each state, making planted means exact."""
    blocks = []
    labels = []
    for state, count in enumerate(counts):
        for attempt in range(100):
            noise = rng.normal(0.0, 0.08, (int(count), 128))
            for module in range(4):
                sl = slice(96 + 8 * module, 104 + 8 * module)
                factor = rng.normal(0.0, 0.55, (int(count), 1))
                noise[:, sl] = factor + rng.normal(0.0, 0.04, (int(count), 8))
            noise -= noise.mean(axis=0)
            block = means[state] + noise
            if np.all(block >= 0):
                break
        else:
            raise RuntimeError("Could not generate nonnegative synthetic population")
        blocks.append(block)
        labels.extend([f"state_{state}"] * int(count))
    return np.vstack(blocks), np.asarray(labels)


def _adata(matrix, *, labels=None, prefix="cell"):
    obs = pd.DataFrame(index=[f"{prefix}_{i:04d}" for i in range(matrix.shape[0])])
    if labels is not None:
        obs["celltype"] = pd.Categorical(labels)
    var = pd.DataFrame(index=[f"synthetic_gene_{i:03d}" for i in range(128)])
    result = ad.AnnData(X=np.asarray(matrix, dtype=np.float64), obs=obs, var=var)
    result.uns["synthetic_only"] = True
    return result


def _hash(matrix):
    return hashlib.sha256(np.asarray(matrix, dtype="<f8").tobytes(order="C")).hexdigest()


def create_case(case="direction", seed=0):
    """Return prediction/target/reference AnnData plus privileged planted metadata.

    Prediction and reference have no cell-type labels. Only target labels are
    provided to the evaluator; state identities in metadata are test truth only.
    """
    if case not in CASES:
        raise ValueError(f"Unknown synthetic case {case!r}; choose from {', '.join(CASES)}")
    if not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    means, delta = _template(int(seed))
    streams = np.random.SeedSequence(int(seed)).spawn(3)
    reference, labels = _draw(means, COUNTS, np.random.default_rng(streams[0]))
    target, target_labels = _draw(means + delta, COUNTS, np.random.default_rng(streams[1]))
    prediction = target.copy()
    prediction_labels = target_labels.copy()
    rng = np.random.default_rng(streams[2])
    planted_genes = []
    planted_states = []
    causes = []
    if case in ("direction", "combined"):
        prediction[:, SELECTED] -= 2.0 * delta[SELECTED]
        planted_genes.extend(SELECTED.tolist())
        causes.append("subset_sign_reversal")
    if case in ("covariance", "combined"):
        for state in range(4):
            rows = np.flatnonzero(prediction_labels == f"state_{state}")
            for gene in COVARIANCE:
                prediction[rows, gene] = prediction[rng.permutation(rows), gene]
        planted_genes.extend(COVARIANCE.tolist())
        causes.append("covariance_scrambling")
    if case in ("missing", "combined"):
        keep = prediction_labels != "state_3"
        prediction = prediction[keep]
        prediction_labels = prediction_labels[keep]
        planted_states.append("state_3")
        causes.append("missing_population")
    elif case == "overrepresented":
        rows = np.concatenate([
            rng.choice(np.flatnonzero(target_labels == f"state_{state}"), count, replace=True)
            for state, count in enumerate((288, 240, 144, 288))
        ])
        prediction = target[rows].copy()
        prediction_labels = target_labels[rows]
        planted_states.append("state_3")
        causes.append("overrepresented_population")
    elif case == "within_state":
        rows = np.flatnonzero(prediction_labels == "state_2")
        prediction[np.ix_(rows, np.arange(96, 104))] += 0.75
        planted_genes.extend(range(96, 104))
        planted_states.append("state_2")
        causes.append("within_state_expression_shift")
    elif case == "scale":
        prediction = 1.8 * prediction + 0.7
        planted_genes.extend(range(128))
        causes.append("global_affine_scale")
    elif case == "magnitude":
        prediction[:, SELECTED] += 0.8 * delta[SELECTED]
        planted_genes.extend(SELECTED.tolist())
        causes.append("subset_response_magnitude")
    elif case == "response_magnitude":
        prediction[:, RESPONSE] += delta[RESPONSE]
        planted_genes.extend(RESPONSE.tolist())
        causes.append("global_response_magnitude")
    elif case == "collapse":
        prediction = np.broadcast_to(target.mean(axis=0), target.shape).copy()
        planted_genes.extend(range(128))
        causes.append("population_collapse")
    elif case == "row_permutation":
        rows = rng.permutation(target.shape[0])
        prediction = target[rows].copy()
        prediction_labels = target_labels[rows]

    if not np.all(np.isfinite(prediction)) or np.any(prediction < 0):
        raise RuntimeError("Synthetic corruption violated finite/nonnegative prerequisite")
    genes = [f"synthetic_gene_{i:03d}" for i in range(128)]
    planted_genes = sorted(set(planted_genes))
    metadata = {
        "synthetic_only": True,
        "protocol": "PLANTED_FAILURE_PROTOCOL.md",
        "case": case,
        "seed": int(seed),
        "causes": causes,
        "planted_genes": [genes[i] for i in planted_genes],
        "planted_gene_indices": planted_genes,
        "sign_gene_indices": SELECTED.tolist() if case in ("direction", "combined") else [],
        "response_gene_indices": RESPONSE.tolist(),
        "covariance_modules": [list(range(96 + 8 * m, 104 + 8 * m)) for m in range(4)],
        "planted_states": planted_states,
        "target_counts": COUNTS.tolist(),
        "target_response": delta.tolist(),
        "precision": "float64",
        "generation": {"background": [3.8, 4.2], "marker_elevation": 2.0,
                       "independent_noise_sd": 0.08, "module_factor_sd": 0.55,
                       "module_residual_sd": 0.04, "within_state_centered": True},
        "shapes": {"prediction": list(prediction.shape), "target": list(target.shape),
                   "reference": list(reference.shape)},
        "matrix_sha256": {"prediction": _hash(prediction), "target": _hash(target),
                          "reference": _hash(reference)},
        "privileged_test_only": {"prediction_state_ids": prediction_labels.tolist()},
        "limitations": ["Expression-like synthetic floats, not a biological assay simulator",
                        "No enforced 10,000 implied library-size normalization"],
    }
    return {
        "prediction": _adata(prediction, prefix="prediction"),
        "target": _adata(target, labels=target_labels, prefix="target"),
        "reference": _adata(reference, prefix="reference"),
        "planted": metadata,
    }


def create_healthy_sample(seed=0, replicate=0):
    """Independent unchanged target sample for false-alarm calibration."""
    if replicate < 0:
        raise ValueError("replicate must be nonnegative")
    case = create_case("healthy", seed)
    means, delta = _template(seed)
    matrix, labels = _draw(means + delta, COUNTS,
                           np.random.default_rng(10_000 + 100 * seed + replicate))
    case["prediction"] = _adata(matrix, prefix="healthy_sample")
    case["planted"]["case"] = "independent_healthy_sample"
    case["planted"]["replicate"] = int(replicate)
    case["planted"]["matrix_sha256"]["prediction"] = _hash(matrix)
    case["planted"]["privileged_test_only"]["prediction_state_ids"] = labels.tolist()
    return case


def write_case(out, case="direction", seed=0):
    """Atomically create three H5ADs and truth JSON in a NEW directory."""
    out = Path(out).expanduser().absolute()
    if out.exists() or out.is_symlink():
        raise FileExistsError(f"Synthetic output already exists: {out}")
    out.parent.mkdir(parents=True, exist_ok=True)
    sample = create_case(case, seed)
    staging = Path(tempfile.mkdtemp(prefix=f".{out.name}.synthetic-", dir=out.parent))
    try:
        for name in ("prediction", "target", "reference"):
            sample[name].write_h5ad(staging / f"{name}.h5ad")
        (staging / "truth.json").write_text(json.dumps(sample["planted"], indent=2) + "\n")
        if out.exists() or out.is_symlink():
            raise FileExistsError(f"Synthetic output appeared during generation: {out}")
        staging.rename(out)
    except BaseException:
        shutil.rmtree(staging)
        raise
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--case", choices=CASES, default="direction")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)
    try:
        out = write_case(args.out, args.case, args.seed)
    except (ValueError, FileExistsError, RuntimeError) as exc:
        parser.exit(2, f"error: {exc}\n")
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
