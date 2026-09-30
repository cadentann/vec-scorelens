"""Bounded, selective HDF5 loading for ScoreLens inputs."""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any

import h5py
import numpy as np
from scipy import sparse

MIN_CELLS = 32
MIN_GENES = 2
MAX_CELLS = 1024
MAX_GENES = 40_000
MAX_MULTITASK_GENES = 2_000
MAX_DENSE_MIB = 256
_MIB = 1024 * 1024
_CSC_CHUNK_NNZ = 262_144
MAX_LABEL_VALUES = 1_000_000


@dataclass(frozen=True)
class Bundle:
    """Immutable selected inputs plus JSON-safe source provenance."""
    prediction: np.ndarray
    target: np.ndarray
    reference: np.ndarray
    genes: list[str]
    target_labels: np.ndarray
    provenance: dict[str, Any]
    task: str = "T1"
    prediction_coords: np.ndarray | None = None
    target_coords: np.ndarray | None = None
    reference_coords: np.ndarray | None = None


class InputPolicyError(ValueError):
    """The supplied file is outside ScoreLens' narrow supported-input policy."""


@dataclass(frozen=True)
class _Meta:
    path: Path
    n_obs: int
    n_vars: int
    dtype: np.dtype
    genes: list[str]
    layout: str
    labels: np.ndarray | None


def require_target_variation(target: np.ndarray) -> None:
    if target.ndim != 2 or target.shape[0] < 2:
        raise InputPolicyError("target must contain at least two selected rows for the pinned T1 scorer")
    if not np.any(target != target[0]):
        raise InputPolicyError("target selected rows are all identical; pinned target-fitted MMD requires within-target variation")


def _fingerprint(path: Path) -> dict[str, Any]:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    stat = path.stat()
    return {"path": str(path.resolve()), "sha256": digest.hexdigest(), "bytes": int(stat.st_size), "mtime_ns": int(stat.st_mtime_ns)}


def _text(value: Any) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def _encoding(obj: h5py.Group | h5py.Dataset) -> str:
    value = obj.attrs.get("encoding-type")
    return _text(value) if value is not None else ""


def _local_member(group: h5py.Group | h5py.File, key: str, context: str) -> h5py.Group | h5py.Dataset | None:
    """Reject links and dataset storage that escape the fingerprinted H5AD file."""
    link = group.get(key, getlink=True)
    if isinstance(link, (h5py.ExternalLink, h5py.SoftLink)):
        raise InputPolicyError(f"{context} must use an in-file hard link; external and soft HDF5 links are unsupported")
    value = group.get(key)
    if isinstance(value, h5py.Dataset):
        if value.is_virtual or value.id.get_create_plist().get_external_count():
            raise InputPolicyError(f"{context} uses external or virtual dataset storage; only in-file storage is supported")
    return value


def _dataset(group: h5py.Group, key: str, context: str) -> h5py.Dataset:
    value = _local_member(group, key, f"{context}.{key}")
    if not isinstance(value, h5py.Dataset):
        raise InputPolicyError(f"{context} requires dataset {key!r}")
    return value


def _real_dtype(dtype: np.dtype, name: str) -> np.dtype:
    dtype = np.dtype(dtype)
    if not np.issubdtype(dtype, np.number) or np.issubdtype(dtype, np.complexfloating) or np.issubdtype(dtype, np.bool_):
        raise InputPolicyError(f"{name}.X has unsupported dtype {dtype}; expected real numeric values")
    return dtype


def _index_dataset(root: h5py.File, axis: str, name: str) -> h5py.Dataset:
    frame = _local_member(root, axis, f"{name}.{axis}")
    if not isinstance(frame, h5py.Group):
        raise InputPolicyError(f"{name}.{axis} must be a supported H5AD dataframe")
    index_name = frame.attrs.get("_index")
    if index_name is None:
        raise InputPolicyError(f"{name}.{axis} has no supported dataframe index")
    index_key = _text(index_name)
    if not index_key or "/" in index_key or index_key in {".", ".."}:
        raise InputPolicyError(f"{name}.{axis} dataframe index must name a direct dataset")
    index = _dataset(frame, index_key, f"{name}.{axis}")
    if index.ndim != 1:
        raise InputPolicyError(f"{name}.{axis} index must be one-dimensional")
    return index


def _strings(dataset: h5py.Dataset, rows: np.ndarray | None = None) -> np.ndarray:
    if dataset.ndim != 1:
        raise InputPolicyError("supported H5AD string array must be one-dimensional")
    if h5py.check_string_dtype(dataset.dtype) is None:
        raise InputPolicyError("supported H5AD string array must use a string dtype")
    values = dataset[:] if rows is None else dataset[rows]
    if np.asarray(values).dtype.fields is not None:
        raise InputPolicyError("structured H5AD string arrays are unsupported")
    return np.asarray([_text(value) for value in values], dtype=object)


def _labels(root: h5py.File, name: str, n_obs: int) -> np.ndarray:
    if n_obs > MAX_LABEL_VALUES:
        raise InputPolicyError(f"{name}.obs['celltype'] has {n_obs} values; ScoreLens metadata maximum is {MAX_LABEL_VALUES}")
    obs = _local_member(root, "obs", f"{name}.obs")
    if not isinstance(obs, h5py.Group) or "celltype" not in obs:
        raise InputPolicyError(f"{name}.obs['celltype'] is required for the T1 frozen probe")
    value = _local_member(obs, "celltype", f"{name}.obs['celltype']")
    if isinstance(value, h5py.Dataset):
        if value.ndim != 1:
            raise InputPolicyError(f"{name}.obs['celltype'] must be one-dimensional")
        if len(value) != n_obs:
            raise InputPolicyError(f"{name}.obs['celltype'] length does not match X rows")
        if h5py.check_string_dtype(value.dtype) is None:
            raise InputPolicyError(f"{name}.obs['celltype'] must use a supported string encoding")
        return np.asarray([_text(item) for item in value[:]], dtype=object)
    if not isinstance(value, h5py.Group) or _encoding(value) != "categorical":
        raise InputPolicyError(f"{name}.obs['celltype'] has unsupported H5AD encoding")
    codes, categories = (_dataset(value, key, f"{name}.obs['celltype']") for key in ("codes", "categories"))
    if codes.ndim != 1 or not np.issubdtype(codes.dtype, np.integer):
        raise InputPolicyError(f"{name}.obs['celltype'] categorical codes must be one-dimensional integers")
    if len(codes) != n_obs:
        raise InputPolicyError(f"{name}.obs['celltype'] categorical codes length does not match X rows")
    if categories.ndim != 1 or len(categories) > MAX_LABEL_VALUES:
        raise InputPolicyError(f"{name}.obs['celltype'] categories must be one-dimensional and within the metadata limit")
    cats, raw = _strings(categories), np.asarray(codes[:])
    if np.any((raw < -1) | (raw >= len(cats))):
        raise InputPolicyError(f"{name}.obs['celltype'] categorical codes are out of range")
    return np.asarray([None if code == -1 else cats[int(code)] for code in raw], dtype=object)


def _matrix(root: h5py.File, name: str) -> tuple[str, tuple[int, int], np.dtype]:
    x = _local_member(root, "X", f"{name}.X")
    if x is None:
        raise InputPolicyError(f"{name}.X is missing; a 2D cells x genes matrix is required")
    if isinstance(x, h5py.Dataset):
        if x.ndim != 2:
            raise InputPolicyError(f"{name}.X must be a 2D cells x genes matrix")
        return "dense", (int(x.shape[0]), int(x.shape[1])), _real_dtype(x.dtype, name)
    if not isinstance(x, h5py.Group) or _encoding(x) not in {"csr_matrix", "csc_matrix"}:
        raise InputPolicyError(f"{name}.X has unsupported H5AD encoding; supported: dense, csr_matrix, csc_matrix")
    shape = x.attrs.get("shape")
    if shape is None or np.asarray(shape).shape != (2,) or not np.issubdtype(np.asarray(shape).dtype, np.integer):
        raise InputPolicyError(f"{name}.X sparse shape must contain two dimensions")
    n_obs, n_vars = (int(v) for v in np.asarray(shape).reshape(-1))
    if n_obs < 1 or n_vars < 1:
        raise InputPolicyError(f"{name}.X sparse shape dimensions must be positive")
    data, indices, indptr = (_dataset(x, key, f"{name}.X") for key in ("data", "indices", "indptr"))
    expected = n_obs + 1 if _encoding(x) == "csr_matrix" else n_vars + 1
    if data.ndim != 1 or indices.ndim != 1 or indptr.ndim != 1 or len(data) != len(indices) or len(indptr) != expected:
        raise InputPolicyError(f"{name}.X sparse arrays must be compatible one-dimensional datasets")
    if not np.issubdtype(indices.dtype, np.integer) or not np.issubdtype(indptr.dtype, np.integer):
        raise InputPolicyError(f"{name}.X sparse index arrays must be integers")
    endpoints = np.asarray(indptr[[0, len(indptr) - 1]])
    if int(endpoints[0]) != 0 or int(endpoints[1]) != len(data):
        raise InputPolicyError(f"{name}.X sparse indptr endpoints are invalid")
    return _encoding(x), (n_obs, n_vars), _real_dtype(data.dtype, name)


def _inspect(path: Path, name: str, *, labels: bool, max_genes: int) -> _Meta:
    try:
        with h5py.File(path, "r") as root:
            layout, (n_obs, n_vars), dtype = _matrix(root, name)
            if n_obs < MIN_CELLS:
                raise InputPolicyError(f"{name} has {n_obs} cells; at least {MIN_CELLS} are required")
            var, obs = _index_dataset(root, "var", name), _index_dataset(root, "obs", name)
            if len(var) != n_vars or len(obs) != n_obs:
                raise InputPolicyError(f"{name} dataframe index length does not match X")
            if n_vars < MIN_GENES: raise InputPolicyError(f"{name} has {n_vars} genes; at least {MIN_GENES} are required")
            if n_vars > max_genes: raise InputPolicyError(f"{name} has {n_vars} genes; ScoreLens maximum is {max_genes}")
            genes = _strings(var).tolist()
            if len(set(genes)) != len(genes): raise InputPolicyError(f"{name}.var_names must be unique")
            target_labels = _labels(root, name, n_obs) if labels else None
            if target_labels is not None and len(target_labels) != n_obs: raise InputPolicyError(f"{name}.obs['celltype'] length does not match X rows")
            return _Meta(path, n_obs, n_vars, dtype, genes, layout, target_labels)
    except OSError as exc:
        raise InputPolicyError(f"{name} is not a readable H5AD/HDF5 input: {path}") from exc


def _select_rows(n_rows: int, count: int, seed: int, fingerprint: str) -> np.ndarray:
    seed_value = int.from_bytes(sha256(f"{seed}:{fingerprint}".encode()).digest()[:8], "little")
    return np.sort(np.random.default_rng(seed_value).choice(n_rows, size=count, replace=False)).astype(np.intp, copy=False)


def _estimated_peak_bytes(specs: list[tuple[int, int, np.dtype]]) -> int:
    return sum(rows * genes * (max(np.dtype(dtype).itemsize, 4) + 4) for rows, genes, dtype in specs)


def _csr_selected(group: h5py.Group, rows: np.ndarray, n_vars: int, name: str) -> sparse.csr_matrix:
    data, indices, indptr = (_dataset(group, key, f"{name}.X") for key in ("data", "indices", "indptr"))
    points = np.unique(np.concatenate((rows, rows + 1)))
    bounds = dict(zip(points.tolist(), np.asarray(indptr[points]).tolist()))
    values: list[np.ndarray] = []; cols: list[np.ndarray] = []; outptr = [0]
    for row in rows:
        start, stop = int(bounds[int(row)]), int(bounds[int(row) + 1])
        if start < 0 or stop < start or stop > len(data): raise InputPolicyError(f"{name}.X CSR indptr is invalid for selected rows")
        if stop - start > n_vars:
            raise InputPolicyError(f"{name}.X selected CSR row has more stored entries than genes; excessive duplicate storage is unsupported")
        col = np.asarray(indices[start:stop])
        if np.any((col < 0) | (col >= n_vars)): raise InputPolicyError(f"{name}.X CSR indices are out of range")
        selected_values = np.asarray(data[start:stop])
        if not np.isfinite(selected_values).all(): raise InputPolicyError(f"{name}.X selected rows contain NaN or infinite values")
        if (selected_values < 0).any(): raise InputPolicyError(f"{name}.X selected rows contain negative values")
        values.append(selected_values); cols.append(col); outptr.append(outptr[-1] + stop - start)
    return sparse.csr_matrix((np.concatenate(values) if values else [], np.concatenate(cols) if cols else [], np.asarray(outptr)), shape=(len(rows), n_vars))


def _csc_selected(group: h5py.Group, rows: np.ndarray, n_obs: int, n_vars: int, name: str) -> sparse.csr_matrix:
    data, indices, indptr = (_dataset(group, key, f"{name}.X") for key in ("data", "indices", "indptr"))
    pointers = np.asarray(indptr[:])
    if np.any(pointers[1:] < pointers[:-1]) or int(pointers[-1]) != len(data): raise InputPolicyError(f"{name}.X CSC indptr is invalid")
    output_rows: list[np.ndarray] = []; output_cols: list[np.ndarray] = []; output_data: list[np.ndarray] = []
    retained_nnz = 0
    for start in range(0, len(data), _CSC_CHUNK_NNZ):
        stop = min(start + _CSC_CHUNK_NNZ, len(data)); idx, values = np.asarray(indices[start:stop]), np.asarray(data[start:stop])
        if np.any((idx < 0) | (idx >= n_obs)): raise InputPolicyError(f"{name}.X CSC indices are out of range")
        position = np.searchsorted(rows, idx)
        keep = (position < len(rows)) & (rows[np.minimum(position, len(rows) - 1)] == idx)
        if np.any(keep):
            retained_nnz += int(keep.sum())
            if retained_nnz > len(rows) * n_vars:
                raise InputPolicyError(f"{name}.X selected CSC entries exceed the dense selected shape; excessive duplicate storage is unsupported")
            selected_values = values[keep]
            if not np.isfinite(selected_values).all(): raise InputPolicyError(f"{name}.X selected rows contain NaN or infinite values")
            if (selected_values < 0).any(): raise InputPolicyError(f"{name}.X selected rows contain negative values")
            output_rows.append(position[keep]); output_cols.append(np.searchsorted(pointers, np.arange(start, stop)[keep], side="right") - 1); output_data.append(values[keep])
    return sparse.coo_matrix((np.concatenate(output_data) if output_data else [], (np.concatenate(output_rows) if output_rows else [], np.concatenate(output_cols) if output_cols else [])), shape=(len(rows), n_vars)).tocsr()


def _read_selected(meta: _Meta, rows: np.ndarray, name: str) -> np.ndarray:
    with h5py.File(meta.path, "r") as root:
        x = _local_member(root, "X", f"{name}.X")
        if meta.layout == "dense": raw: Any = np.asarray(x[rows, :])
        elif meta.layout == "csr_matrix": raw = _csr_selected(x, rows, meta.n_vars, name).toarray()
        else: raw = _csc_selected(x, rows, meta.n_obs, meta.n_vars, name).toarray()
    raw = np.asarray(raw)
    if not np.isfinite(raw).all(): raise InputPolicyError(f"{name}.X selected rows contain NaN or infinite values")
    if (raw < 0).any(): raise InputPolicyError(f"{name}.X selected rows contain negative values")
    result = np.asarray(raw, dtype=np.float32)
    if not np.isfinite(result).all(): raise InputPolicyError(f"{name}.X selected rows are not representable as finite float32")
    result.setflags(write=False); return result


def _rows(meta: _Meta, selected: np.ndarray, name: str) -> list[dict[str, Any]]:
    with h5py.File(meta.path, "r") as root:
        names = _strings(_index_dataset(root, "obs", name), selected)
    return [{"position": int(row), "obs_name": str(obs_name)} for row, obs_name in zip(selected, names)]


def _coords(meta: _Meta, rows: np.ndarray, name: str) -> np.ndarray:
    with h5py.File(meta.path, "r") as root:
        obsm = _local_member(root, "obsm", f"{name}.obsm")
        value = _local_member(obsm, "spatial_3D", f"{name}.obsm['spatial_3D']") if isinstance(obsm, h5py.Group) else None
        if not isinstance(value, h5py.Dataset) or value.ndim != 2 or value.shape[0] != meta.n_obs or value.shape[1] < 3:
            raise InputPolicyError(f"{name}.obsm['spatial_3D'] must be row-aligned with at least three columns")
        if _encoding(value) not in {"array", ""} or not np.issubdtype(value.dtype, np.number) or np.issubdtype(value.dtype, np.complexfloating):
            raise InputPolicyError(f"{name}.obsm['spatial_3D'] has unsupported encoding or dtype")
        result = np.asarray(value[rows, :3], dtype=np.float32)
    if not np.isfinite(result).all(): raise InputPolicyError(f"{name}.obsm['spatial_3D'] selected rows contain NaN or infinite values")
    result.setflags(write=False); return result


def _valid_label(value: object) -> bool:
    if value is None: return False
    try: return bool(str(value).strip()) and not bool(np.asarray(value != value))
    except (TypeError, ValueError): return False


def load_task_bundle(*, task: str, prediction: str | Path, target: str | Path, reference: str | Path | None = None, wt: str | Path | None = None, max_cells: int = 256, seed: int = 0, max_dense_mib: int = MAX_DENSE_MIB) -> Bundle:
    task = str(task).upper()
    if task not in {"T1", "T2", "T3"}: raise InputPolicyError("task must be one of T1, T2, T3")
    if task == "T3":
        if reference is not None or wt is None: raise InputPolicyError("T3 requires wt and does not accept reference")
        third, contrast_role = wt, "wt"
    else:
        if reference is None or wt is not None: raise InputPolicyError(f"{task} requires reference and does not accept wt")
        third, contrast_role = reference, "reference"
    for option, value in (("max_cells", max_cells), ("max_dense_mib", max_dense_mib), ("seed", seed)):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.number)) or not np.isfinite(value):
            raise InputPolicyError(f"{option} must be a finite numeric value")
    if int(max_cells) != max_cells or not MIN_CELLS <= max_cells <= MAX_CELLS: raise InputPolicyError(f"max_cells must be an integer between {MIN_CELLS} and {MAX_CELLS}")
    if not 1 <= max_dense_mib <= MAX_DENSE_MIB: raise InputPolicyError(f"max_dense_mib must be between 1 and {MAX_DENSE_MIB}")
    if int(seed) != seed or seed < 0: raise InputPolicyError("seed must be a nonnegative integer")
    paths = {"prediction": Path(prediction), "target": Path(target), "reference": Path(third)}
    for name, path in paths.items():
        if not path.is_file(): raise FileNotFoundError(f"{name} file does not exist: {path}")
    fingerprints = {name: _fingerprint(path) for name, path in paths.items()}
    max_genes = MAX_GENES if task == "T1" else MAX_MULTITASK_GENES
    metas = {name: _inspect(path, name, labels=(name == "target"), max_genes=max_genes) for name, path in paths.items()}
    genes = metas["prediction"].genes
    for name in ("target", "reference"):
        if metas[name].genes != genes: raise InputPolicyError(f"{name}.var_names must exactly match prediction genes and order")
    selected = {name: _select_rows(meta.n_obs, min(int(max_cells), meta.n_obs), seed, fingerprints[name]["sha256"]) for name, meta in metas.items()}
    estimate = _estimated_peak_bytes([(len(selected[name]), len(genes), metas[name].dtype) for name in paths])
    if estimate > int(max_dense_mib) * _MIB: raise InputPolicyError(f"selected read/conversion peak estimate {estimate / _MIB:.1f} MiB exceeds max_dense_mib={max_dense_mib} MiB")
    arrays = {name: _read_selected(metas[name], selected[name], name) for name in paths}
    coords = {name: _coords(metas[name], selected[name], name) for name in paths} if task != "T1" else {name: None for name in paths}
    all_labels = metas["target"].labels; assert all_labels is not None
    valid = [str(value) for value in all_labels if _valid_label(value)]; full_classes = sorted(set(valid))
    if len(full_classes) < 2: raise InputPolicyError("target needs at least two nonmissing celltype classes")
    chosen = all_labels[selected["target"]]
    if not all(_valid_label(value) for value in chosen): raise InputPolicyError("target selected rows contain missing/blank celltype labels")
    labels = np.asarray([str(value) for value in chosen], dtype=str); selected_classes = sorted(set(labels.tolist()))
    if len(selected_classes) < 2: raise InputPolicyError("uniform target sampling retained fewer than two celltype classes; increase max_cells or choose a different seed")
    labels.setflags(write=False); require_target_variation(arrays["target"])
    lost = sorted(set(full_classes) - set(selected_classes)); warnings = [] if not lost else ["Uniform target sampling omitted valid full-file celltype classes: " + ", ".join(lost)]
    label_provenance: dict[str, Any] = {"full_valid_classes": full_classes, "selected_classes": selected_classes, "full_missing_or_blank_count": int(len(all_labels) - len(valid))}
    provenance_names = {"reference": "wt"} if task == "T3" else {}
    public_name = lambda name: provenance_names.get(name, name)
    coordinate_provenance = {} if task == "T1" else {public_name(name): {"source": "obsm/spatial_3D", "selected_shape": list(coords[name].shape)} for name in paths}
    provenance: dict[str, Any] = {"scope": "sampled_subset_only", "scope_note": "Only selected rows were numerically scanned and scored; no full-file equivalence is claimed.", "task": task, "contrast_role": contrast_role, "selection": {"method": "uniform_without_replacement_sorted_positions", "requested_max_cells": int(max_cells), "seed": int(seed), "rows": {public_name(name): _rows(metas[name], selected[name], name) for name in paths}}, "inputs": {public_name(name): {"fingerprint": fingerprints[name], "original_shape": [metas[name].n_obs, metas[name].n_vars], "dtype": str(metas[name].dtype), "layout": metas[name].layout, "selected_shape": list(arrays[name].shape)} for name in paths}, "coordinates": coordinate_provenance, "genes": {"count": len(genes), "ordered_sha256": sha256("\0".join(genes).encode()).hexdigest()}, "target_celltypes": label_provenance, "resource": {"max_dense_mib": int(max_dense_mib), "estimated_read_conversion_peak_mib": round(estimate / _MIB, 4), "selected_float32_arrays_mib": round(sum(value.nbytes for value in arrays.values()) / _MIB, 4), "selected_coordinate_arrays_mib": round(sum(value.nbytes for value in coords.values() if value is not None) / _MIB, 4)}, "warnings": warnings}
    return Bundle(arrays["prediction"], arrays["target"], arrays["reference"], list(genes), labels, provenance, task, coords["prediction"], coords["target"], coords["reference"])


def load_bundle(prediction: str | Path, target: str | Path, reference: str | Path, *, max_cells: int = 256, seed: int = 0, max_dense_mib: int = MAX_DENSE_MIB) -> Bundle:
    """Backward-compatible T1 wrapper around the selective task loader."""
    return load_task_bundle(task="T1", prediction=prediction, target=target, reference=reference, max_cells=max_cells, seed=seed, max_dense_mib=max_dense_mib)
