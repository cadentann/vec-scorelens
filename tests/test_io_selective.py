"""Synthetic acceptance tests for the narrow direct-HDF5 loader."""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pytest
from scipy import sparse

from vec_scorelens.io import InputPolicyError, load_bundle, load_task_bundle


def _write(path: Path, layout: str, *, coords: bool = False) -> None:
    values = (np.arange(32 * 4, dtype=np.float32).reshape(32, 4) / 100) + .1
    matrix = {"dense": values, "csr": sparse.csr_matrix(values), "csc": sparse.csc_matrix(values)}[layout]
    data = ad.AnnData(matrix)
    data.var_names = ["g0", "g1", "g2", "g3"]
    data.obs_names = [f"cell-{i}" for i in range(32)]
    data.obs["celltype"] = np.where(np.arange(32) % 2, "b", "a")
    if coords:
        data.obsm["spatial_3D"] = np.column_stack((np.arange(32), np.arange(32) + 1, np.arange(32) + 2, np.arange(32) + 3)).astype(np.float32)
    data.write_h5ad(path)


@pytest.mark.parametrize("layout", ["dense", "csr", "csc"])
def test_selective_loader_supports_all_x_layouts(tmp_path: Path, layout: str) -> None:
    paths = [tmp_path / f"{name}.h5ad" for name in ("prediction", "target", "reference")]
    for path in paths: _write(path, layout)
    bundle = load_bundle(*paths, max_cells=32, seed=17)
    assert bundle.task == "T1"
    assert bundle.prediction.shape == (32, 4)
    assert bundle.provenance["inputs"]["prediction"]["layout"] == ("dense" if layout == "dense" else f"{layout}_matrix")
    assert not bundle.prediction.flags.writeable


def test_loader_never_opens_unused_ann_data_sidecars(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    paths = [tmp_path / f"{name}.h5ad" for name in ("prediction", "target", "reference")]
    for path in paths:
        _write(path, "csr")
        with h5py.File(path, "r+") as handle:
            for group, key, shape in (("uns", "unused_payload", (4_000_000,)), ("obsm", "unused_embedding", (32, 250_000)), ("layers", "unused_layer", (32, 4))):
                value = handle[group].create_dataset(key, shape=shape, dtype="f4", chunks=True, fillvalue=0)
                value.attrs["encoding-type"] = "array"; value.attrs["encoding-version"] = "0.2.0"
    original = h5py.Dataset.__getitem__
    def guarded(dataset, key):
        if any(part in dataset.name for part in ("/uns/", "/obsm/unused", "/layers/unused", "/raw/")):
            raise AssertionError(f"unused sidecar was read: {dataset.name}")
        return original(dataset, key)
    monkeypatch.setattr(h5py.Dataset, "__getitem__", guarded)
    assert load_bundle(*paths, max_cells=32).prediction.shape == (32, 4)


def test_multitask_requires_row_aligned_spatial_and_t3_wt_role(tmp_path: Path) -> None:
    pred, target, wt = (tmp_path / f"{name}.h5ad" for name in ("prediction", "target", "wt"))
    for path in (pred, target, wt): _write(path, "dense", coords=True)
    bundle = load_task_bundle(task="T3", prediction=pred, target=target, wt=wt, max_cells=32)
    assert bundle.task == "T3"
    assert bundle.reference_coords is not None and bundle.reference_coords.shape == (32, 3)
    assert bundle.provenance["contrast_role"] == "wt"
    assert "wt" in bundle.provenance["inputs"] and "reference" not in bundle.provenance["inputs"]
    with h5py.File(pred, "r+") as handle:
        del handle["obsm/spatial_3D"]
        value = handle["obsm"].create_dataset("spatial_3D", shape=(31, 3), dtype="f4")
        value.attrs["encoding-type"] = "array"; value.attrs["encoding-version"] = "0.2.0"
    with pytest.raises(InputPolicyError, match="row-aligned"):
        load_task_bundle(task="T3", prediction=pred, target=target, wt=wt, max_cells=32)


def test_task_roles_and_multitask_gene_cap_are_explicit(tmp_path: Path) -> None:
    pred, target, ref = (tmp_path / f"{name}.h5ad" for name in ("prediction", "target", "reference"))
    for path in (pred, target, ref): _write(path, "dense", coords=True)
    with pytest.raises(InputPolicyError, match="requires wt"):
        load_task_bundle(task="T3", prediction=pred, target=target, reference=ref)
    with pytest.raises(InputPolicyError, match="does not accept wt"):
        load_task_bundle(task="T2", prediction=pred, target=target, reference=ref, wt=ref)
    with h5py.File(target, "r+") as handle:
        del handle["obs/celltype"]
        handle["obs"].create_dataset("celltype", data=np.arange(32, dtype=np.float32))
    with pytest.raises(InputPolicyError, match="supported string encoding"):
        load_task_bundle(task="T2", prediction=pred, target=target, reference=ref)


@pytest.mark.parametrize("categorical", [False, True])
def test_malformed_label_length_rejected_before_payload_read(tmp_path, monkeypatch, categorical):
    paths = [tmp_path / f"{name}.h5ad" for name in ("prediction", "target", "reference")]
    for path in paths:
        _write(path, "dense")
    with h5py.File(paths[1], "r+") as handle:
        del handle["obs/celltype"]
        if categorical:
            field = handle["obs"].create_group("celltype")
            field.attrs["encoding-type"] = "categorical"
            field.create_dataset("codes", shape=(10_000_000,), dtype="i4", chunks=True)
            field.create_dataset("categories", data=np.asarray(["a"], dtype=h5py.string_dtype()))
        else:
            handle["obs"].create_dataset("celltype", shape=(10_000_000,), dtype="S1", chunks=True)
    original = h5py.Dataset.__getitem__
    def guarded(dataset, key):
        if dataset.name.startswith("/obs/celltype"):
            raise AssertionError("malformed label payload was read before length validation")
        return original(dataset, key)
    monkeypatch.setattr(h5py.Dataset, "__getitem__", guarded)
    with pytest.raises(InputPolicyError, match="length does not match"):
        load_bundle(*paths, max_cells=32)
