"""A reported input SHA-256 must cover every consumed HDF5 value."""
from __future__ import annotations

import h5py
import numpy as np
import pytest

from vec_scorelens.io import InputPolicyError, load_task_bundle


STRING = h5py.string_dtype("utf-8")


def _input(path, *, target=False, coords=False):
    with h5py.File(path, "w") as f:
        x = np.ones((40, 4), dtype="f4")
        x[:, 0] += np.arange(40) % 5
        f.create_dataset("X", data=x)
        for axis, values in (("obs", [f"cell_{i}" for i in range(40)]), ("var", list("abcd"))):
            frame = f.create_group(axis)
            frame.attrs["_index"] = "_index"
            frame.create_dataset("_index", data=np.asarray(values, dtype=object), dtype=STRING)
        if target:
            f["obs"].create_dataset("celltype", data=np.asarray(["A" if i % 2 else "B" for i in range(40)], dtype=object), dtype=STRING)
        if coords:
            f.create_group("obsm").create_dataset("spatial_3D", data=np.tile([1., 2., 3.], (40, 1)))


@pytest.fixture
def files(tmp_path):
    paths = [tmp_path / f"{role}.h5ad" for role in ("prediction", "target", "reference")]
    for index, path in enumerate(paths):
        _input(path, target=index == 1, coords=True)
    return paths


def _load(paths, task="T1"):
    return load_task_bundle(task=task, prediction=paths[0], target=paths[1], reference=paths[2])


@pytest.mark.parametrize("location,role,task", [
    ("X", 0, "T1"), ("obs/_index", 0, "T1"), ("var/_index", 0, "T1"),
    ("obs/celltype", 1, "T1"), ("obsm/spatial_3D", 0, "T2"),
])
def test_external_link_on_consumed_path_is_rejected(files, tmp_path, location, role, task):
    sidecar = tmp_path / "sidecar.h5"
    with h5py.File(files[role], "r") as source, h5py.File(sidecar, "w") as external:
        source.copy(location, external, name="value")
    with h5py.File(files[role], "r+") as source:
        del source[location]
        source[location] = h5py.ExternalLink(sidecar.name, "/value")
    with pytest.raises(InputPolicyError, match="external and soft HDF5 links"):
        _load(files, task)


def test_soft_link_to_external_link_is_rejected(files, tmp_path):
    sidecar = tmp_path / "sidecar.h5"
    with h5py.File(sidecar, "w") as external:
        external.create_dataset("X", data=np.ones((40, 4), dtype="f4"))
    with h5py.File(files[0], "r+") as source:
        del source["X"]
        source["external"] = h5py.ExternalLink(sidecar.name, "/X")
        source["X"] = h5py.SoftLink("/external")
    with pytest.raises(InputPolicyError, match="external and soft HDF5 links"):
        _load(files)


def test_virtual_dataset_is_rejected(files, tmp_path):
    sidecar = tmp_path / "sidecar.h5"
    with h5py.File(sidecar, "w") as external:
        external.create_dataset("X", data=np.ones((40, 4), dtype="f4"))
    layout = h5py.VirtualLayout(shape=(40, 4), dtype="f4")
    layout[:, :] = h5py.VirtualSource(str(sidecar), "X", shape=(40, 4))
    with h5py.File(files[0], "r+") as source:
        del source["X"]
        source.create_virtual_dataset("X", layout)
    with pytest.raises(InputPolicyError, match="external or virtual dataset storage"):
        _load(files)


def test_external_raw_storage_is_rejected(files, tmp_path):
    backing = tmp_path / "backing.bin"
    with h5py.File(files[0], "r+") as source:
        del source["X"]
        source.create_dataset("X", shape=(40, 4), dtype="f4", external=[(str(backing), 0, 40 * 4 * 4)])
    with pytest.raises(InputPolicyError, match="external or virtual dataset storage"):
        _load(files)


def test_dataframe_index_cannot_traverse_external_intermediate(files, tmp_path):
    sidecar = tmp_path / "sidecar.h5"
    with h5py.File(sidecar, "w") as external:
        external.create_dataset("index", data=np.asarray([f"cell_{i}" for i in range(40)], dtype=object), dtype=STRING)
    with h5py.File(files[0], "r+") as source:
        source["obs"]["external"] = h5py.ExternalLink(sidecar.name, "/")
        source["obs"].attrs["_index"] = "external/index"
    with pytest.raises(InputPolicyError, match="direct dataset"):
        _load(files)
