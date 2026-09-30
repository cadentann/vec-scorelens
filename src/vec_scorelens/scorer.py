"""Source-pinned execution of upstream T1/T2/T3 v2 panels on bounded inputs."""
from __future__ import annotations

from hashlib import sha256
from importlib import import_module
from importlib.metadata import PackageNotFoundError, distributions
import json
from pathlib import Path
import subprocess
import sys
from types import ModuleType
from typing import Any
from urllib.parse import unquote, urlparse

import numpy as np

from .io import Bundle, InputPolicyError, require_target_variation


PINNED_COMMIT = "46d41e63f42a9aab815db20b742feeccd249cb17"
PINNED_VERSION = "0.1.1"
_PINNED_HASHES = {
    "score_h5ad.py": "52034554f03aec10193cb09baa2c78a1de04218ef678f327eb63acc58fa28add",
    "veckit.py": "b991dc27e1b169b732721f4d79ff183c29879f6c4fee9ceb7e897a63172b4108",
    "T1/metrics.py": "c6252f964e782b4bb20f8b66352eb43ee244de402dafea5d3a0864ad5ed2fab8",
    "T1/metrics_v2.py": "8d49a59cb853da0763e48236a07f2b141e485437314e69812f27bf8aa53ba313",
    "T2/metrics.py": "2c8aeb720aa44661d0d05bfd6d389093b8d879601a41ac7280a1da506a64bd13",
    "T2/metrics_v2.py": "5d06e68f6d6848f64eb7b798ea42ca31d8c1f9a61084c174a4f50eaf88bdd878",
    "T3/metrics.py": "b338d3ff6e977f093922ec5de2dd05c07faa472b29d44d901ce9a7cb432d1376",
    "T3/metrics_v2.py": "ebc0882a1c1d3a83f05b577ef0084e5fb2e9075faad065304d094df1727f67a2",
    "common/core_metrics.py": "3be7099a0c9a7ad5b609f86078871ed8290ed0bd8916b0cb3e14fb9831909f31",
    "common/shape_metrics.py": "bd80ef2705f9be93fe1538e77445d4eb217d4b955c6cd38c69b26b0671be91c0",
}


class SourcePinError(RuntimeError):
    """The imported veckit source does not match ScoreLens' accepted scorer."""


def _file_hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _module_path(module: ModuleType, label: str) -> Path:
    module_file = getattr(module, "__file__", None)
    if not module_file:
        raise SourcePinError(f"{label} has no importable source file")
    return Path(module_file).resolve()


def _require_module_path(module: ModuleType, expected: Path, label: str) -> None:
    actual = _module_path(module, label)
    if actual != expected.resolve():
        raise SourcePinError(f"{label} was imported from {actual}, expected {expected.resolve()}")


def _direct_url_metadata(dist: Any) -> dict[str, Any]:
    raw = dist.read_text("direct_url.json")
    if not raw:
        return {}
    try:
        metadata = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise SourcePinError("veckit direct_url.json is not valid JSON") from exc
    if not isinstance(metadata, dict):
        raise SourcePinError("veckit direct_url.json must contain an object")
    return metadata


def _direct_url_commit(metadata: dict[str, Any]) -> str | None:
    commit = metadata.get("vcs_info", {}).get("commit_id")
    return str(commit) if commit else None


def _editable_root(metadata: dict[str, Any]) -> Path | None:
    if not metadata.get("dir_info", {}).get("editable"):
        return None
    url = metadata.get("url")
    if not isinstance(url, str):
        return None
    parsed = urlparse(url)
    if parsed.scheme != "file":
        return None
    return Path(unquote(parsed.path)).resolve()


def _runtime_distribution(score_path: Path) -> tuple[Any, dict[str, Any]]:
    """Choose metadata that attests to the *currently imported* scorer file.

    Upstream's dynamic loader inserts its source root into ``sys.path``.  A
    later plain ``distribution('veckit')`` can then select a source-side
    ``.egg-info`` without PEP 610 provenance instead of the installed editable
    distribution.  Match metadata to the loaded runtime, not path precedence.
    """
    matches: list[tuple[Any, dict[str, Any]]] = []
    for candidate in distributions(name="veckit"):
        metadata = _direct_url_metadata(candidate)
        candidate_score_path = Path(candidate.locate_file("score_h5ad.py")).resolve()
        if candidate_score_path.is_file() and candidate_score_path == score_path:
            matches.append((candidate, metadata))
        elif _editable_root(metadata) == score_path.parent:
            matches.append((candidate, metadata))
    if not matches:
        raise SourcePinError(f"no installed veckit metadata attests to imported {score_path}")
    # Prefer a PEP 610 record: it carries the VCS commit used by a wheel and
    # remains valid after the scorer has mutated sys.path.
    matches.sort(key=lambda item: bool(item[1]), reverse=True)
    return matches[0]


def verify_source_pin() -> dict[str, Any]:
    """Fail closed unless the imported scorer is exactly the reviewed source."""
    # Import and bind every runtime module that carries scorer semantics.  The
    # distribution path check prevents an older globally-importable module from
    # satisfying hashes rooted in a different installation.
    try:
        veckit_module = import_module("veckit")
        score_module = import_module("score_h5ad")
        core_module = import_module("common.core_metrics")
        shape_module = import_module("common.shape_metrics")
        legacy_module = import_module("T1.metrics")
    except ImportError as exc:
        raise SourcePinError("veckit runtime modules are unavailable") from exc
    score_path = _module_path(score_module, "score_h5ad")
    root = score_path.parent
    try:
        dist, direct_url = _runtime_distribution(score_path)
    except PackageNotFoundError as exc:
        raise SourcePinError("veckit package metadata is unavailable") from exc
    installed_version = dist.version
    if installed_version != PINNED_VERSION:
        raise SourcePinError(f"expected veckit {PINNED_VERSION}, found {installed_version}")
    distribution_score_path = Path(dist.locate_file("score_h5ad.py")).resolve()
    if distribution_score_path.is_file() and score_path != distribution_score_path:
        raise SourcePinError(
            f"score_h5ad was imported from {score_path}, not installed veckit distribution {distribution_score_path}"
        )
    if not distribution_score_path.is_file():
        editable_root = _editable_root(direct_url)
        if editable_root is None or root != editable_root:
            raise SourcePinError(
                f"score_h5ad was imported from {score_path}; installed veckit metadata does not attest that root"
            )
    _require_module_path(veckit_module, root / "veckit.py", "veckit")
    _require_module_path(core_module, root / "common" / "core_metrics.py", "common.core_metrics")
    _require_module_path(shape_module, root / "common" / "shape_metrics.py", "common.shape_metrics")
    _require_module_path(legacy_module, root / "T1" / "metrics.py", "T1.metrics")
    missing = [relative for relative in _PINNED_HASHES if not (root / relative).is_file()]
    if missing:
        raise SourcePinError("veckit source root is incomplete: " + ", ".join(missing))

    # A regular wheel can live underneath an unrelated project's .git (for
    # example `.venv/site-packages`).  Prefer its PEP 610 VCS provenance and
    # never let Git's parent-directory discovery attest that wheel.  Git is an
    # allowed fallback only for an editable install rooted at the actual
    # checked-out scorer.
    direct_commit = _direct_url_commit(direct_url)
    if direct_commit:
        commit = direct_commit
        commit_provenance = "direct_url.vcs_info.commit_id"
    else:
        editable_root = _editable_root(direct_url)
        if editable_root is None or editable_root != root:
            raise SourcePinError("veckit needs VCS direct_url provenance or a verified editable Git root")
        top_level = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "--show-toplevel"],
            check=False,
            capture_output=True,
            text=True,
        )
        if top_level.returncode != 0 or Path(top_level.stdout.strip()).resolve() != root:
            raise SourcePinError("editable veckit root is not its own Git worktree")
        head = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            check=False,
            capture_output=True,
            text=True,
        )
        commit = head.stdout.strip()
        if head.returncode != 0:
            raise SourcePinError("editable veckit Git commit is unavailable")
        commit_provenance = "git_worktree"
    if commit != PINNED_COMMIT:
        raise SourcePinError(f"expected veckit commit {PINNED_COMMIT}, found {commit or 'unavailable'}")
    actual_hashes = {relative: _file_hash(root / relative) for relative in _PINNED_HASHES}
    mismatched = [relative for relative, expected in _PINNED_HASHES.items() if actual_hashes[relative] != expected]
    if mismatched:
        raise SourcePinError("pinned scorer hash mismatch: " + ", ".join(mismatched))
    return {
        "package": "veckit",
        "version": installed_version,
        "commit": commit,
        "commit_provenance": commit_provenance,
        "source_root": str(root),
        "source_hashes": actual_hashes,
    }


def _loaded_task_modules(task: str, root: Path, legacy: ModuleType, panel: ModuleType) -> ModuleType | None:
    """Attest every module loaded through upstream's clean reuse stand-in.

    The pinned loader intentionally never imports a real ``_reuse.py`` or a
    task's ``data.py``. Its stand-in must expose only the pinned dependencies.
    Check the bound functions too: changing tasks replaces generic module names
    in sys.modules, but each panel must retain the right task's function objects.
    """
    _require_module_path(legacy, root / task / "metrics.py", f"loaded {task} metrics")
    _require_module_path(panel, root / task / "metrics_v2.py", f"loaded {task} metrics_v2")
    if sys.modules.get("metrics") is not legacy or sys.modules.get("metrics_v2") is not panel:
        raise SourcePinError("upstream task loader did not bind the selected panel")
    reuse = sys.modules.get("_reuse")
    if not isinstance(reuse, ModuleType) or getattr(reuse, "__file__", None):
        raise SourcePinError("upstream task loader needs its clean in-memory _reuse stand-in")
    expected_reuse = {"t1_metrics"} | ({"t2_metrics"} if task != "T1" else set())
    public_reuse = {key for key in vars(reuse) if not key.startswith("__")}
    if public_reuse != expected_reuse:
        raise SourcePinError("upstream _reuse dependencies do not match the selected task")
    t1 = reuse.t1_metrics
    _require_module_path(t1, root / "T1" / "metrics.py", "_reuse.t1_metrics")
    if sys.modules.get("t1_metrics") is not t1:
        raise SourcePinError("upstream _reuse.t1_metrics binding is inconsistent")
    spatial_legacy = None
    if task != "T1":
        spatial_legacy = reuse.t2_metrics
        _require_module_path(spatial_legacy, root / "T2" / "metrics.py", "_reuse.t2_metrics")
        if sys.modules.get("t2_metrics") is not spatial_legacy:
            raise SourcePinError("upstream _reuse.t2_metrics binding is inconsistent")
        if spatial_legacy.train_frozen_probe is not t1.train_frozen_probe:
            raise SourcePinError("upstream T2 reuse binding is inconsistent")
        # T2 loads its source twice (once for the reuse stand-in, once as
        # ``metrics``), so its two neighborhood functions are distinct objects.
        if task == "T3" and legacy.neighborhood_mmd is not spatial_legacy.neighborhood_mmd:
            raise SourcePinError("upstream spatial reuse binding is inconsistent")
    if panel.train_frozen_probe is not legacy.train_frozen_probe or panel.composition_jsd is not legacy.composition_jsd:
        raise SourcePinError("upstream panel is bound to a different task's legacy metrics")
    return spatial_legacy


def _require_spatial_support(bundle: Bundle, spatial_legacy: ModuleType) -> None:
    for role in ("prediction", "target", "reference"):
        coords = getattr(bundle, f"{role}_coords", None)
        expression = getattr(bundle, role)
        if coords is None or coords.ndim != 2 or coords.shape != (expression.shape[0], 3):
            raise InputPolicyError(f"{role} needs selected spatial_3D coordinates with shape cells x 3")
        if not np.isfinite(coords).all():
            raise InputPolicyError(f"{role} selected spatial_3D coordinates contain NaN or infinite values")
        if role != "reference":
            centered = np.asarray(coords, dtype=np.float64) - np.asarray(coords, dtype=np.float64).mean(0)
            if float(np.sqrt(np.mean(np.sum(centered ** 2, axis=1)))) < 1e-12:
                raise InputPolicyError(f"{role} selected spatial_3D coordinates have zero spatial extent")
    # Upstream's neighborhood MMD fits its bandwidth to the target neighborhood
    # means. Varying raw rows alone is insufficient when those means collapse.
    neighborhoods = spatial_legacy._knn_neighborhood_pb(bundle.target, bundle.target_coords)
    if not np.any(neighborhoods != neighborhoods[0]):
        raise InputPolicyError(
            "target selected neighborhood pseudobulks are all identical; pinned neighborhood MMD requires within-target variation"
        )


def run_official(bundle: Bundle, seed: int = 0) -> dict[str, Any]:
    """Run the actual upstream task-v2 functions on immutable selected arrays.

    This deliberately bypasses file loading only because ``Bundle`` already
    records a bounded selected subset.  It does not reimplement any metric.
    """
    task = bundle.task
    if task not in ("T1", "T2", "T3"):
        raise InputPolicyError(f"unsupported task {task!r}; choose T1, T2, or T3")
    require_target_variation(bundle.target)
    source_provenance = verify_source_pin()
    score_h5ad = import_module("score_h5ad")

    legacy, metrics_v2 = score_h5ad._load_task_metrics(task)
    root = Path(source_provenance["source_root"])
    spatial_legacy = _loaded_task_modules(task, root, legacy, metrics_v2)
    core = import_module("common.core_metrics")
    _require_module_path(core, root / "common" / "core_metrics.py", "loaded common.core_metrics")
    shape = import_module("common.shape_metrics")
    _require_module_path(shape, root / "common" / "shape_metrics.py", "loaded common.shape_metrics")
    if metrics_v2.mmd_unbiased is not core.mmd_unbiased:
        raise SourcePinError("upstream panel is bound to a different common.core_metrics")
    if task != "T1":
        if metrics_v2.d2_distance is not shape.d2_distance:
            raise SourcePinError("upstream panel is bound to a different common.shape_metrics")
        _require_spatial_support(bundle, spatial_legacy)
    probe = legacy.train_frozen_probe(bundle.target, bundle.target_labels)
    prediction_labels = np.full(bundle.prediction.shape[0], "NA", dtype=str)
    if task == "T1":
        metrics = metrics_v2.score_task1_v2(
            bundle.prediction, prediction_labels, bundle.target, bundle.target_labels,
            bundle.reference, probe=probe, seed=seed,
        )
    else:
        score = metrics_v2.score_task2_v2 if task == "T2" else metrics_v2.score_task3_v2
        metrics = score(
            bundle.prediction, bundle.prediction_coords, prediction_labels,
            bundle.target, bundle.target_coords, bundle.target_labels,
            bundle.reference, probe=probe, seed=seed,
        )
    return {
        "metrics": metrics,
        "probe": probe,
        "core": core,
        "legacy": legacy,
        "task": task,
        "shape": shape,
        "spatial_legacy": spatial_legacy,
        "provenance": {
            **source_provenance,
            "task": task,
            "seed": int(seed),
            "scope": bundle.provenance.get("scope"),
            "scope_note": "Official v2 formulas were run exactly on the bounded selected Bundle, not on full input files.",
        },
    }
