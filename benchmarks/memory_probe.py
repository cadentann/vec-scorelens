#!/usr/bin/env python3
"""Portable bounded-memory acceptance probe for ScoreLens.

The probe writes three *synthetic-only* sparse H5AD inputs with a large logical
shape, then runs the normal ``load_bundle -> run_official -> diagnose ->
write_report`` path in a fresh subprocess.  It deliberately never needs an
input data file.  The worker guards sparse ``toarray`` calls so a regression
that densifies an unselected backing matrix fails rather than merely consuming
more memory. The default 2 GiB full-workflow test budget allows observed
process-RSS variability for 128 selected cells; it is not a total-memory
guarantee for other inputs, cell caps, or environments.

Examples (run from ``03_scorelens``)::

    python benchmarks/memory_probe.py
    python benchmarks/memory_probe.py --selected-cells 64 --result memory-result.json

``--result`` is optional and writes only portable, path-free evidence.  The
temporary H5ADs and private diagnostic report are removed after each run.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from typing import Any, Iterator

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

try:  # ``resource`` is unavailable on Windows; psutil covers that case.
    import resource
except ImportError:  # pragma: no cover - exercised on Windows only
    resource = None  # type: ignore[assignment]

try:  # Optional because POSIX can obtain a child peak from ``resource``.
    import psutil
except ImportError:  # pragma: no cover - normal on minimal POSIX installs
    psutil = None  # type: ignore[assignment]


LOGICAL_ROWS = 10_000
LOGICAL_GENES = 32_285
DENSITY = 0.001
DEFAULT_SELECTED_CELLS = 128
DEFAULT_MEMORY_LIMIT_BYTES = 2 << 30
DEFAULT_RUNTIME_LIMIT_SECONDS = 180.0
SEED = 61_903


@dataclass(frozen=True)
class Fingerprint:
    """Content and metadata evidence for one immutable synthetic source."""

    sha256: str
    bytes: int
    mtime_ns: int


def _fingerprint(path: Path) -> Fingerprint:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    stat = path.stat()
    return Fingerprint(digest.hexdigest(), int(stat.st_size), int(stat.st_mtime_ns))


def _sparse_matrix(*, rows: int, genes: int, density: float, seed: int) -> sparse.csr_matrix:
    """Make a deterministic nonnegative CSR matrix without a dense logical array."""
    nnz = int(rows * genes * density)
    rng = np.random.default_rng(seed)
    # Generator.choice has an O(nnz) implementation for this sparse sample;
    # it avoids allocating an array of ``rows * genes`` logical positions.
    positions = rng.choice(rows * genes, size=nnz, replace=False)
    values = rng.uniform(0.05, 3.0, size=nnz).astype(np.float32)
    matrix = sparse.csr_matrix(
        (values, (positions // genes, positions % genes)), shape=(rows, genes), dtype=np.float32
    )
    matrix.sort_indices()
    assert matrix.nnz == nnz
    return matrix


def write_sparse_inputs(directory: Path, *, rows: int, genes: int, density: float) -> dict[str, Path]:
    """Create sparse, non-biological inputs with two valid target classes."""
    if not directory.is_dir():
        raise ValueError("input directory must already exist")
    names = [f"synthetic_gene_{index:05d}" for index in range(genes)]
    var = pd.DataFrame(index=names)
    paths: dict[str, Path] = {}
    for offset, name in enumerate(("prediction", "target", "reference")):
        obs = pd.DataFrame(index=[f"synthetic_{name}_{index:05d}" for index in range(rows)])
        if name == "target":
            # Alternation makes every uniform sample of two or more adjacent
            # rows representative, while retaining a genuinely two-class probe.
            obs["celltype"] = pd.Categorical(np.where(np.arange(rows) % 2, "type_b", "type_a"))
        adata = ad.AnnData(
            X=_sparse_matrix(rows=rows, genes=genes, density=density, seed=SEED + offset),
            obs=obs,
            var=var.copy(),
        )
        adata.uns["synthetic_only"] = True
        adata.uns["memory_probe"] = {
            "logical_shape": [rows, genes], "density": density, "seed": SEED + offset,
        }
        path = directory / f"{name}.h5ad"
        adata.write_h5ad(path)
        paths[name] = path
    return paths


@contextmanager
def _forbid_large_sparse_toarray(selected_cap: int) -> Iterator[None]:
    """Reject accidental backing-matrix densification above the selected cap."""
    originals: dict[type, Any] = {}

    def guarded(original: Any):
        def toarray(self: sparse.spmatrix, *args: Any, **kwargs: Any) -> np.ndarray:
            if self.shape[0] > selected_cap:
                raise AssertionError(
                    f"refusing sparse .toarray() for {self.shape}; selected cap is {selected_cap} rows"
                )
            return original(self, *args, **kwargs)

        return toarray

    # H5AD backed slicing normally yields CSR, but cover common sparse formats
    # so a storage-format change remains protected by the acceptance probe.
    classes = (sparse.csr_matrix, sparse.csc_matrix, sparse.coo_matrix, sparse.lil_matrix, sparse.dok_matrix)
    try:
        for cls in classes:
            originals[cls] = cls.toarray
            cls.toarray = guarded(cls.toarray)
        yield
    finally:
        for cls, original in originals.items():
            cls.toarray = original


def _worker(args: argparse.Namespace) -> dict[str, Any]:
    """Run the actual public workflow and return portable evidence."""
    from vec_scorelens.diagnostics import diagnose
    from vec_scorelens.io import load_bundle
    from vec_scorelens.report import write_report
    from vec_scorelens.scorer import run_official

    work = Path(args.work).resolve()
    sources = {name: work / f"{name}.h5ad" for name in ("prediction", "target", "reference")}
    before = {name: _fingerprint(path) for name, path in sources.items()}
    report = work / "private_report"
    started = time.perf_counter()
    with _forbid_large_sparse_toarray(args.selected_cells):
        bundle = load_bundle(
            sources["prediction"], sources["target"], sources["reference"],
            max_cells=args.selected_cells, seed=SEED, max_dense_mib=args.max_dense_mib,
        )
        scored = run_official(bundle, seed=SEED)
        diagnostics = diagnose(bundle, scored, seed=SEED)
        final_report = write_report(bundle, scored, diagnostics, report)
    elapsed = time.perf_counter() - started
    after = {name: _fingerprint(path) for name, path in sources.items()}
    if before != after:
        raise AssertionError("source H5AD fingerprint, size, or mtime changed during ScoreLens execution")
    if not final_report.is_dir():
        raise AssertionError("full report workflow did not create its report directory")
    selected_shapes = {
        name: list(bundle.provenance["inputs"][name]["selected_shape"])
        for name in ("prediction", "target", "reference")
    }
    if any(shape[0] != args.selected_cells or shape[1] != args.genes for shape in selected_shapes.values()):
        raise AssertionError(f"unexpected selected shapes: {selected_shapes}")
    if len(bundle.provenance["target_celltypes"]["selected_classes"]) < 2:
        raise AssertionError("selected target rows did not retain two valid cell types")
    return {
        "worker_elapsed_seconds": elapsed,
        "worker_peak_rss_bytes": _resource_peak_rss_bytes(),
        "logical_shape": [args.rows, args.genes],
        "density": args.density,
        "declared_selected_shapes": selected_shapes,
        "selected_target_classes": bundle.provenance["target_celltypes"]["selected_classes"],
        "source_before": {name: asdict(value) for name, value in before.items()},
        "source_after": {name: asdict(value) for name, value in after.items()},
        "source_unchanged": True,
        "large_sparse_toarray_guard": "passed",
        "report_files": sorted(path.name for path in final_report.iterdir()),
    }


def _resource_peak_rss_bytes() -> int:
    """Normalize ru_maxrss across macOS (bytes) and Linux (KiB)."""
    if resource is None:
        return 0
    peak = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return peak if sys.platform == "darwin" else peak * 1024


def _run_parent(args: argparse.Namespace) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="scorelens-memory-probe-") as temporary:
        work = Path(temporary)
        write_sparse_inputs(work, rows=args.rows, genes=args.genes, density=args.density)
        command = [
            sys.executable, str(Path(__file__).resolve()), "--worker", "--work", str(work),
            "--rows", str(args.rows), "--genes", str(args.genes), "--density", repr(args.density),
            "--selected-cells", str(args.selected_cells), "--max-dense-mib", str(args.max_dense_mib),
        ]
        started = time.perf_counter()
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        child = psutil.Process(process.pid) if psutil is not None else None
        sampled_peak = 0
        while process.poll() is None:
            if child is not None:
                try:
                    sampled_peak = max(sampled_peak, int(child.memory_info().rss))
                except psutil.Error:
                    pass
            time.sleep(0.02)
        stdout, stderr = process.communicate()
        elapsed = time.perf_counter() - started
        if process.returncode != 0:
            raise RuntimeError(f"memory probe worker failed (exit {process.returncode}): {stderr.strip() or stdout.strip()}")
        worker = json.loads(stdout)
    peak = max(sampled_peak, int(worker["worker_peak_rss_bytes"]))
    if peak <= 0:
        raise RuntimeError("memory probe needs psutil or a platform with resource.getrusage")
    result = {
        "format": "vec-scorelens-memory-probe-v1",
        "synthetic_only": True,
        "logical_shape": worker["logical_shape"],
        "density": worker["density"],
        "declared_selected_shapes": worker["declared_selected_shapes"],
        "selected_target_classes": worker["selected_target_classes"],
        "source_before": worker["source_before"],
        "source_after": worker["source_after"],
        "source_unchanged": worker["source_unchanged"],
        "large_sparse_toarray_guard": worker["large_sparse_toarray_guard"],
        "report_files": worker["report_files"],
        "wall_seconds": elapsed,
        "worker_elapsed_seconds": worker["worker_elapsed_seconds"],
        "peak_rss_bytes": peak,
        "memory_limit_bytes": args.memory_limit_bytes,
        "runtime_limit_seconds": args.runtime_limit_seconds,
    }
    if peak > args.memory_limit_bytes:
        raise AssertionError(f"peak RSS {peak} exceeds full-workflow test budget {args.memory_limit_bytes} bytes")
    if elapsed > args.runtime_limit_seconds:
        raise AssertionError(f"wall time {elapsed:.2f}s exceeds limit {args.runtime_limit_seconds:.2f}s")
    return result


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=int, default=LOGICAL_ROWS)
    parser.add_argument("--genes", type=int, default=LOGICAL_GENES)
    parser.add_argument("--density", type=float, default=DENSITY)
    parser.add_argument("--selected-cells", type=int, default=DEFAULT_SELECTED_CELLS)
    parser.add_argument("--max-dense-mib", type=int, default=256)
    parser.add_argument("--memory-limit-bytes", type=int, default=DEFAULT_MEMORY_LIMIT_BYTES)
    parser.add_argument("--runtime-limit-seconds", type=float, default=DEFAULT_RUNTIME_LIMIT_SECONDS)
    parser.add_argument("--result", type=Path, help="optional path for a portable JSON evidence record")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--work", type=Path, help=argparse.SUPPRESS)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.rows < args.selected_cells or args.selected_cells < 32 or args.genes < 2:
        raise ValueError("rows must be >= selected-cells >= 32 and genes must be >= 2")
    if not 0 < args.density <= 1:
        raise ValueError("density must be in (0, 1]")
    if args.worker:
        if args.work is None:
            raise ValueError("--worker requires --work")
        print(json.dumps(_worker(args), allow_nan=False, sort_keys=True))
        return 0
    result = _run_parent(args)
    rendered = json.dumps(result, allow_nan=False, indent=2, sort_keys=True) + "\n"
    if args.result is not None:
        args.result.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
