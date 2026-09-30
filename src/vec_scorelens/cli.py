"""Small command-line front door for supplied local files and synthetic examples."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import warnings

from threadpoolctl import threadpool_limits


def run(prediction, target, reference, out, *, seed=0, max_cells=256, max_dense_mib=256):
    """Preserved T1 entrypoint for supplied local files."""
    return run_task(task="T1", prediction=prediction, target=target, reference=reference,
                    out=out, seed=seed, max_cells=max_cells, max_dense_mib=max_dense_mib)


def _validate_roles(task, reference, wt):
    if task not in ("T1", "T2", "T3"):
        raise ValueError(f"unsupported task {task!r}; choose T1, T2, or T3")
    if task == "T3":
        if reference is not None:
            raise ValueError("T3 accepts --wt only; --reference is for T1/T2")
        if wt is None:
            raise ValueError("T3 requires an explicit --wt file")
    else:
        if wt is not None:
            raise ValueError(f"{task} accepts --reference only; --wt is for T3")
        if reference is None:
            raise ValueError(f"{task} requires an explicit --reference file")


def run_task(*, task, prediction, target, out, reference=None, wt=None,
             seed=0, max_cells=256, max_dense_mib=256):
    """Inspect one bounded task and atomically write a new report."""
    _validate_roles(task, reference, wt)
    from .io import load_task_bundle
    from .scorer import run_official
    from .diagnostics import diagnose
    from .report import write_report

    destination = Path(out)
    if destination.exists() or destination.is_symlink():
        raise ValueError(f"Output already exists: {destination}. Choose a new output directory.")
    with threadpool_limits(limits=1), warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        bundle = load_task_bundle(task=task, prediction=prediction, target=target,
                                  reference=reference, wt=wt, max_cells=max_cells,
                                  seed=seed, max_dense_mib=max_dense_mib)
        scored = run_official(bundle, seed=seed)
        diagnostics = diagnose(bundle, scored, seed=seed)
    diagnostics.setdefault("warnings", []).extend(
        sorted({f"Runtime {w.category.__name__}: {w.message}" for w in caught}))
    return write_report(bundle, scored, diagnostics, destination)


def main(argv=None):
    parser = argparse.ArgumentParser(description=(
        "Local T1/T2/T3 prediction diagnostics against supplied permissible pseudo-targets. "
        "Default: deterministic sample of at most 256 cells per file; all genes retained."))
    parser.add_argument("--task", choices=["T1", "T2", "T3"], default="T1")
    parser.add_argument("--prediction", required=True, type=Path)
    parser.add_argument("--target", required=True, type=Path)
    parser.add_argument("--reference", type=Path, help="Required preceding-stage reference for T1/T2")
    parser.add_argument("--wt", type=Path, help="Required matched wild-type reference for T3")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-cells", type=int, default=256)
    parser.add_argument("--max-dense-mib", type=float, default=256)
    args = parser.parse_args(argv)
    try:
        # Validate roles before opening any input or checking an output path.
        _validate_roles(args.task, args.reference, args.wt)
        report = run_task(task=args.task, prediction=args.prediction, target=args.target,
                          reference=args.reference, wt=args.wt, out=args.out,
                          seed=args.seed, max_cells=args.max_cells, max_dense_mib=args.max_dense_mib)
        report_warnings = json.loads(
            (report / "diagnostics.json").read_text(encoding="utf-8")
        )["warnings"]
        if not isinstance(report_warnings, list) or not all(
            isinstance(warning, str) for warning in report_warnings
        ):
            raise ValueError("generated report warning list is invalid")
    except (ValueError, OSError, TypeError, RuntimeError, KeyError) as exc:
        print(f"vec-scorelens: {exc}", file=sys.stderr)
        return 2
    print(f"Wrote local {args.task} diagnostics: {report}")
    print("Read summary.md for sample scope, diagnostic definitions, findings, and limitations.")
    for warning in report_warnings:
        print(f"vec-scorelens warning: {json.dumps(warning, ensure_ascii=True)}", file=sys.stderr)
    return 0


def demo_main(argv=None):
    parser = argparse.ArgumentParser(description="Create synthetic local ScoreLens inputs in a new directory.")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--task", choices=["T1", "T2", "T3"], default="T1")
    parser.add_argument("--case")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)
    try:
        if args.task == "T1":
            from .synthetic import write_case
            write_case(args.out, case=args.case or "direction", seed=args.seed)
        else:
            from .synthetic_multitask import write_task_case
            write_task_case(args.out, task=args.task, case=args.case, seed=args.seed)
    except (ValueError, OSError, TypeError, RuntimeError, KeyError) as exc:
        print(f"vec-scorelens-demo: {exc}", file=sys.stderr)
        return 2
    print(f"Wrote synthetic inputs: {args.out}")
    print("These fixtures are local experiments, not challenge submissions.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
