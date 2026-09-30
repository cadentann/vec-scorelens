"""Create data-free public HTML/Markdown/figure examples from synthetic inputs.

Run from the installed repository: python examples/generate_examples.py NEW_DIR
This copies no H5ADs, row identities, source provenance or JSON/CSV ledgers.
"""
import argparse
from pathlib import Path
import shutil
import tempfile

from vec_scorelens.cli import run_task
from vec_scorelens.synthetic import write_case
from vec_scorelens.synthetic_multitask import write_task_case


EXAMPLE_NOTE = (
    "Synthetic public layout example. JSON and CSV ledgers are omitted from this "
    "display copy; a full locally generated report includes them."
)
FULL_REPORT_FOOTER = (
    "Local static report · no remote assets · full ledgers in diagnostics.json "
    "and CSV files · provenance in provenance.json"
)
EXAMPLE_FOOTER = (
    "Synthetic public layout example · no remote assets · "
    "full locally generated reports include JSON and CSV ledgers"
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("new_directory", type=Path)
    out = parser.parse_args().new_directory.resolve()
    if out.exists():
        raise SystemExit("Choose a new example output directory.")
    with tempfile.TemporaryDirectory(prefix="scorelens-examples-") as temporary:
        root = Path(temporary)
        for task in ("T1", "T2", "T3"):
            inputs = root / f"{task}-inputs"
            if task == "T1":
                write_case(inputs, case="direction")
            else:
                write_task_case(inputs, task=task, case="mixed")
            report = run_task(task=task, prediction=inputs / "prediction.h5ad",
                              target=inputs / "target.h5ad", out=root / f"{task}-report",
                              **{"wt" if task == "T3" else "reference": inputs / "reference.h5ad"})
            destination = out / task
            destination.mkdir(parents=True)
            for name in ("report.html", "summary.md"):
                text = (report / name).read_text()
                if temporary in text:
                    raise RuntimeError("Example includes temporary source identities")
                if name == "report.html":
                    if text.count(FULL_REPORT_FOOTER) != 1 or text.count("<body>") != 1:
                        raise RuntimeError("Full report HTML layout changed; review example note placement")
                    text = text.replace(
                        "<body>",
                        "<body><aside class=notice role=note>" + EXAMPLE_NOTE + "</aside>",
                        1,
                    ).replace(FULL_REPORT_FOOTER, EXAMPLE_FOOTER, 1)
                else:
                    if "\n\n" not in text:
                        raise RuntimeError("Full report Markdown layout changed; review example note placement")
                    text = text.replace("\n\n", "\n\n> " + EXAMPLE_NOTE + "\n\n", 1)
                (destination / name).write_text(text)
            shutil.copytree(report / "figures", destination / "figures")
            print(f"{task}: synthetic display example generated")


if __name__ == "__main__":
    main()
