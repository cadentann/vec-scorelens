"""Participant entrypoints, notebook execution, and protected source files."""
from pathlib import Path
import base64
import hashlib
import json
import os
import subprocess
import sys

import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
import pytest


ROOT = Path(__file__).resolve().parents[1]


def signature(path):
    stat = path.stat()
    return hashlib.sha256(path.read_bytes()).hexdigest(), stat.st_size, stat.st_mtime_ns


def test_cli_demo_report_rerun_and_no_source_mutation(tmp_path):
    env = dict(os.environ)
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    demo = subprocess.run([sys.executable, "-c", "from vec_scorelens.cli import demo_main; raise SystemExit(demo_main())",
                           "--out", str(tmp_path / "inputs"), "--seed", "0"],
                          capture_output=True, text=True, env=env, timeout=120)
    assert demo.returncode == 0, demo.stderr
    files = [tmp_path / "inputs" / f"{name}.h5ad" for name in ("prediction", "target", "reference")]
    before = [signature(p) for p in files]
    command = [sys.executable, "-m", "vec_scorelens.cli", "--task", "T1",
               "--prediction", str(files[0]), "--target", str(files[1]), "--reference", str(files[2]),
               "--out", str(tmp_path / "report")]
    result = subprocess.run(command, capture_output=True, text=True, env=env, timeout=180)
    assert result.returncode == 0, result.stderr
    assert [signature(p) for p in files] == before
    for name in ("summary.md", "diagnostics.json", "provenance.json", "gene_errors.csv", "composition_errors.csv", "covariance_errors.csv"):
        assert (tmp_path / "report" / name).is_file()
    diagnostics = json.loads((tmp_path / "report" / "diagnostics.json").read_text(),
                             parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
    terminal_warnings = [
        json.loads(line.removeprefix("vec-scorelens warning: "))
        for line in result.stderr.splitlines()
        if line.startswith("vec-scorelens warning: ")
    ]
    assert terminal_warnings == diagnostics["warnings"]
    rerun = subprocess.run(command, capture_output=True, text=True, env=env, timeout=30)
    assert rerun.returncode == 2
    assert "new output directory" in rerun.stderr.lower()
    assert "Traceback" not in rerun.stderr
    assert [signature(p) for p in files] == before


def test_cli_terminal_warning_escapes_control_characters(tmp_path, monkeypatch, capsys):
    from vec_scorelens import cli

    report = tmp_path / "report"
    report.mkdir()
    warnings = ["line one\nline two", "color \x1b[31m escape"]
    (report / "diagnostics.json").write_text(json.dumps({"warnings": warnings}))
    monkeypatch.setattr(cli, "run_task", lambda **kwargs: report)
    result = cli.main([
        "--prediction", str(tmp_path / "prediction.h5ad"),
        "--target", str(tmp_path / "target.h5ad"),
        "--reference", str(tmp_path / "reference.h5ad"),
        "--out", str(report),
    ])
    terminal = capsys.readouterr()
    lines = [line for line in terminal.err.splitlines()
             if line.startswith("vec-scorelens warning: ")]
    assert result == 0 and len(lines) == len(warnings)
    assert [json.loads(line.removeprefix("vec-scorelens warning: ")) for line in lines] == warnings
    assert "\x1b" not in terminal.err


@pytest.mark.parametrize("name", ["scorelens_synthetic.ipynb", "scorelens_t2.ipynb", "scorelens_t3.ipynb"])
def test_notebook_clean_schema_and_execution(tmp_path, name):
    notebook = nbformat.read(ROOT / "notebooks" / name, as_version=4)
    nbformat.validate(notebook)
    for cell in notebook.cells:
        if cell.cell_type == "code":
            assert cell.execution_count is None and cell.outputs == []
    manager = KernelManager(kernel_name="python3")
    manager.kernel_spec.argv = [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"]
    NotebookClient(notebook, timeout=300, km=manager,
                   resources={"metadata": {"path": str(tmp_path)}}).execute()
    for cell in notebook.cells:
        assert not any(output.output_type == "error" for output in cell.get("outputs", []))
    view = next(cell for cell in notebook.cells if cell.get("id") in ("view", "t2-summary", "t3-summary"))
    images = [output.data["image/png"] for output in view.outputs
              if output.output_type == "display_data" and "image/png" in output.data]
    assert len(images) == {"scorelens_synthetic.ipynb": 3, "scorelens_t2.ipynb": 7,
                           "scorelens_t3.ipynb": 8}[name]
    assert all(base64.b64decode(image).startswith(b"\x89PNG\r\n\x1a\n") for image in images)
    assert all("](figures/" not in output.data.get("text/markdown", "")
               for output in view.outputs if output.output_type == "display_data")
