"""Short sentinel for the portable sparse bounded-memory acceptance probe."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys


def test_sparse_logical_matrix_is_scored_without_large_densification(tmp_path: Path) -> None:
    """Exercise the complete workflow on a smaller selected subset for CI speed."""
    probe = Path(__file__).resolve().parents[1] / "benchmarks" / "memory_probe.py"
    result_path = tmp_path / "memory-result.json"
    completed = subprocess.run(
        [
            sys.executable, str(probe), "--selected-cells", "64", "--result", str(result_path),
            "--runtime-limit-seconds", "180",
        ],
        capture_output=True,
        text=True,
        timeout=210,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr or completed.stdout
    result = json.loads(result_path.read_text(encoding="utf-8"))
    assert result["synthetic_only"] is True
    assert result["logical_shape"] == [10_000, 32_285]
    assert result["large_sparse_toarray_guard"] == "passed"
    assert result["source_unchanged"] is True
    assert set(result["selected_target_classes"]) == {"type_a", "type_b"}
    assert all(shape == [64, 32_285] for shape in result["declared_selected_shapes"].values())
    assert result["peak_rss_bytes"] <= result["memory_limit_bytes"]
    assert result["wall_seconds"] <= result["runtime_limit_seconds"]
    for name in ("prediction", "target", "reference"):
        assert result["source_before"][name] == result["source_after"][name]
