"""Render a report illustration; optional install: pip install '.[figures]'.

Usage: python examples/render_synthetic_figure.py REPORT_DIRECTORY NEW_OUTPUT.png
Use only synthetic reports for public examples. This script is a presentation
helper; it does not compute metrics or interpret gene function.
"""
from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("report", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output exists; choose a new filename")
    genes = pd.read_csv(args.report / "gene_errors.csv")
    top = genes.sort_values(["absolute_error", "gene"], ascending=[False, True]).head(8)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    axes[0].scatter(genes.change_truth, genes.change_pred, s=18, color="#8799ae", label="All synthetic genes")
    axes[0].scatter(top.change_truth, top.change_pred, s=36, color="#c24735", label="Eight largest mean errors")
    bound = max(np.abs(genes.change_truth).max(), np.abs(genes.change_pred).max()) * 1.12
    axes[0].plot([-bound, bound], [-bound, bound], "--", color="#64748b", lw=1)
    axes[0].axhline(0, color="#cbd5e1", lw=.7)
    axes[0].axvline(0, color="#cbd5e1", lw=.7)
    axes[0].set(xlabel="Target mean shift from reference", ylabel="Prediction mean shift from reference",
                title="Selected synthetic cells: change direction")
    axes[0].legend(fontsize=8, loc="upper left")
    y = np.arange(len(top))
    axes[1].barh(y - .18, top.change_truth, height=.35, label="Target", color="#327da8")
    axes[1].barh(y + .18, top.change_pred, height=.35, label="Prediction", color="#c24735")
    axes[1].set_yticks(y, top.gene)
    axes[1].invert_yaxis()
    axes[1].axvline(0, color="#64748b", lw=.8)
    axes[1].set(xlabel="Mean log-expression shift (synthetic units)", title="Named discrepancies in the report")
    axes[1].legend(fontsize=8)
    fig.suptitle("ScoreLens synthetic direction example — descriptive evidence, not causal attribution", fontsize=12)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
