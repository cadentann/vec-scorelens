"""Measured, local static figures for the selected-cell diagnostic report."""
from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np


def _rows(value: Any) -> list[Mapping]:
    return [row for row in value if isinstance(row, Mapping)] if isinstance(value, list) else []


def _finite(row: Mapping, name: str) -> float | None:
    try:
        value = float(row.get(name))
    except (TypeError, ValueError, OverflowError):
        return None
    return value if np.isfinite(value) else None


def _unavailable(ax, text: str) -> None:
    ax.text(0.5, 0.5, text, ha="center", va="center", wrap=True, transform=ax.transAxes)
    ax.set_axis_off()


def write_figures(bundle: Any, diagnostics: Mapping, out: Path, *, task: str, seed: int = 0) -> list[dict]:
    """Write figures from actual ledgers/arrays; unavailable evidence is explicit.

    The optional visual D2 draw is a separate seeded 10,000-pair diagnostic
    illustration, never a replacement for the official 200,000-pair score.
    No full matrix copies, network assets, or plotting state are retained.
    """
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    out.mkdir(mode=0o700)
    manifest: list[dict] = []
    colors = ("#335b9a", "#bd5a36")

    def figure(title: str, *, panels: int = 1):
        fig = Figure(figsize=(8.5, 4.8 if panels == 1 else 4.3), layout="constrained")
        FigureCanvasAgg(fig)
        axes = np.atleast_1d(fig.subplots(1, panels))
        fig.suptitle(f"{task}: {title}\nSelected supplied cells", fontsize=12)
        return fig, axes

    def save(fig, name: str, caption: str, status: str, evidence_class: str, **extra):
        path = out / name
        # Gene/class labels are untrusted literal text, not mathematical markup.
        for label in fig.findobj(match=matplotlib.text.Text):
            label.set_parse_math(False)
        with matplotlib.rc_context({"text.parse_math": False}):
            fig.savefig(path, dpi=150, facecolor="white")
        manifest.append({"path": f"figures/{name}", "caption": caption, "status": status,
                         "evidence_class": evidence_class, "scope": "selected supplied cells", **extra})
        fig.clear()

    genes = _rows(diagnostics.get("genes"))
    meaningful = [row for row in genes if _finite(row, "change_pred") is not None and _finite(row, "change_truth") is not None]
    fig, (ax,) = figure("Reference-relative gene changes" if task != "T3" else "WT-relative gene changes")
    if meaningful:
        x = np.array([float(row["change_truth"]) for row in meaningful])
        y = np.array([float(row["change_pred"]) for row in meaningful])
        ax.scatter(x, y, s=10, alpha=.35, color=colors[0])
        lo, hi = float(min(x.min(), y.min())), float(max(x.max(), y.max()))
        padding = max((hi - lo) * .05, .05)
        ax.plot([lo-padding, hi+padding], [lo-padding, hi+padding], "--", color="gray", linewidth=1, label="Equal changes")
        top = sorted(meaningful, key=lambda row: (-abs(float(row["change_pred"]) - float(row["change_truth"])), str(row.get("gene"))))[:5]
        for row, color in zip(top, ("#bd5a36", "#448454", "#7d559c", "#a17b2b", "#446b70")):
            ax.scatter([float(row["change_truth"])], [float(row["change_pred"])], s=22, color=color,
                       label=str(row.get("gene")), zorder=3)
        ax.set_xlabel("Target mean − contrast mean (transformed units)")
        ax.set_ylabel("Prediction mean − contrast mean (transformed units)")
        ax.legend(fontsize=8)
    else:
        _unavailable(ax, "Gene changes are unavailable; no values were invented.")
    save(fig, "gene_changes.png", "Per-gene measured mean changes. Contrast is the supplied reference (T1/T2) or matched WT (T3). Named genes have the largest literal change residuals; these are not biological importance rankings.", "defined" if meaningful else "unavailable", "A")

    comp = [row for row in _rows(diagnostics.get("composition")) if _finite(row, "target_fraction") is not None and _finite(row, "pred_soft_fraction") is not None]
    fig, (ax,) = figure("Target labels and frozen-probe composition")
    if comp:
        comp.sort(key=lambda row: (-abs(float(row["pred_soft_fraction"])-float(row["target_fraction"])), str(row.get("celltype"))))
        shown = comp[:12]
        pos = np.arange(len(shown))
        ax.barh(pos-.18, [float(row["target_fraction"]) for row in shown], height=.36, color=colors[0], label="Target label fraction")
        ax.barh(pos+.18, [float(row["pred_soft_fraction"]) for row in shown], height=.36, color=colors[1], label="Prediction soft fraction")
        ax.set_yticks(pos, [str(row.get("celltype")) for row in shown], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("Selected-sample proportion")
        ax.legend(fontsize=8)
    else:
        _unavailable(ax, "Target/probe composition evidence is unavailable.")
    save(fig, "composition.png", "Up to 12 classes with largest measured soft-fraction gaps. Prediction labels are ignored. A low predicted fraction is not biological absence; omitted classes remain explicit report warnings.", "defined" if comp else "unavailable", "A")

    pair_rows = [row for row in _rows(diagnostics.get("pairs")) if _finite(row, "variogram_squared_error") is not None]
    fig, (ax,) = figure("Sampled pair discrepancies")
    if pair_rows:
        shown = sorted(pair_rows, key=lambda row: (-float(row["variogram_squared_error"]), str(row.get("gene_i")), str(row.get("gene_j"))))[:10]
        pos = np.arange(len(shown))
        ax.barh(pos, [float(row["variogram_squared_error"]) for row in shown], color=colors[0])
        ax.set_yticks(pos, [f"{row.get('gene_i')} / {row.get('gene_j')}" for row in shown], fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("Squared fractional-moment residual (per unique ordered pair)")
    else:
        _unavailable(ax, "Sampled variogram evidence is unavailable.")
    save(fig, "pair_discrepancies.png", "Largest measured unique ordered-pair residuals. Official reconciliation uses retained occurrence weights, not the displayed unweighted bar heights. Covariance remains separate descriptive evidence, not regulatory interaction.", "defined" if pair_rows else "unavailable", "B", seed=seed)

    if task == "T3":
        response = diagnostics.get("response", {})
        response = response if isinstance(response, Mapping) else {}
        response_rows = [row for row in _rows(response.get("genes")) if row.get("truth_de_eligible") and _finite(row, "truth_delta") is not None and _finite(row, "pred_delta") is not None]
        summary = response.get("summary", {})
        summary = summary if isinstance(summary, Mapping) else {}
        fig, (ax,) = figure("Eligible WT-relative response and guarded fit")
        if response_rows:
            x = np.array([float(row["truth_delta"]) for row in response_rows])
            y = np.array([float(row["pred_delta"]) for row in response_rows])
            ax.scatter(x, y, s=14, alpha=.6, color=colors[0], label="Operational target-DE genes")
            lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
            ax.plot([lo, hi], [lo, hi], "--", color="gray", label="Equal response")
            beta = _finite(summary, "beta")
            if beta is not None:
                ax.plot([lo, hi], [beta*lo, beta*hi], color=colors[1], label=f"Descriptive through-origin β={beta:.3g}")
            ax.set_xlabel("Target − WT mean (transformed units)")
            ax.set_ylabel("Prediction − WT mean (transformed units)")
            ax.legend(fontsize=7)
            status = summary.get("status", "undefined")
            gate = "Source low across-gene variation guard\n(can include large uniform offsets)" if status == "no_predicted_response" else f"Gate: {status}"
            severity_label = "Official severity sentinel" if status not in ("accepted", "undefined") and _finite(summary, "official_logbeta") is not None else "Official log severity"
            # Reserve layout space above the axes; an opaque in-axes label can
            # conceal all points in a uniform-offset response.
            ax.set_title(f"{gate}\n{severity_label}: {summary.get('official_logbeta', 'undefined')}",
                         loc="left", fontsize=8, pad=10)
        else:
            _unavailable(ax, "Eligible WT-relative response evidence is unavailable or undefined.")
        save(fig, "response.png", "Measured WT-relative deltas on operational target-DE genes. A line is descriptive even when the official reliability gate rejects it; the displayed gate determines official severity interpretation. No gene residual is an allocation of log severity.", "defined" if response_rows else "unavailable", "A")

    if task in ("T2", "T3"):
        clouds = []
        for name in ("prediction", "target"):
            value = getattr(bundle, f"{name}_coords", None)
            try:
                coords = np.asarray(value, dtype=float)
                if coords.ndim != 2 or coords.shape[1] != 3 or len(coords) == 0 or not np.isfinite(coords).all():
                    coords = None
            except (TypeError, ValueError):
                coords = None
            clouds.append(coords)
        fig, axes = figure("Selected coordinate clouds in their native frames", panels=2)
        for ax, coords, label, color in zip(axes, clouds, ("Prediction", "Target"), colors):
            if coords is None:
                _unavailable(ax, f"{label} selected coordinates unavailable.")
            else:
                centered = coords - coords.mean(0)
                ax.scatter(centered[:, 0], centered[:, 1], s=7, alpha=.5, color=color)
                ax.set_xlabel("Native x − mean x (supplied units)")
                ax.set_ylabel("Native y − mean y (supplied units)")
                ax.set_aspect("equal", adjustable="box")
                ax.set_title(f"{label}: {len(coords)} selected points", fontsize=9)
        save(fig, "spatial_clouds.png", "Independent centered native XY projections. Axes are arbitrary supplied frames, unregistered and unpaired; visual differences are not pointwise errors or laterality evidence.", "defined" if all(value is not None for value in clouds) else "unavailable", "D")

        fig, (ax,) = figure("Illustrative normalized within-cloud distance distributions")
        valid = all(value is not None for value in clouds)
        distance_sets = []
        if valid:
            rng = np.random.default_rng(seed)
            for coords, label, color in zip(clouds, ("Prediction", "Target"), colors):
                centered = coords - coords.mean(0)
                rms = float(np.sqrt(np.mean(np.sum(centered**2, axis=1))))
                if rms <= 0:
                    valid = False
                    continue
                i, j = rng.integers(0, len(coords), size=(2, 10000))
                keep = i != j
                dist = np.linalg.norm(centered[i[keep]] - centered[j[keep]], axis=1) / rms
                distance_sets.append((dist, label, color))
            extent = max(6.0, max((float(dist.max()) for dist, _, _ in distance_sets if len(dist)), default=6.0))
            for dist, label, color in distance_sets:
                ax.hist(dist, bins=np.linspace(0, extent, 61), density=True, histtype="step", color=color, linewidth=1.5, label=label)
            ax.set_xlabel("Within-cloud point distance / RMS radius")
            ax.set_ylabel("Sampled density")
            ax.legend(fontsize=8)
            if not valid:
                ax.text(.5, .9, "One cloud has zero RMS radius; its distribution is undefined.", transform=ax.transAxes, ha="center", fontsize=8)
        else:
            _unavailable(ax, "Selected-coordinate distance evidence is unavailable.")
        save(fig, "shape_distances.png", "Separate illustrative draw of 10,000 ordered pairs per selected cloud; self-pairs excluded. Normalization removes size and arbitrary rigid frames. This is sampled evidence, distinct from the official D2 draw, and blind to reflection. Shared histogram bins include the complete visual draw.", "defined" if valid else "unavailable", "B", seed=seed, ordered_pair_draws_per_cloud=10000)

        spatial = diagnostics.get("spatial", {})
        spatial = spatial if isinstance(spatial, Mapping) else {}
        growth = spatial.get("growth", {})
        growth = growth if isinstance(growth, Mapping) else {}
        scale, count = growth.get("scale", {}), growth.get("count", {})
        scale = scale if isinstance(scale, Mapping) else {}
        count = count if isinstance(count, Mapping) else {}
        fig, axes = figure("Measured selected size and original/selected counts", panels=2)
        radii = [_finite(scale, name) for name in ("prediction_rms", "target_rms")]
        if all(value is not None for value in radii):
            axes[0].bar(["Prediction", "Target"], radii, color=colors)
            axes[0].set_ylabel("Selected RMS radius (supplied units)")
        else:
            _unavailable(axes[0], "Selected RMS-radius summary unavailable.")
        original = [_finite(count, name) for name in ("original_prediction_count", "original_target_count")]
        selected = [_finite(count, name) for name in ("prediction_count", "target_count")]
        if all(value is not None for value in original + selected):
            axes[1].bar(np.arange(2)-.18, original, width=.36, color=colors[0], label="Original metadata")
            axes[1].bar(np.arange(2)+.18, selected, width=.36, color=colors[1], label="Selected")
            axes[1].set_xticks([0, 1], ["Prediction", "Target"])
            axes[1].set_ylabel("Cell count")
            axes[1].legend(fontsize=7)
        else:
            _unavailable(axes[1], "Original/selected counts unavailable.")
        save(fig, "spatial_size_counts.png", "Direct size and count summaries. The official count scalar uses selected counts; original counts are descriptive metadata only. Neither identifies biological growth or proliferation.", "defined" if all(value is not None for value in radii + original + selected) else "unavailable", "D")
        local = spatial.get("local", {})
        local = local if isinstance(local, Mapping) else {}
        neighborhood = local.get("neighborhood_mmd", {})
        neighborhood = neighborhood if isinstance(neighborhood, Mapping) else {}
        anchors = [row for row in _rows(neighborhood.get("anchor_rows")) if _finite(row, "allocation") is not None
                   and isinstance(row.get("coordinates"), (list, tuple)) and len(row["coordinates"]) == 3]
        fig, axes = figure("Local neighborhood signed kernel terms", panels=2)
        bound = max(max((abs(float(row["allocation"])) for row in anchors), default=0.0), 1e-12)
        for ax, side in zip(axes, ("prediction", "truth")):
            side_rows = [row for row in anchors if row.get("side") == side]
            if side_rows:
                coords = np.array([row["coordinates"] for row in side_rows], dtype=float)
                terms = np.array([float(row["allocation"]) for row in side_rows])
                centered = coords - coords.mean(0)
                points = ax.scatter(centered[:, 0], centered[:, 1], c=terms, s=14, cmap="coolwarm", vmin=-bound, vmax=bound)
                fig.colorbar(points, ax=ax, label="Signed MMD kernel allocation", shrink=.75)
                ax.set_xlabel("Native x − mean x")
                ax.set_ylabel("Native y − mean y")
                ax.set_title(f"{side}: {len(side_rows)} selected anchors", fontsize=9)
                ax.set_aspect("equal", adjustable="box")
            else:
                _unavailable(ax, f"{side} signed anchor terms unavailable.")
        save(fig, "local_neighborhood.png", "Signed finite-sample kernel terms in independent native XY frames. Both prediction and truth terms can be negative; they are not per-cell causal errors or anatomical correspondence. No local-priority threshold is adopted; effective upstream MMD seed is 0.", "defined" if anchors else "unavailable", "A", effective_seed=0)
    return manifest
