"""Write a private, portable ScoreLens diagnostic report.

The report is deliberately a *presentation* boundary: diagnostic calculations live
elsewhere and this module never attempts to recompute or rank official metrics.
"""

from __future__ import annotations

import csv
import html
import json
import math
import os
import re
import shutil
import uuid
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any


_FORMULA_PREFIXES = ("=", "+", "-", "@")
_GENE_COLUMNS = (
    "gene", "mean_pred", "mean_truth", "mean_reference", "change_pred",
    "change_truth", "change_error", "absolute_error", "squared_error",
    "pb_squared_allocation", "variance_pred", "variance_truth", "variance_ratio",
    "truth_de_direction", "pred_local_de_direction", "local_de_hit",
    "local_de_miss", "local_de_false_positive", "official_rank", "official_slot",
    "official_slot_hit", "literal_pred_direction", "literal_sign_agreement",
    "slot_hit_wrong_literal_sign", "de_raw_allocation", "de_score_allocation",
)
_COMPOSITION_COLUMNS = (
    "celltype", "target_count", "target_fraction", "pred_soft_fraction",
    "proportion_error", "smoothed_target", "smoothed_pred", "jsd_term",
)
_PAIR_COLUMNS = (
    "gene_i", "gene_j", "sampled_occurrences", "variogram_pred",
    "variogram_truth", "variogram_error", "variogram_squared_error",
    "variogram_contribution", "covariance_pred", "covariance_truth",
    "covariance_error",
)


def _json_safe(value: Any, *, path: str, notes: dict[str, int]) -> Any:
    """Convert values to strict JSON and record undefined numeric values."""
    # Numpy scalar support without a hard numpy dependency at import time.
    item = getattr(value, "item", None)
    if callable(item) and not isinstance(value, (str, bytes, bytearray)):
        try:
            return _json_safe(item(), path=path, notes=notes)
        except (TypeError, ValueError):
            pass
    if isinstance(value, float):
        if not math.isfinite(value):
            # A gene-level undefined value can legitimately recur tens of
            # thousands of times. Keep a field-level count instead of bloating
            # diagnostics.json and summary.md with one warning per row.
            location = re.sub(r"\[\d+\]", "[]", path)
            notes[location] = notes.get(location, 0) + 1
            return None
        return value
    if isinstance(value, Mapping):
        return {
            str(key): _json_safe(child, path=f"{path}.{key}", notes=notes)
            for key, child in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [_json_safe(child, path=f"{path}[{index}]", notes=notes) for index, child in enumerate(value)]
    if value is None:
        return None
    if isinstance(value, (str, int, bool)):
        return value
    # Diagnostics are contractually JSON serialisable, but preserve a useful
    # representation rather than letting a report fail on an incidental scalar.
    return str(value)


def _csv_cell(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        # Spreadsheet applications may evaluate a formula after leading spaces or
        # tabs, so examine the stripped value but retain the visible label.
        if value.lstrip(" \t\r\n").startswith(_FORMULA_PREFIXES):
            return "'" + value
    return value


def _write_csv(path: Path, rows: Iterable[Any], preferred_columns: tuple[str, ...]) -> None:
    materialized = [row for row in rows if isinstance(row, Mapping)]
    extras = sorted({str(key) for row in materialized for key in row} - set(preferred_columns))
    columns = list(preferred_columns) + extras
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([_csv_cell(column) for column in columns])
        for row in materialized:
            writer.writerow([_csv_cell(row.get(column)) for column in columns])


def _markdown(value: Any) -> str:
    """Make an untrusted label safe in normal Markdown table text."""
    text = html.escape("" if value is None else str(value), quote=False)
    for char in ("\\", "`", "*", "_", "[", "]", "#", "+", "-", "!", "|"):
        text = text.replace(char, "\\" + char)
    return text.replace("\r", "").replace("\n", "<br>")


def _value(value: Any) -> str:
    if value is None:
        return "null (undefined/not applicable)"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def _contains_null(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, Mapping):
        return any(_contains_null(child) for child in value.values())
    if isinstance(value, list):
        return any(_contains_null(child) for child in value)
    return False


def _warning_strings(value: Any) -> list[str]:
    """Normalize warning inputs and retain first-seen order for de-duplication."""
    if not isinstance(value, list):
        return []
    return list(dict.fromkeys(str(item) for item in value))


def _top(rows: Any, field: str, count: int = 5) -> list[Mapping[str, Any]]:
    if not isinstance(rows, list):
        return []

    def key(row: Any) -> float:
        if not isinstance(row, Mapping):
            return float("-inf")
        value = row.get(field)
        return abs(float(value)) if isinstance(value, (int, float)) and math.isfinite(float(value)) else float("-inf")

    return [row for row in sorted((row for row in rows if isinstance(row, Mapping)), key=key, reverse=True)[:count]]


def _table(lines: list[str], headers: tuple[str, ...], rows: Iterable[Mapping[str, Any]]) -> None:
    rows = list(rows)
    if not rows:
        lines.append("No defined rows were available for this section.")
        return
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join("---" for _ in headers) + " |")
    for row in rows:
        lines.append("| " + " | ".join(_markdown(_value(row.get(header))) for header in headers) + " |")


def _summary(bundle: Any, scored: Mapping[str, Any], diagnostics: Mapping[str, Any]) -> str:
    metrics = scored.get("metrics", {}) if isinstance(scored.get("metrics", {}), Mapping) else {}
    genes = diagnostics.get("genes", [])
    composition = diagnostics.get("composition", [])
    pairs = diagnostics.get("pairs", [])
    gene_summary = diagnostics.get("gene_summary", {}) if isinstance(diagnostics.get("gene_summary", {}), Mapping) else {}
    composition_summary = diagnostics.get("composition_summary", {}) if isinstance(diagnostics.get("composition_summary", {}), Mapping) else {}
    pair_summary = diagnostics.get("pair_summary", {}) if isinstance(diagnostics.get("pair_summary", {}), Mapping) else {}
    distribution = diagnostics.get("distribution", {}) if isinstance(diagnostics.get("distribution", {}), Mapping) else {}
    task = str(diagnostics.get("task", getattr(bundle, "task", "T1")))
    contrast = "matched WT" if task == "T3" else "supplied reference"

    lines = [
        f"# VEC ScoreLens {task} private diagnostic report — selected subset",
        "",
        "This report may identify supplied data. Keep it private and do not treat selected-subset results as full-file certification.",
        "",
        "## Diagnostic follow-up priorities",
        "",
        *_priority_markdown(diagnostics.get("priorities", {})),
        "",
        "## Assessed input scope",
        "",
        f"The supplied ordered gene panel and independently selected cells are assessed. The contrast is the {contrast}; equal cell counts do not create paired cells.",
        "",
        *_scope_markdown(diagnostics.get("scope", {})),
        "",
        "## Official score output",
        "",
        "The values below are the supplied official score output for the selected subset. They are different scales and definitions, so this report intentionally does not label any one metric the ‘weakest’. ",
    ]
    if metrics:
        _table(lines, ("metric", "value"), [{"metric": key, "value": value} for key, value in metrics.items()])
    else:
        lines.append("No official metrics were supplied.")

    lines.extend(["", "## Concrete gene follow-up", "", "Genes below have the largest defined per-gene pseudobulk squared residuals in the selected data. Check their input alignment, normalization, and target/reference contrast before changing a model."])
    lines.append(f"A — direct selected-sample evidence. `change_truth` is target mean minus {contrast} mean; `change_pred` is prediction mean minus {contrast} mean, in supplied transformed expression units.")
    _table(lines, ("gene", "change_truth", "change_pred", "absolute_error", "squared_error", "truth_de_direction", "pred_local_de_direction", "official_slot_hit"), _top(genes, "squared_error"))
    lines.append("")
    lines.append("Local significant-set hit/miss/false-positive fields and literal mean-log-expression signs are diagnostic aids. A false-positive means a prediction-only local DE call, not a known biological false response or an FDR conclusion; target nonsignificance does not mean no change. These fields are not the official target-sized ranked slots or official partial rank correlation. Rank-based direction can change with float32 roundoff near zero shifts; tiny literal sign differences are not evidence of a meaningful response.")
    if gene_summary:
        reconciliation = gene_summary.get("pb_reconciliation_abs_error", gene_summary.get("reconciliation_abs_error"))
        lines.append(f"Pseudobulk allocation reconciliation error: {_markdown(_value(reconciliation))}. Guard state: {_markdown(_value(gene_summary.get('de_guard')))}.")
        count_keys = ("n_truth_up", "n_truth_down", "n_local_hits", "n_local_misses", "n_local_false_positives", "n_slot_hits_wrong_literal_sign")
        if all(key in gene_summary for key in count_keys):
            lines.append(
                "Target local DE calls: "
                f"{_markdown(_value(gene_summary.get('n_truth_up')))} up and {_markdown(_value(gene_summary.get('n_truth_down')))} down; "
                f"local signed matches {_markdown(_value(gene_summary.get('n_local_hits')))}, misses {_markdown(_value(gene_summary.get('n_local_misses')))}, "
                f"and prediction-only calls {_markdown(_value(gene_summary.get('n_local_false_positives')))}. "
                f"Target rank-slot hits with a wrong literal sign: {_markdown(_value(gene_summary.get('n_slot_hits_wrong_literal_sign')))}."
            )
    if isinstance(distribution, Mapping) and "variance_ratio" in distribution:
        lines.append(f"Distribution variance-ratio evidence (selected supplied cells): {_markdown(_value(distribution.get('variance_ratio')))}. A low value is descriptive collapse evidence, not a biological conclusion.")

    lines.extend(["", "## Concrete composition follow-up", "", "Classes below have the largest defined target-versus-probe soft-proportion deltas. Inspect target labels and probe-domain fit; prediction labels are intentionally not used. Bounded sampling can omit rare classes, and a low predicted proportion is not biological absence."])
    _table(lines, ("celltype", "target_fraction", "pred_soft_fraction", "proportion_error", "jsd_term"), _top(composition, "proportion_error"))
    if composition_summary:
        lines.append(f"JSD term reconciliation error: {_markdown(_value(composition_summary.get('reconciliation_abs_error')))}.")

    lines.extend(["", "## Sampled variogram and descriptive covariance follow-up", "", "These rows aggregate sampled ordered pairs (duplicates retained as occurrences). Variogram errors are sampled fractional-moment errors; covariance residuals are descriptive companion columns, not a claim that variogram is pure covariance. Covariance uses population ddof=0 on the supplied transformed values, so its units follow those values rather than a biological scale."])
    _table(lines, ("gene_i", "gene_j", "sampled_occurrences", "variogram_squared_error", "covariance_error"), _top(pairs, "variogram_squared_error"))
    if pair_summary:
        lines.append(f"Sampled pair count: {_markdown(_value(pair_summary.get('n_sampled_pairs')))}; unique ordered pairs: {_markdown(_value(pair_summary.get('n_unique_ordered_pairs')))}; reconciliation error: {_markdown(_value(pair_summary.get('reconciliation_abs_error')))}.")
        lines.append(f"Possible ordered nonself pairs: {_markdown(_value(pair_summary.get('total_possible_ordered_pairs')))}; unique ordered-pair coverage fraction: {_markdown(_value(pair_summary.get('unique_ordered_pair_coverage')))}. {_markdown(pair_summary.get('coverage_interpretation', 'Unsampled pairs have no diagnostic evidence.'))}")

    sensitivity = distribution.get("mmd_sensitivity", distribution.get("mmd", {}))
    lines.extend(["", "## MMD sensitivity (noncausal)", ""])
    if isinstance(sensitivity, Mapping):
        affected = sensitivity.get("affected_genes", sensitivity.get("genes", []))
        lines.append(
            "One nonnegative-clipped target mean-shift sensitivity was run: "
            f"before {_markdown(_value(sensitivity.get('before')))}, after {_markdown(_value(sensitivity.get('after')))}, "
            f"change {_markdown(_value(sensitivity.get('change')))}; affected genes: {_markdown(', '.join(map(str, affected)) if isinstance(affected, list) else affected)}."
        )
        lines.append("The listed delta is after minus before. This target-aware, selected-and-evaluated-on-the-same-target artificial diagnostic is not a validated model correction, causal or additive attribution, percent contribution, or expected training improvement; it can worsen the score.")
    else:
        lines.append("No MMD sensitivity result was supplied.")

    lines.extend(_task_markdown(task, diagnostics))
    lines.extend(["", "## Measured figures", ""])
    for figure in diagnostics.get("figures", []):
        if isinstance(figure, Mapping):
            lines.append(f"![{_markdown(figure.get('caption'))}]({figure.get('path')})")
            lines.append(f"{_markdown(figure.get('evidence_class'))} — {_markdown(figure.get('status'))}: {_markdown(figure.get('caption'))}")
            lines.append("")
    lines.extend(["", "## Evidence taxonomy and limits", "", "A: direct or algebraic selected-sample evidence; B: sampled evidence; C: artificial noncausal sensitivity; D: descriptive companion statistics. Exact/direct does not mean exact biology. Gene and composition calculations are exact for the bounded selected cells and supplied probe within stated numerical reconciliation tolerance; pair calculations are sampled. The source-provided gene panel is used as supplied and this report is not an official leaderboard or board-validation result. The selected subset is explicitly not evidence of full-file equivalence. Probe labels are limited to its domain, so a low predicted proportion is not biological absence. Variogram/covariance diagnostics do not establish gene regulation. Any max-dense input limit controls selected input conversion only; it is not a total process peak-memory guarantee."])
    warnings = diagnostics.get("warnings", [])
    if isinstance(warnings, list) and warnings:
        lines.extend(["", "## Warnings", ""])
        lines.extend("- " + _markdown(warning) for warning in warnings)
    return "\n".join(lines) + "\n"


def _scope_rows(scope: Any) -> list[dict]:
    scope = scope if isinstance(scope, Mapping) else {}
    selected = scope.get("selected_counts", {})
    original = scope.get("original_counts", {})
    selected = selected if isinstance(selected, Mapping) else {}
    original = original if isinstance(original, Mapping) else {}
    return [{"input": role, "original_cells": original.get(role), "selected_cells": selected.get(role)}
            for role in dict.fromkeys([*selected, *original])]


def _scope_markdown(scope: Any) -> list[str]:
    lines: list[str] = []
    _table(lines, ("input", "original_cells", "selected_cells"), _scope_rows(scope))
    return lines


def _priority_markdown(value: Any) -> list[str]:
    priority = value if isinstance(value, Mapping) else {}
    assessment = str(priority.get("assessment", "not_evaluated")).replace("_", " ")
    lines = [f"Assessment: **{assessment}**. Coverage: {_markdown(priority.get('coverage', 'not evaluated'))}.",
             "These explicit heuristic conventions order diagnostic follow-up. They are not statistically calibrated, have no guaranteed false-alarm rate or FDR, and do not establish biological health or expected model improvement."]
    profile = priority.get("profile", {})
    if isinstance(profile, Mapping):
        lines.append(f"Profile: {_markdown(profile.get('id', 'none'))}. At most {_markdown(profile.get('max_findings', 6))} findings are displayed in fixed domain order; all evaluated thresholds are recorded in diagnostics.json.")
    for finding in priority.get("findings", []):
        if not isinstance(finding, Mapping):
            continue
        lines.extend(["", f"### {_markdown(finding.get('priority_order'))}. {_markdown(finding.get('title'))}", "",
                      f"Evidence {_markdown(finding.get('evidence_class'))}; {_markdown(finding.get('scope'))}. Triggered rows: {_markdown(finding.get('n_triggered'))}.",
                      _markdown(finding.get("suggested_check"))])
        _table(lines, ("entity", "measured_value", "comparator", "threshold", "units"), finding.get("evidence", []))
    if priority.get("n_findings_omitted"):
        lines.append(f"{_markdown(priority['n_findings_omitted'])} additional triggered rule(s) are retained in the full evaluation ledger.")
    unavailable = priority.get("not_evaluated", [])
    if unavailable:
        lines.extend(["", "Coverage exclusions:"])
        lines.extend(f"- {_markdown(row.get('rule_id'))}: {_markdown(row.get('reason'))}" for row in unavailable if isinstance(row, Mapping))
    lines.append("Within listed thresholds means only that evaluated listed rules did not trigger on supplied selected cells. Undefined evidence and unassessed domains remain explicit; this is not a full-file all-clear.")
    return lines


def _task_markdown(task: str, diagnostics: Mapping) -> list[str]:
    lines: list[str] = []
    if task == "T3":
        response = diagnostics.get("response", {})
        response = response if isinstance(response, Mapping) else {}
        summary = response.get("summary", {})
        summary = summary if isinstance(summary, Mapping) else {}
        lines.extend(["", "## WT-relative response", "", "A — direct operational response evidence. WT matching, target DE eligibility, literal signs, magnitude ratios, and the guarded through-origin severity fit are separate concepts."])
        _table(lines, ("field", "value"), [{"field": key, "value": summary.get(key)} for key in
               ("status", "n_truth_de", "n_sign_reversals", "beta", "official_logbeta", "official_r2", "magnitude_median_ratio", "n_magnitude_eligible", "reported_logbeta_origin")])
        _table(lines, ("gene", "truth_delta", "pred_delta", "literal_sign_status", "magnitude_ratio", "magnitude_status", "slope_residual"),
               _top(response.get("genes", []), "direct_error"))
        lines.append("The n_sign_reversals count covers opposite tolerance-based signs across all supplied genes, including genes without a target DE call; it is not a count of meaningful DE response reversals. A failure sentinel is not a measured log amplitude. The source no_predicted_response low across-gene variation guard compares standard deviations of all deltas, so a constant nonzero uniform offset can also trigger it; it does not mean every literal delta is zero. Source-guarded, inversion/unrelated response, and undefined truth branches must be read through the displayed gate. Gene slope residuals do not sum to log severity, and these selected-cell operational calls are not confirmed biological responses.")
    if task in ("T2", "T3"):
        spatial = diagnostics.get("spatial", {})
        spatial = spatial if isinstance(spatial, Mapping) else {}
        lines.extend(["", "## Spatial form, size, counts, and local organization", "",
                      "Coordinates are selected with expression row identities. Cells are unpaired and coordinate frames arbitrary. D2 is reflection blind; no calibrated laterality claim is possible. SW/Dice retain a source-frame caveat and receive no priority threshold. Local organization is measured without a calibrated follow-up threshold."])
        for family_name in ("shape", "growth", "local"):
            family = spatial.get(family_name, {})
            if not isinstance(family, Mapping):
                continue
            for section_name, section in family.items():
                if not isinstance(section, Mapping):
                    continue
                official_label = "not part of the official T3 panel (extra diagnostic only)" if task == "T3" and section.get("official_panel_includes_metric") is False else _value(section.get("official_value"))
                lines.extend(["", f"### {_markdown(family_name)} / {_markdown(section_name)}", "",
                              f"{_markdown(section.get('definition', 'Undefined diagnostic definition'))}",
                              f"Status: {_markdown(section.get('status', 'undefined'))}; official: {_markdown(official_label)}, replay {_markdown(_value(section.get('replayed_value')))}, reconciliation {_markdown(_value(section.get('reconciliation_abs_error')))}."])
                for limitation in section.get("limitations", []):
                    lines.append("- " + _markdown(limitation))
                if section_name == "neighborhood_mmd":
                    lines.extend(["", "Signed finite-sample two-sided anchor kernel allocations, ordered by absolute term magnitude within this table only. These are neither causal cell errors nor diagnostic priorities; coordinates remain in each side's own unregistered frame."])
                    _table(lines, ("side", "row_position", "coordinates", "allocation"), _top(section.get("anchor_rows", []), "allocation", 10))
        growth = spatial.get("growth", {})
        growth = growth if isinstance(growth, Mapping) else {}
        count = growth.get("count", {})
        count = count if isinstance(count, Mapping) else {}
        lines.append(f"Original prediction/target counts: {_markdown(_value(count.get('original_prediction_count')))} / {_markdown(_value(count.get('original_target_count')))}; selected: {_markdown(_value(count.get('prediction_count')))} / {_markdown(_value(count.get('target_count')))}. Original count log ratio: {_markdown(_value(count.get('original_log_ratio')))} (D — metadata only). The official count_log_ratio uses selected counts and can be zero after both files reach the cap. Do not infer biological growth or proliferation.")
    return lines


def _html_table(headers: tuple[str, ...], rows: Iterable[Mapping]) -> str:
    rows = [row for row in rows if isinstance(row, Mapping)]
    if not rows:
        return "<p class=muted>No defined rows were available.</p>"
    head = "".join(f"<th scope=col>{html.escape(key.replace('_', ' '))}</th>" for key in headers)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(_value(row.get(key)), quote=True)}</td>" for key in headers) + "</tr>" for row in rows)
    return f"<div class=table-wrap><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>"


def _html_report(bundle: Any, scored: Mapping, diagnostics: Mapping) -> str:
    task = str(diagnostics.get("task", "T1"))
    esc = lambda value: html.escape(str(value), quote=True)
    priority = diagnostics.get("priorities", {})
    priority = priority if isinstance(priority, Mapping) else {}
    profile = priority.get("profile", {})
    profile = profile if isinstance(profile, Mapping) else {}
    title = f"VEC ScoreLens {task} private diagnostic report"
    content = [f"<header><p class=eyebrow>{esc(task)} · selected supplied cells · private</p><h1>{esc(title)}</h1><p>This report assesses a bounded selected subset and the supplied gene panel. Keep it private; it does not certify unscanned rows or hidden-target performance.</p></header>",
               "<nav aria-label='Report sections'><a href='#priorities'>Follow-up</a><a href='#scope'>Scope</a><a href='#official'>Official output</a><a href='#genes'>Genes</a><a href='#composition'>Composition</a><a href='#pairs'>Pairs</a><a href='#task'>Task evidence</a><a href='#figures'>Figures</a><a href='#warnings'>Warnings</a></nav>",
               f"<section id=priorities><h2>Diagnostic follow-up</h2><div class=notice><strong>{esc(str(priority.get('assessment', 'not_evaluated')).replace('_', ' '))}</strong><p>Coverage: {esc(priority.get('coverage', 'not evaluated'))}. These heuristic conventions are not statistically calibrated, carry no guaranteed false-alarm rate or FDR, and do not establish biological health or model improvement.</p><p>Profile {esc(profile.get('id', 'none'))}; fixed domain ordering, at most six displayed findings. Unlike official scalar scales are never compared.</p></div>"]
    for finding in priority.get("findings", []):
        if isinstance(finding, Mapping):
            content.extend([f"<article class=finding><h3>{esc(finding.get('priority_order'))}. {esc(finding.get('title'))}</h3><p><span class=badge>{esc(finding.get('evidence_class'))}</span> {esc(finding.get('scope'))}; {esc(finding.get('n_triggered'))} triggered row(s).</p>",
                            _html_table(("entity", "measured_value", "comparator", "threshold", "units"), finding.get("evidence", [])),
                            f"<p>{esc(finding.get('suggested_check'))}</p></article>"])
    if priority.get("n_findings_omitted"):
        content.append(f"<p>{esc(priority['n_findings_omitted'])} further triggered rule(s) remain in the JSON evaluation ledger.</p>")
    content.append("<details open><summary>Unassessed or unavailable follow-up rules</summary><ul>")
    for row in priority.get("not_evaluated", []):
        if isinstance(row, Mapping):
            content.append(f"<li><strong>{esc(row.get('rule_id'))}</strong>: {esc(row.get('reason'))}</li>")
    content.extend(["</ul></details><p>Within listed thresholds means only that evaluated rules did not trigger. It is not a biological or full-file all-clear. Threshold calculations for every assessed row are in diagnostics.json.</p></section>",
                    "<section id=scope><h2>Assessed scope</h2><p>Independent deterministic selection preserves expression/coordinate row correspondence within each file; it does not establish paired cells across files.</p>",
                    _html_table(("input", "original_cells", "selected_cells"), _scope_rows(diagnostics.get("scope"))), "</section>",
                    "<section id=official><h2>Official scalar output — selected subset</h2><p>These values are supplied unchanged by the pinned task v2 scorer. Different metrics have different scales and definitions; no weakest metric or combined priority score is inferred.</p>",
                    _html_table(("metric", "value"), [{"metric": key, "value": value} for key, value in scored.get("metrics", {}).items()]), "</section>",
                    "<section id=genes><h2>Gene follow-up <span class=badge>A</span></h2><p>Changes are mean transformed expression relative to the supplied reference (T1/T2) or matched WT (T3), not conventional fold changes. Operational DE calls, ranked slots, and literal signs remain separate. Target nonsignificance does not prove no biological response.</p>",
                    _html_table(("gene", "change_truth", "change_pred", "absolute_error", "truth_de_direction", "pred_local_de_direction", "official_slot_hit"), _top(diagnostics.get("genes"), "squared_error", 10)),
                    "<p>Local false-positive means a prediction-only significant signed call, including a reversal; it is not known false biology or an FDR estimate. Full rows are in gene_errors.csv.</p></section>",
                    "<section id=composition><h2>Frozen-probe composition <span class=badge>A</span></h2><p>Prediction labels are ignored. Fractions use the selected target-fitted official probe; a low fraction does not establish biological absence.</p>",
                    _html_table(("celltype", "target_fraction", "pred_soft_fraction", "proportion_error", "jsd_term"), _top(diagnostics.get("composition"), "proportion_error", 10)), "</section>",
                    "<section id=pairs><h2>Sampled pairs and descriptive covariance <span class=badge>B / D</span></h2><p>Ordered-pair occurrences and multiplicities are retained for variogram reconciliation. Fractional moments are not pure covariance; population ddof=0 covariance companion columns do not identify regulatory edges.</p>",
                    _html_table(("gene_i", "gene_j", "sampled_occurrences", "variogram_squared_error", "covariance_error"), _top(diagnostics.get("pairs"), "variogram_squared_error", 10)), "</section>",
                    ""])
    pair_summary = diagnostics.get("pair_summary", {})
    if isinstance(pair_summary, Mapping):
        # Insert coverage beside pair evidence rather than leaving it only in JSON.
        content.insert(-2, _html_table(("field", "value"), [{"field": key, "value": pair_summary.get(key)} for key in
            ("n_sampled_pairs", "n_unique_ordered_pairs", "total_possible_ordered_pairs", "unique_ordered_pair_coverage", "observed_gene_coverage", "reconciliation_abs_error")]) + f"<p>{esc(pair_summary.get('coverage_interpretation', 'Unsampled pairs have no diagnostic evidence.'))}</p>")
    content.append("<section id=sensitivity><h2>Target-aware MMD sensitivity <span class=badge>C</span></h2><p>One clipped mean-shift intervention is selected and evaluated on the same target. It is noncausal and nonadditive, may worsen MMD, and is not a repair recommendation or expected training improvement.</p>")
    dist = diagnostics.get("distribution", {})
    dist = dist if isinstance(dist, Mapping) else {}
    sensitivity = dist.get("mmd_sensitivity", {})
    sensitivity = sensitivity if isinstance(sensitivity, Mapping) else {}
    content.append(_html_table(("field", "value"), [{"field": key, "value": sensitivity.get(key)} for key in ("before", "after", "change", "affected_genes")]))
    content.extend(["<p>Change = after − before; negative unbiased MMD values remain valid. Requested shifts, achieved shifts, and clipping are in diagnostics.json.</p></section><section id=task><h2>Task-specific evidence</h2>"])
    if task == "T3":
        response = diagnostics.get("response", {})
        response = response if isinstance(response, Mapping) else {}
        summary = response.get("summary", {})
        summary = summary if isinstance(summary, Mapping) else {}
        content.extend(["<h3>WT-relative response <span class=badge>A</span></h3><p>The official severity gate distinguishes accepted response, source-guarded response, inversion/unrelated response, and undefined truth. The n_sign_reversals count covers opposite tolerance-based signs across all supplied genes, including genes without a target DE call; it is not a count of meaningful DE response reversals. The source no_predicted_response low across-gene variation guard compares standard deviations of deltas and can trigger for a constant nonzero uniform offset; it does not establish that every literal delta is zero. A failure sentinel is not an amplitude estimate. Per-gene residuals do not sum to log severity.</p>",
                        _html_table(("field", "value"), [{"field": key, "value": summary.get(key)} for key in ("status", "n_truth_de", "n_sign_reversals", "beta", "official_logbeta", "official_r2", "magnitude_median_ratio", "n_magnitude_eligible", "reported_logbeta_origin")]),
                        _html_table(("gene", "truth_delta", "pred_delta", "literal_sign_status", "magnitude_ratio", "magnitude_status"), _top(response.get("genes"), "direct_error", 10))])
    if task in ("T2", "T3"):
        content.append("<h3>Spatial form, size, and local organization</h3><div class=notice>Cells are unpaired; source frames arbitrary. D2 is reflection blind. SW/Dice retain a source-frame caveat and have no adopted priority thresholds. Laterality is unassessed. Local organization is measured but has no calibrated follow-up threshold.</div>")
        spatial = diagnostics.get("spatial", {})
        spatial = spatial if isinstance(spatial, Mapping) else {}
        for family_name in ("shape", "growth", "local"):
            family = spatial.get(family_name, {})
            if isinstance(family, Mapping):
                for section_name, section in family.items():
                    if not isinstance(section, Mapping):
                        continue
                    content.extend([f"<article><h4>{esc(family_name)} / {esc(section_name)}</h4><p>{esc(section.get('definition', 'Undefined definition'))}</p>",
                                    _html_table(("field", "value"), [{"field": key, "value": "not part of the official T3 panel (extra diagnostic only)" if key == "official_value" and task == "T3" and section.get("official_panel_includes_metric") is False else section.get(key)} for key in ("status", "official_value", "replayed_value", "reconciliation_abs_error")])])
                    rows = section.get("rows", [])
                    if isinstance(rows, list) and rows:
                        headers = tuple(dict.fromkeys(key for row in rows[:5] if isinstance(row, Mapping) for key in row))
                        content.append(_html_table(headers, rows[:5]))
                    if section_name == "neighborhood_mmd":
                        content.append("<h5>Signed anchor kernel terms</h5><p>Strongest absolute terms, ordered within this table only. Finite-sample terms on both sides may be negative. They are not causal per-cell errors, anatomical matches, or follow-up priorities. Coordinates remain in each side's own unregistered frame.</p>")
                        content.append(_html_table(("side", "row_position", "original_row_position", "obs_name", "coordinates", "allocation"), _top(section.get("anchor_rows", []), "allocation", 10)))
                        allocation = section.get("anchor_allocation", {})
                        if isinstance(allocation, Mapping):
                            content.append(_html_table(("field", "value"), [{"field": key, "value": allocation.get(key)} for key in ("sum", "reconciliation_abs_error", "tolerance", "passed", "effective_seed", "n_prediction_anchors", "n_truth_anchors")]))
                    content.append("<ul>" + "".join(f"<li>{esc(item)}</li>" for item in section.get("limitations", [])) + "</ul></article>")
        content.append("<p>The official count_log_ratio uses selected counts and may become zero when both files reach the cap. Original counts are metadata only; these comparisons do not establish biological growth or proliferation.</p>")
    if task == "T1":
        content.append("<p>T1 has no spatial or WT-response task section. Shared gene, composition, pair, and distribution diagnostics are shown above.</p>")
    content.append("</section><section id=figures><h2>Measured local figures</h2>")
    for figure in diagnostics.get("figures", []):
        if isinstance(figure, Mapping):
            content.append(f"<figure><img loading=lazy src='{esc(figure['path'])}' alt='{esc(figure['caption'])}'><figcaption><span class=badge>{esc(figure['evidence_class'])}</span> {esc(figure['status'])}: {esc(figure['caption'])}</figcaption></figure>")
    content.extend(["</section><section id=warnings><h2>Warnings and limitations</h2><ul>",
                    *[f"<li>{esc(warning)}</li>" for warning in _warning_strings(diagnostics.get("warnings"))],
                    "</ul><p>A: direct/algebraic selected-sample evidence; B: sampled evidence; C: noncausal sensitivity; D: descriptive companion values. Undefined values are null, never implied zero. Selected-array memory limits do not bound whole-process memory. Repeated seed behavior is sampling evidence, not biological uncertainty.</p></section>",
                    "<footer>Local static report · no remote assets · full ledgers in diagnostics.json and CSV files · provenance in provenance.json</footer>"])
    css = "body{font-family:system-ui,sans-serif;color:#192636;max-width:1120px;margin:0 auto;padding:1.5rem;line-height:1.55}h1,h2,h3{line-height:1.2}header{padding:1rem 0;border-bottom:3px solid #335b9a}.eyebrow,.muted{color:#52647b}.eyebrow{font-size:.85rem}nav{display:flex;flex-wrap:wrap;gap:1rem;padding:1rem 0}a{color:#244f89}section{margin:2.5rem 0}.notice{background:#f0f4fa;border-left:4px solid #335b9a;padding:1rem}.finding{border:1px solid #c7d3e1;border-radius:6px;padding:1rem;margin:1rem 0}.badge{font-size:.75rem;background:#e6edf6;border-radius:3px;padding:.15rem .4rem;white-space:nowrap}table{border-collapse:collapse;width:100%;font-size:.87rem}th,td{text-align:left;padding:.55rem;border-bottom:1px solid #d6e0ec;vertical-align:top;overflow-wrap:anywhere}th{background:#eef3f8}.table-wrap{overflow:auto;margin:1rem 0}figure{margin:1.5rem 0;border:1px solid #d6e0ec;padding:.75rem}img{max-width:100%;height:auto;display:block;margin:auto}figcaption,footer{font-size:.85rem;color:#52647b}details{margin:1rem 0}summary{cursor:pointer;font-weight:600}footer{border-top:1px solid #c7d3e1;padding:1rem 0}@media print{nav{display:none}section,figure{break-inside:avoid}body{padding:0}}"
    return f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'><title>{esc(title)}</title><style>{css}</style></head><body>" + "".join(content) + "</body></html>"


def _spatial_csv_rows(spatial: Any) -> list[dict]:
    if not isinstance(spatial, Mapping):
        return []
    rows = []
    for family_name in ("shape", "growth", "local"):
        family = spatial.get(family_name, {})
        if isinstance(family, Mapping):
            for name, section in family.items():
                if isinstance(section, Mapping):
                    for row in section.get("rows", []):
                        if isinstance(row, Mapping):
                            rows.append({"diagnostic": f"{family_name}.{name}", "status": section.get("status"), **row})
    return rows


def _source_paths(value: Any, *, key: str = "") -> Iterable[Path]:
    """Find provenance paths conservatively, to avoid overwriting an input tree."""
    if isinstance(value, Mapping):
        for child_key, child in value.items():
            yield from _source_paths(child, key=str(child_key).lower())
    elif isinstance(value, (list, tuple)):
        for child in value:
            yield from _source_paths(child, key=key)
    elif isinstance(value, str) and any(token in key for token in ("path", "file", "input")):
        yield Path(value).expanduser()


def _validate_destination(bundle: Any, out: str | Path) -> Path:
    requested = Path(out).expanduser()
    if not requested.name or requested.name in (".", ".."):
        raise ValueError("out must name a new report directory")
    if requested.parent.is_symlink():
        raise ValueError("report output parent must not be a symlink")
    try:
        parent = requested.parent.resolve(strict=True)
    except FileNotFoundError as exc:
        raise ValueError("report output parent must already exist") from exc
    if not parent.is_dir() or parent.is_symlink():
        raise ValueError("report output parent must be a real directory")
    final = parent / requested.name
    try:
        final.lstat()
    except FileNotFoundError:
        pass
    else:
        raise FileExistsError(f"refusing to replace existing report destination: {final}")

    provenance = getattr(bundle, "provenance", {})
    for source in _source_paths(provenance):
        source = source.resolve(strict=False)
        if source == final or source in final.parents or final in source.parents:
            raise ValueError(f"report destination overlaps a source path: {source}")
    return final


def _write_text(path: Path, contents: str) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(contents)
        handle.flush()
        os.fsync(handle.fileno())


def write_report(bundle: Any, scored: Mapping[str, Any], diagnostics: Mapping[str, Any], out: str | Path) -> Path:
    """Write report files into a newly-created directory and return its path.

    A sibling staging directory makes partial reports invisible. The destination
    must be absent; reports are never appended to, overwritten, or written inside
    a provenance-declared source path.
    """
    final = _validate_destination(bundle, out)
    notes: dict[str, int] = {}
    safe_diagnostics = _json_safe(diagnostics, path="diagnostics", notes=notes)
    safe_scored = _json_safe(scored, path="scored", notes=notes)
    safe_provenance = _json_safe(getattr(bundle, "provenance", {}), path="bundle.provenance", notes=notes)
    if not isinstance(safe_diagnostics, dict):
        raise TypeError("diagnostics must be a mapping")
    if not isinstance(safe_scored, dict):
        raise TypeError("scored must be a mapping")
    # Preserve supplied official scalar values alongside diagnostic ledgers for
    # machine consumers.  The report does not synthesize or compare them.
    safe_diagnostics["official_metrics"] = safe_scored.get("metrics", {})
    task = str(safe_diagnostics.get("task", getattr(bundle, "task", safe_scored.get("provenance", {}).get("task", "T1"))))
    if task not in ("T1", "T2", "T3"):
        raise ValueError("report task must be T1, T2, or T3")
    safe_diagnostics["task"] = task
    safe_diagnostics["schema_version"] = "vec-scorelens-diagnostics-v2"
    if not isinstance(safe_diagnostics.get("scope"), Mapping):
        inputs = safe_provenance.get("inputs", {}) if isinstance(safe_provenance, Mapping) else {}
        selected, original = {}, {}
        for role, source in inputs.items() if isinstance(inputs, Mapping) else []:
            if isinstance(source, Mapping):
                label = "wt" if task == "T3" and role == "reference" else role
                for destination, field in ((selected, "selected_shape"), (original, "original_shape")):
                    shape = source.get(field)
                    if isinstance(shape, list) and shape:
                        destination[label] = shape[0]
        safe_diagnostics["scope"] = {"assessment": "sampled_subset_only", "gene_panel": "supplied_ordered_panel",
                                     "contrast_role": "wt" if task == "T3" else "reference",
                                     "selected_counts": selected, "original_counts": original}
    if not isinstance(safe_diagnostics.get("priorities"), Mapping):
        from .priority import prioritize
        safe_diagnostics["priorities"] = prioritize(safe_diagnostics, task=task)
    input_warnings = safe_provenance.get("warnings", []) if isinstance(safe_provenance, Mapping) else []
    merged_warnings = _warning_strings(safe_diagnostics.get("warnings"))
    merged_warnings.extend(_warning_strings(input_warnings))
    for section in ("response", "spatial"):
        block = safe_diagnostics.get(section, {})
        if isinstance(block, Mapping):
            merged_warnings.extend(_warning_strings(block.get("warnings")))
    merged_warnings.extend(
        f"{count} non-finite value(s) at {location} were written as null."
        for location, count in sorted(notes.items())
    )
    if notes or _contains_null(safe_diagnostics) or _contains_null(safe_scored):
        merged_warnings.append(
            "Null values in this report mean undefined, not applicable, or non-finite input values; see each diagnostic definition and limitation for its reason."
        )
    if merged_warnings:
        safe_diagnostics["warnings"] = list(dict.fromkeys(merged_warnings))

    stage = final.parent / f".{final.name}.stage-{uuid.uuid4().hex}"
    try:
        stage.mkdir(mode=0o700)
        from .figures import write_figures
        safe_diagnostics["figures"] = write_figures(bundle, safe_diagnostics, stage / "figures", task=task,
                                                    seed=int(safe_scored.get("provenance", {}).get("seed", 0)))
        unavailable_figures = [figure["path"] for figure in safe_diagnostics["figures"] if figure["status"] != "defined"]
        if unavailable_figures:
            safe_diagnostics.setdefault("warnings", []).append("Measured figure evidence unavailable: " + ", ".join(unavailable_figures) + ". Missing values are shown explicitly, without invented measurements.")
        summary = _summary(bundle, safe_scored, safe_diagnostics)
        _write_text(stage / "summary.md", summary)
        _write_text(stage / "diagnostics.json", json.dumps(safe_diagnostics, allow_nan=False, indent=2, sort_keys=True) + "\n")
        provenance = {
            "bundle": safe_provenance,
            "scorer": safe_scored.get("provenance", {}),
            "report": {"format": "vec-scorelens-private-report-v2", "task": task, "selected_subset_only": True,
                       "priority_profile": safe_diagnostics["priorities"].get("profile", {}).get("id"),
                       "figures": safe_diagnostics["figures"]},
        }
        _write_text(stage / "provenance.json", json.dumps(provenance, allow_nan=False, indent=2, sort_keys=True) + "\n")
        _write_csv(stage / "gene_errors.csv", safe_diagnostics.get("genes", []), _GENE_COLUMNS)
        _write_csv(stage / "composition_errors.csv", safe_diagnostics.get("composition", []), _COMPOSITION_COLUMNS)
        _write_csv(stage / "covariance_errors.csv", safe_diagnostics.get("pairs", []), _PAIR_COLUMNS)
        if task == "T3":
            response = safe_diagnostics.get("response", {})
            _write_csv(stage / "response_errors.csv", response.get("genes", []) if isinstance(response, Mapping) else [],
                       ("gene", "truth_delta", "pred_delta", "direct_error", "truth_de_eligible", "literal_sign_status", "magnitude_ratio", "magnitude_status", "slope_residual"))
        if task in ("T2", "T3"):
            _write_csv(stage / "spatial_errors.csv", _spatial_csv_rows(safe_diagnostics.get("spatial")), ("diagnostic", "status"))
            spatial = safe_diagnostics.get("spatial", {})
            local = spatial.get("local", {}) if isinstance(spatial, Mapping) else {}
            neighborhood = local.get("neighborhood_mmd", {}) if isinstance(local, Mapping) else {}
            if isinstance(neighborhood, Mapping) and "anchor_rows" in neighborhood:
                _write_csv(stage / "neighborhood_anchors.csv", neighborhood["anchor_rows"],
                           ("side", "row_position", "original_row_position", "obs_name", "coordinates", "allocation"))
        _write_text(stage / "report.html", _html_report(bundle, safe_scored, safe_diagnostics))
        os.rename(stage, final)
    except BaseException:
        if stage.exists() or stage.is_symlink():
            shutil.rmtree(stage, ignore_errors=True)
        raise
    return final
