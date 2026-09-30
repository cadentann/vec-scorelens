"""T1 acceptance controls with independently sampled, uncentered populations.

These deliberately differ from the exact, centered synthetic fixtures.  They
exercise the public T1 scorer and companion on ordinary finite samples, where
sample means and mixtures genuinely differ.  The simple priority thresholds
are exposed heuristics, not false-discovery-rate estimates.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import pytest

from vec_scorelens.diagnostics import diagnose
from vec_scorelens.io import Bundle
from vec_scorelens.scorer import run_official


_SMALL_GENES = 128
_FULL_GENES = 32_285
_SMALL_CELLS = 192
_MODULES = (tuple(range(32_000, 32_008)), tuple(range(32_008, 32_016)))


def _readonly(matrix: np.ndarray) -> np.ndarray:
    result = np.asarray(matrix, dtype=np.float32)
    result.setflags(write=False)
    return result


def _mixture_counts(rng: np.random.Generator, cells: int, states: int = 4) -> np.ndarray:
    """Random, non-identical mixture counts rather than fixed fixture counts."""
    # A high-but-finite Dirichlet concentration retains ordinary sampling
    # variation without turning a small control into a rare-class stress test.
    counts = rng.multinomial(cells, rng.dirichlet(np.full(states, 80.0)))
    if np.any(counts == 0):  # Defensive only; preserves every probe class.
        return _mixture_counts(rng, cells, states)
    return counts


def _ordinary_draw(
    rng: np.random.Generator, means: np.ndarray, counts: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Draw uncentered cells with cell-scale and gene-level sampling noise."""
    rows: list[np.ndarray] = []
    labels: list[str] = []
    gene_loading = np.linspace(-0.18, 0.22, means.shape[1])
    for state, count in enumerate(counts):
        cell_scale = rng.normal(0.0, 0.16, (int(count), 1))
        residual = rng.normal(0.0, 0.18, (int(count), means.shape[1]))
        matrix = means[state] + cell_scale * gene_loading + residual
        # The construction's positive baseline makes this clipping inactive in
        # ordinary draws; it simply honors the adapter's nonnegative contract.
        rows.append(np.maximum(matrix, 0.0))
        labels.extend([f"state_{state}"] * int(count))
    return np.vstack(rows), np.asarray(labels, dtype=str)


def _ordinary_bundle(seed: int) -> tuple[Bundle, dict[str, np.ndarray]]:
    """Three independently sampled populations with an actual target response."""
    rng = np.random.default_rng(seed)
    baseline = rng.uniform(2.7, 3.7, _SMALL_GENES)
    state_means = np.broadcast_to(baseline, (4, _SMALL_GENES)).copy()
    for state in range(4):
        state_means[state, state * 12 : (state + 1) * 12] += 1.4
    response = np.zeros(_SMALL_GENES)
    response[48:64] = np.linspace(0.24, 0.42, 16)
    response[64:80] = -np.linspace(0.20, 0.38, 16)

    reference_counts = _mixture_counts(rng, _SMALL_CELLS)
    target_counts = _mixture_counts(rng, _SMALL_CELLS)
    prediction_counts = _mixture_counts(rng, _SMALL_CELLS)
    reference, _ = _ordinary_draw(rng, state_means, reference_counts)
    target, target_labels = _ordinary_draw(rng, state_means + response, target_counts)
    prediction, _ = _ordinary_draw(rng, state_means + response, prediction_counts)
    bundle = Bundle(
        prediction=_readonly(prediction), target=_readonly(target), reference=_readonly(reference),
        genes=[f"ordinary_gene_{index:03d}" for index in range(_SMALL_GENES)],
        target_labels=target_labels,
        provenance={"scope": "independently sampled uncentered synthetic control", "synthetic_only": True},
    )
    return bundle, {
        "prediction_counts": prediction_counts,
        "target_counts": target_counts,
        "reference_counts": reference_counts,
    }


@lru_cache(maxsize=30)
def _ordinary_evaluation(seed: int) -> tuple[Bundle, dict[str, np.ndarray], dict, dict]:
    bundle, counts = _ordinary_bundle(seed)
    scored = run_official(bundle, seed=seed)
    return bundle, counts, scored, diagnose(bundle, scored, seed=seed)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_independently_sampled_uncentered_controls_are_not_perfect_means(seed: int) -> None:
    """Use the real scorer/diagnostics without exact centering or copied cells."""
    bundle, counts, scored, diagnostic = _ordinary_evaluation(seed)
    assert len({tuple(values) for values in counts.values()}) > 1
    assert max(row["absolute_error"] for row in diagnostic["genes"]) > 0.0
    assert any(abs(row["proportion_error"]) > 0.0 for row in diagnostic["composition"])
    assert 0.5 <= diagnostic["distribution"]["variance_ratio"] <= 2.0
    assert np.isfinite(list(scored["metrics"].values())).all()
    # The target's response is present but neither prediction nor reference is
    # an exact realization of it; this is not an all-zero-error control.
    assert diagnostic["gene_summary"]["pb_rel_err"] > 0.0
    assert diagnostic["gene_summary"]["n_local_false_positives"] >= 0


def test_priority_heuristic_has_no_alerts_on_held_out_ordinary_controls() -> None:
    """Exercise the actual exposed heuristic, without claiming FDR control."""
    observed: list[dict] = []
    for seed in (0, 1, 2):
        _, _, _, diagnostic = _ordinary_evaluation(seed)
        priority = diagnostic["priorities"]
        assert priority["profile"]["calibration"].startswith("heuristic")
        assert {entry["rule_id"] for entry in priority["evaluations"]} >= {
            "mean_residual", "composition_gap", "variance_ratio"
        }
        observed.append({
            "assessment": priority["assessment"],
            "rules": [finding["rule_id"] for finding in priority["findings"]],
        })
    # This checks the threshold actually exposed to users.  Three ordinary
    # controls support only this bounded regression claim, not a calibration or
    # FDR statement.
    assert observed == [
        {"assessment": "within_listed_thresholds", "rules": []},
        {"assessment": "within_listed_thresholds", "rules": []},
        {"assessment": "within_listed_thresholds", "rules": []},
    ]


def test_thirty_independent_uncentered_controls_screen_exposed_priority_rules() -> None:
    """A larger finite-sample screen; it is not biological/FDR calibration.

    The assertions intentionally preserve any observed alerts as test evidence
    rather than selecting a favorable subset or tuning the held-out generator.
    The companion closeout record reports the measured count, rate, and error
    ranges from this exact seed sequence.
    """
    assessed: list[str] = []
    triggered_rules: list[str] = []
    max_gene_errors: list[float] = []
    max_composition_gaps: list[float] = []
    variance_ratios: list[float] = []
    for seed in range(30):
        _, _, _, diagnostic = _ordinary_evaluation(seed)
        priority = diagnostic["priorities"]
        assessed.append(priority["assessment"])
        triggered_rules.extend(finding["rule_id"] for finding in priority["findings"])
        max_gene_errors.append(max(float(row["absolute_error"]) for row in diagnostic["genes"]))
        max_composition_gaps.append(max(abs(float(row["proportion_error"])) for row in diagnostic["composition"]))
        variance_ratios.append(float(diagnostic["distribution"]["variance_ratio"]))

    assert len(assessed) == 30
    assert set(assessed) <= {"within_listed_thresholds", "thresholds_exceeded"}
    assert set(triggered_rules) <= {"opposite_sign", "mean_residual", "composition_gap", "variance_ratio"}
    # Confirm that this remains an ordinary, non-perfect sample screen even
    # where no threshold triggers; numerical ranges are recorded privately.
    assert min(max_gene_errors) > 0.0
    assert min(max_composition_gaps) > 0.0
    assert all(ratio > 0.0 for ratio in variance_ratios)


def _full_width_module_bundle(seed: int = 914) -> Bundle:
    """A 32,285-gene, 64-cell test-only panel with two prespecified modules."""
    rng = np.random.default_rng(seed)
    genes = [f"fullwidth_gene_{index:05d}" for index in range(_FULL_GENES)]
    baseline = rng.uniform(2.2, 3.0, _FULL_GENES)
    state_offset = np.zeros((_FULL_GENES,), dtype=float)
    state_offset[:32] = 0.9

    def draw(*, target: bool) -> np.ndarray:
        labels = np.repeat(np.arange(2), 32)
        matrix = baseline + rng.normal(0.0, 0.12, (64, _FULL_GENES))
        matrix += (labels[:, None] == 0) * state_offset
        if target:
            matrix[:, 96:112] += 0.20
        for module in _MODULES:
            matrix[:, module] += rng.normal(0.0, 0.60, (64, 1))
        return np.maximum(matrix, 0.0)

    reference = draw(target=False)
    target = draw(target=True)
    prediction = draw(target=True)
    # Test-only planted covariance loss: each module's marginals remain, but
    # its within-module dependence is independently permuted by gene.
    for module in _MODULES:
        for gene in module:
            prediction[:, gene] = prediction[rng.permutation(64), gene]
    return Bundle(
        prediction=_readonly(prediction), target=_readonly(target), reference=_readonly(reference),
        genes=genes, target_labels=np.repeat(np.array(["type_a", "type_b"]), 32),
        provenance={"scope": "full-width selected synthetic module coverage control", "synthetic_only": True},
    )


def test_full_width_pair_sampling_reports_coverage_without_module_discovery_claim() -> None:
    """20k random ordered pairs cannot establish coverage of tiny modules."""
    bundle = _full_width_module_bundle()
    truth_dependence = []
    prediction_dependence = []
    for module in _MODULES:
        triangle = np.triu_indices(len(module), k=1)
        truth_dependence.extend(np.abs(np.corrcoef(bundle.target[:, module], rowvar=False)[triangle]))
        prediction_dependence.extend(np.abs(np.corrcoef(bundle.prediction[:, module], rowvar=False)[triangle]))
    # Establish that the test-only construction has a meaningful dependence
    # loss before testing its sampling coverage.  This is not a product module
    # score and does not ask diagnostics to infer these prespecified sets.
    assert np.median(truth_dependence) > 0.80
    assert np.median(prediction_dependence) < 0.30
    scored = run_official(bundle, seed=914)
    diagnostic = diagnose(bundle, scored, seed=914)
    summary = diagnostic["pair_summary"]
    total = _FULL_GENES * (_FULL_GENES - 1)
    assert summary["total_possible_ordered_pairs"] == total
    assert summary["n_sampled_pairs"] <= 20_000
    assert summary["unique_ordered_pair_coverage"] == pytest.approx(summary["n_unique_ordered_pairs"] / total)
    assert summary["unique_ordered_pair_coverage"] < 3e-5

    observed_gene_coverage = summary.get(
        "observed_gene_coverage", summary["n_genes_observed_in_pairs"] / _FULL_GENES
    )
    assert observed_gene_coverage == pytest.approx(
        len({row["gene_i"] for row in diagnostic["pairs"]} | {row["gene_j"] for row in diagnostic["pairs"]}) / _FULL_GENES
    )

    module_pairs = {
        (bundle.genes[left], bundle.genes[right])
        for module in _MODULES
        for left in module
        for right in module
        if left != right
    }
    sampled = {(row["gene_i"], row["gene_j"]) for row in diagnostic["pairs"]}
    # This known-module check is test truth only.  It establishes absence of
    # sampled pair evidence for this seed; it neither asks the product to
    # discover a module nor turns unsampled corruption into a negative finding.
    assert not sampled.intersection(module_pairs)
    assert "expected_module_pair_hits" not in summary
    assert "module" not in summary["coverage_interpretation"].lower() or "not modules" in summary["coverage_interpretation"].lower()
