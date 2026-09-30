# Methods and diagnostic conventions

## Pipeline and scope

ScoreLens reads the supplied AnnData files in bounded, read-only form, deterministically selects rows, and converts only selected expression arrays to float32. It runs the task-specific public `veckit` v2 entry point on those arrays and preserves its scalar dictionary. ScoreLens does not reimplement authoritative scalar formulas or translate scalar scales into a cross-task score.

The source role is `reference` for T1/T2 and `wt` for T3. For T2/T3, selected `obsm['spatial_3D']` coordinates are read with the same row indices as expression. No pointwise correspondence, registration, tissue identity, or laterality is assumed between clouds.

Reports use four visibly separate evidence classes:

| Class | Meaning | Examples |
| --- | --- | --- |
| A — direct | Deterministic selected-array measurement or algebraic allocation. | Pseudobulk squared allocations, frozen-probe JSD terms, WT-relative means. |
| B — sampled | Seeded sampled measurement, reproducible but not a census. | Ordered-pair variogram ledger; sampled distance summaries. |
| C — sensitivity | Defined artificial intervention with same-sample before/after measurement. | Nonnegative-clipped MMD mean-shift replay. |
| D — descriptive | Companion statistic, not an official additive contribution or calibrated biological inference. | Covariance residuals, original count ratio, selected-cloud summaries. |

“Direct” means exact for its stated selected-array computation subject to floating-point reconciliation tolerance. It never means exact biology or full-file scoring.

## Shared expression and composition evidence

For each gene, ScoreLens exposes selected prediction, target, and contrast means; contrast-relative changes; pseudobulk squared residual allocations; and variance summaries. When defined, squared allocations sum to the official pseudobulk relative-error square. Mean transformed-expression shifts are not conventional fold changes.

DE reporting intentionally separates: (1) upstream target/prediction DE calls and official target-sized ranked slots, (2) local significant-set hit/miss/prediction-only calls, and (3) literal signs of contrast-relative mean changes. A literal-sign tally does not reproduce the official partial rank correlation. Guards for no change, degenerate truth, and undefined quantities are reported rather than coerced to zero.

The frozen target probe is fitted only to selected target cells and their `celltype` labels. It yields target class fractions, prediction soft fractions, deltas, and upstream count-versus-fraction-smoothed JSD terms. Prediction labels are ignored. Probe output is limited to the target label domain and is not evidence of biological absence.

The pair ledger replays the selected-cell draw and 20,000 ordered gene-pair draws used for the variogram view, retains duplicate occurrences, and removes self-pairs. It is computed in bounded blocks and reconciles against the corresponding unrounded upstream function within a numeric tolerance. Covariance is a D-class companion residual: a variogram is not pure covariance and neither supports gene-regulation claims.

The MMD sensitivity selects up to five high-pseudobulk-error genes, applies a nonnegative-clipped target mean-shift to a copy of the selected prediction, and recomputes upstream MMD with the same target and seed. Before/after/change, requested and achieved shifts, clipping, and genes are recorded. It is C-class, noncausal, nonadditive, target-aware, and may worsen MMD.

## T2 spatial evidence

T2 records selected coordinate counts/dimensions, radius and count summaries, sampled D2 distance-distribution evidence, shape/occupancy context, and local spatial-expression summaries. Coordinates are centered or canonicalized only for the named comparison; these are not registration operations and do not pair cells. The Moran ledger selects the target-variance panel once, then records centered per-gene products whose sum reconciles to its Pearson agreement when defined; raw per-gene I differences remain descriptive and are not location matches or biological explanations. The neighborhood-MMD ledger also records signed, two-sided prediction and target anchor kernel terms. Both sides are required for the terms to sum to the native MMD-u scalar within the recorded native-dtype reconciliation tolerance; individual anchors can be negative and are neither matched cells nor causal errors.

Spatial figures use each selected cloud's own supplied coordinate frame (or separately centered native projection). They are not matched-cell, registered-anatomy, or laterality visualizations.

## T3 matched-WT response evidence

T3 records target DE eligibility; per-gene `target − WT` and `prediction − WT` mean transformed-expression deltas; slope numerator, denominator, and beta when defined; response residuals; R²/guard branch; and the unchanged pinned severity scalar. It separates complete no-response, inverted/unrelated response, too-few-DE, and degenerate-truth branches. The upstream no-response standard-deviation guard can also trigger for a constant nonzero predicted offset, so its guard label is not evidence that a response is literally absent. Per-gene response residuals do not decompose the log severity metric and are not raw-count logFC values.

## Priority ordering

The priority block ranks diagnostic follow-up only. It is not a biological-health verdict, a calibrated false-discovery rate, an official metric rank, or model-repair advice. Rules are deterministic and versioned; every triggered row names its evidence, operation, threshold, units, class, selected scope, limitation, and suggested check. The implemented order is literal response reversal, mean residual, eligible T3 response magnitude, frozen-probe composition gap, selected variance ratio, selected spatial size, and original-count metadata. Input/scope warnings are shown separately; shape/local organization, sampled pair/covariance evidence, and target-aware sensitivity have no adopted priority threshold. Rules cannot compare unrelated raw metric scales.

In particular, slicewise Wasserstein (SW) and Dice evidence cannot independently label a file a failure: upstream coordinate-frame ambiguity can affect even identical point clouds. Their values are diagnostic context, not a biological-registration or laterality test.
