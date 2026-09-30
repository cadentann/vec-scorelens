# Limitations and appropriate use

ScoreLens is a local diagnostic for user-supplied files. It is not an official submission validator, a hidden-target service, a biological-health assessment, or a recommendation engine. Success only means the bounded local computation completed under its input policies.

## Sampling, counts, and resources

- Official metrics and all selected-array evidence apply only to independently selected rows. The report distinguishes original file counts from selected counts; an original-count ratio is D-class descriptive metadata, not an official scalar or a full-file expression assessment.
- T1 retains all supplied genes up to 40,000. T2/T3 limit supplied panels to 2,000 genes as a companion resource policy because the spatial scorer forms across-gene neighborhoods. Neither bound is an official Challenge rule, and the dense-memory limit is not a total-process peak-memory guarantee.
- Target `celltype` metadata is capped at 1,000,000 values as a ScoreLens resource policy. Sparse matrices with malformed or excessive duplicate stored entries are rejected rather than normalized silently. These are supported-input limits, not statements about an input's biology or the Challenge schema.
- Sampling can omit rare populations, weak DE signals, spatial structures, or nonfinite values outside the selected rows. A report cannot certify full-file form, size, density, or biological state.

## Interpretation boundaries

- A direct diagnostic is exact only for its declared computation on selected data. B-class pair/distance outputs are sampled; C-class MMD output is a target-aware intervention; D-class summaries are descriptive. Only explicitly reconciled terms are additive allocations of the named formula. No diagnostic establishes biological causality.
- Frozen-probe proportions depend on target labels, selected cells, and probe domain. Low predicted mass is not biological absence.
- Local DE calls, ranked slots, and literal mean-change signs have different definitions. Mean log-expression deltas are not raw-count log fold changes. T3 response values are WT-relative transformed-expression deltas, not raw-count logFC.
- The upstream no-response standard-deviation guard can trigger for a constant nonzero predicted offset. Its branch labels a scorer guard condition, not proof that the predicted response is absent or biologically null.
- The T3 sign-reversal count includes all supplied genes with opposite tolerance-based signs, including non-DE genes with tiny changes. It is a descriptive count, not a count of significant or biologically meaningful response reversals. Per-gene magnitude ratios may remain defined when the global severity guard fails; the priority profile still excludes response-magnitude alerts in that branch.
- Variogram/covariance diagnostics do not identify regulation or causation. T3 residuals do not identify a perturbation mechanism or decompose severity.
- The variogram ledger samples ordered gene pairs. It does not localize every possible pair in large T1 panels.
- Spatial clouds are not assumed registered or cell-matched. ScoreLens makes no biological registration, anatomical identity, or laterality inference. D2 is reflection-blind; SW and Dice are frame-sensitive contextual diagnostics. A low Dice voxel-resolution value weakens interpretation; it does not prove a shape defect.
- The pinned upstream spatial implementation has an open source-frame/rotation concern. Therefore SW/Dice cannot independently trigger a priority finding, including for apparently identical clouds.
- T2/T3 reject a selected prediction or target cloud with zero spatial extent before producing a report. This is a ScoreLens supported-input limit: upstream can compute some other spatial terms for a collapsed cloud, but its scale term is undefined. ScoreLens does not currently report those partial diagnostics.

## Data integrity and privacy

- Reports may contain supplied paths, fingerprints, selected row identities, and data-derived tables/figures. Treat private input reports as private.
- Do not use held-out measured data, protected public-stage copies, or leaderboard responses to reconstruct target properties. A locally created pseudo-target is not organizer validation.
- The tool does not verify normalization, biological provenance, authorization, laterality, coordinate frame, or external-data eligibility. The user remains responsible for those decisions and for Challenge/data terms.

See [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md) for concrete implementation and upstream caveats.
