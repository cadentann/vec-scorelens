# Known limitations

- ScoreLens scores a bounded selected subset, not whole files. Original and selected counts are deliberately separate; neither count metadata nor a selected spatial cloud establishes whole-file form or expression quality.
- T1 supports up to 40,000 supplied genes. T2/T3 use a 2,000-gene spatial companion cap to bound 15-neighbor across-gene work. These are ScoreLens policies, not Challenge requirements.
- Target `celltype` metadata is capped at 1,000,000 values. Sparse inputs with excessive duplicate stored entries are unsupported rather than silently coalesced; both are resource/supported-input policies.
- T2/T3 require finite `obsm['spatial_3D']` coordinates aligned with selected expression rows. The first three columns are used. Coordinate frames are not registered, and no cell correspondence, anatomy, or laterality is inferred.
- The pinned upstream spatial scorer has known source-frame/rotation ambiguity. SW and Dice are displayed as contextual evidence but cannot independently cause a priority/failure label. D2 is reflection-blind; Dice depends on voxel resolution and nearest-neighbor spacing.
- T3 compares matched-WT mean transformed expression. Its delta and response ledger are not raw-count log fold changes, causal perturbation attribution, or a replacement severity metric.
- A constant nonzero prediction offset can meet the upstream no-response standard-deviation guard; that status is a scoring guard, not a biological absence finding.
- Frozen-probe composition depends on supplied target labels; prediction labels are ignored. A zero or low soft fraction is not a biological-absence conclusion.
- Variogram pairs are sampled, with duplicate draws retained; covariance is descriptive and not evidence of regulation.
- The pair table is not exhaustive all-gene-pair localization for large transcriptome panels.
- The MMD mean-shift is a deliberately target-aware sensitivity on the same selected target. It may worsen the metric and must not be used as a prediction edit or expected-improvement claim.
- Priority thresholds are disclosed diagnostic conventions. They are not calibrated FDRs, health scores, biological-normality statements, or training recommendations.
- Consumed HDF5 fields must use in-file hard links and storage. External links, soft links, external raw storage and virtual datasets are unsupported so scored values cannot silently come from unfingerprinted sidecars. Unused sidecars are not scanned or assessed.
- An abrupt interruption can leave a hidden `.<report-name>.stage-*` sibling with partial private data, though the final report directory remains absent. Inspect and remove that exact folder only when no run is active and its contents are no longer needed; never treat it as a successful report.
