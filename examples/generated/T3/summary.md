# VEC ScoreLens T3 private diagnostic report — selected subset

> Synthetic public layout example. JSON and CSV ledgers are omitted from this display copy; a full locally generated report includes them.

This report may identify supplied data. Keep it private and do not treat selected-subset results as full-file certification.

## Diagnostic follow-up priorities

Assessment: **thresholds exceeded**. Coverage: listed rules evaluated.
These explicit heuristic conventions order diagnostic follow-up. They are not statistically calibrated, have no guaranteed false-alarm rate or FDR, and do not establish biological health or expected model improvement.
Profile: scorelens\-diagnostic\-conventions\-v1. At most 6 findings are displayed in fixed domain order; all evaluated thresholds are recorded in diagnostics.json.

### 1. Meaningful literal response reversals

Evidence A; selected supplied cells. Triggered rows: 16.
Inspect named gene contrasts, input alignment, and normalization; operational target DE is not confirmed biological regulation.
| entity | measured_value | comparator | threshold | units |
| --- | --- | --- | --- | --- |
| synthetic\_gene\_119 | 1.03549 | opposite sign and absolute predicted delta &gt; | 0.25 | supplied transformed expression |
| synthetic\_gene\_087 | \-1.03548 | opposite sign and absolute predicted delta &gt; | 0.25 | supplied transformed expression |
| synthetic\_gene\_118 | 1.01613 | opposite sign and absolute predicted delta &gt; | 0.25 | supplied transformed expression |
| synthetic\_gene\_086 | \-1.01613 | opposite sign and absolute predicted delta &gt; | 0.25 | supplied transformed expression |
| synthetic\_gene\_085 | \-0.996775 | opposite sign and absolute predicted delta &gt; | 0.25 | supplied transformed expression |

### 2. Large selected\-sample mean residuals

Evidence A; selected supplied cells. Triggered rows: 32.
Inspect named gene means, input alignment, normalization, and the supplied target/reference contrast.
| entity | measured_value | comparator | threshold | units |
| --- | --- | --- | --- | --- |
| synthetic\_gene\_119 | 2.07097 | &gt; | 0.25 | supplied transformed expression |
| synthetic\_gene\_087 | 2.07097 | &gt; | 0.25 | supplied transformed expression |
| synthetic\_gene\_086 | 2.03226 | &gt; | 0.25 | supplied transformed expression |
| synthetic\_gene\_118 | 2.03226 | &gt; | 0.25 | supplied transformed expression |
| synthetic\_gene\_117 | 1.99355 | &gt; | 0.25 | supplied transformed expression |

Coverage exclusions:
- shape\_and\_local\_priority: Shape and local organization are measured separately but have no adopted priority thresholds; SW/Dice have a source\-frame caveat.
- laterality: No calibrated anatomical laterality diagnostic is computed.
- pair\_priority: Sampled variogram and descriptive covariance remain evidence tables without an adopted priority threshold.
- sensitivity\_priority: Target\-aware sensitivity is noncausal and does not determine follow\-up priority.
Within listed thresholds means only that evaluated listed rules did not trigger on supplied selected cells. Undefined evidence and unassessed domains remain explicit; this is not a full-file all-clear.

## Assessed input scope

The supplied ordered gene panel and independently selected cells are assessed. The contrast is the matched WT; equal cell counts do not create paired cells.

| input | original_cells | selected_cells |
| --- | --- | --- |
| prediction | 256 | 256 |
| target | 256 | 256 |
| wt | 256 | 256 |

## Official score output

The values below are the supplied official score output for the selected subset. They are different scales and definitions, so this report intentionally does not label any one metric the ‘weakest’. 
| metric | value |
| --- | --- |
| de\_score | 0.6885 |
| de\_direction | 0.8054 |
| severity\_slope | \-0.3839 |
| energy\_distance | 5.37151 |
| mmd\_u | 0.01005 |
| variogram | 0.02486 |
| pb\_rel\_err | 0.065 |
| library\_size\_ratio | 1.042 |
| variance\_ratio | 1 |
| composition\_JSD | 0 |
| pseudobulk\_pearson | 0.7316 |
| \_de\_raw | 0.7031 |
| \_de\_chance | 0.0469 |
| \_n\_up | 32 |
| \_n\_dn | 32 |
| \_slope\_r2 | 0.464 |
| d2\_shape | 0.00151 |
| sliced\_wasserstein | 0 |
| occupancy\_dice | 1 |
| scale\_log\_ratio | 0 |
| count\_log\_ratio | 0 |
| neighborhood\_mmd | 0.0125 |
| \_sw\_flip\_spread | 0.21287 |
| \_dice\_voxel\_over\_nn | 5.5 |

## Concrete gene follow-up

Genes below have the largest defined per-gene pseudobulk squared residuals in the selected data. Check their input alignment, normalization, and target/reference contrast before changing a model.
A — direct selected-sample evidence. `change_truth` is target mean minus matched WT mean; `change_pred` is prediction mean minus matched WT mean, in supplied transformed expression units.
| gene | change_truth | change_pred | absolute_error | squared_error | truth_de_direction | pred_local_de_direction | official_slot_hit |
| --- | --- | --- | --- | --- | --- | --- | --- |
| synthetic\_gene\_119 | \-1.03548 | 1.03549 | 2.07097 | 4.28891 | down | up | False |
| synthetic\_gene\_087 | 1.03548 | \-1.03548 | 2.07097 | 4.2889 | up | down | False |
| synthetic\_gene\_086 | 1.01613 | \-1.01613 | 2.03226 | 4.13009 | up | down | False |
| synthetic\_gene\_118 | \-1.01613 | 1.01613 | 2.03226 | 4.13008 | down | up | False |
| synthetic\_gene\_117 | \-0.996775 | 0.996775 | 1.99355 | 3.97424 | down | up | False |

Local significant-set hit/miss/false-positive fields and literal mean-log-expression signs are diagnostic aids. A false-positive means a prediction-only local DE call, not a known biological false response or an FDR conclusion; target nonsignificance does not mean no change. These fields are not the official target-sized ranked slots or official partial rank correlation. Rank-based direction can change with float32 roundoff near zero shifts; tiny literal sign differences are not evidence of a meaningful response.
Pseudobulk allocation reconciliation error: 3.04789e\-11. Guard state: none.
Target local DE calls: 32 up and 32 down; local signed matches 48, misses 16, and prediction-only calls 32. Target rank-slot hits with a wrong literal sign: 0.
Distribution variance-ratio evidence (selected supplied cells): 1. A low value is descriptive collapse evidence, not a biological conclusion.

## Concrete composition follow-up

Classes below have the largest defined target-versus-probe soft-proportion deltas. Inspect target labels and probe-domain fit; prediction labels are intentionally not used. Bounded sampling can omit rare classes, and a low predicted proportion is not biological absence.
| celltype | target_fraction | pred_soft_fraction | proportion_error | jsd_term |
| --- | --- | --- | --- | --- |
| state\_0 | 0.375 | 0.375294 | 0.000294417 | 4.15343e\-08 |
| state\_1 | 0.28125 | 0.280962 | \-0.000288486 | 5.34236e\-08 |
| state\_2 | 0.21875 | 0.21872 | \-3.03388e\-05 | 7.52191e\-10 |
| state\_3 | 0.125 | 0.125024 | 2.43038e\-05 | 8.87834e\-10 |
JSD term reconciliation error: 0.

## Sampled variogram and descriptive covariance follow-up

These rows aggregate sampled ordered pairs (duplicates retained as occurrences). Variogram errors are sampled fractional-moment errors; covariance residuals are descriptive companion columns, not a claim that variogram is pure covariance. Covariance uses population ddof=0 on the supplied transformed values, so its units follow those values rather than a biological scale.
| gene_i | gene_j | sampled_occurrences | variogram_squared_error | covariance_error |
| --- | --- | --- | --- | --- |
| synthetic\_gene\_087 | synthetic\_gene\_128 | 1 | 1.34836 | \-7.36357e\-10 |
| synthetic\_gene\_094 | synthetic\_gene\_119 | 1 | 1.34245 | \-1.18115e\-10 |
| synthetic\_gene\_127 | synthetic\_gene\_087 | 1 | 1.33188 | \-6.25103e\-10 |
| synthetic\_gene\_119 | synthetic\_gene\_122 | 1 | 1.30491 | 3.39041e\-10 |
| synthetic\_gene\_121 | synthetic\_gene\_085 | 1 | 1.30165 | \-7.22775e\-11 |
Sampled pair count: 19962; unique ordered pairs: 19213; reconciliation error: 2.47047e\-10.
Possible ordered nonself pairs: 249500; unique ordered-pair coverage fraction: 0.077006. Coverage is of ordered nonself pairs, not modules. Unsampled pairs have no diagnostic evidence; a small gene module may have no sampled within\-module pair even when every gene appears somewhere.

## MMD sensitivity (noncausal)

One nonnegative-clipped target mean-shift sensitivity was run: before 0.0100527, after 0.010922, change 0.000869369; affected genes: synthetic\_gene\_119, synthetic\_gene\_087, synthetic\_gene\_086, synthetic\_gene\_118, synthetic\_gene\_117.
The listed delta is after minus before. This target-aware, selected-and-evaluated-on-the-same-target artificial diagnostic is not a validated model correction, causal or additive attribution, percent contribution, or expected training improvement; it can worsen the score.

## WT-relative response

A — direct operational response evidence. WT matching, target DE eligibility, literal signs, magnitude ratios, and the guarded through-origin severity fit are separate concepts.
| field | value |
| --- | --- |
| status | accepted |
| n\_truth\_de | 64 |
| n\_sign\_reversals | 16 |
| beta | 0.681205 |
| official\_logbeta | \-0.383892 |
| official\_r2 | 0.464 |
| magnitude\_median\_ratio | 1 |
| n\_magnitude\_eligible | 48 |
| reported\_logbeta\_origin | accepted\_log\_beta |
| gene | truth_delta | pred_delta | literal_sign_status | magnitude_ratio | magnitude_status | slope_residual |
| --- | --- | --- | --- | --- | --- | --- |
| synthetic\_gene\_119 | \-1.03548 | 1.03549 | opposite | null (undefined/not applicable) | not\_interpretable\_opposite | 1.74086 |
| synthetic\_gene\_087 | 1.03548 | \-1.03548 | opposite | null (undefined/not applicable) | not\_interpretable\_opposite | \-1.74086 |
| synthetic\_gene\_086 | 1.01613 | \-1.01613 | opposite | null (undefined/not applicable) | not\_interpretable\_opposite | \-1.70832 |
| synthetic\_gene\_118 | \-1.01613 | 1.01613 | opposite | null (undefined/not applicable) | not\_interpretable\_opposite | 1.70832 |
| synthetic\_gene\_117 | \-0.996775 | 0.996775 | opposite | null (undefined/not applicable) | not\_interpretable\_opposite | 1.67578 |
The n_sign_reversals count covers opposite tolerance-based signs across all supplied genes, including genes without a target DE call; it is not a count of meaningful DE response reversals. A failure sentinel is not a measured log amplitude. The source no_predicted_response low across-gene variation guard compares standard deviations of all deltas, so a constant nonzero uniform offset can also trigger it; it does not mean every literal delta is zero. Source-guarded, inversion/unrelated response, and undefined truth branches must be read through the displayed gate. Gene slope residuals do not sum to log severity, and these selected-cell operational calls are not confirmed biological responses.

## Spatial form, size, counts, and local organization

Coordinates are selected with expression row identities. Cells are unpaired and coordinate frames arbitrary. D2 is reflection blind; no calibrated laterality claim is possible. SW/Dice retain a source-frame caveat and receive no priority threshold. Local organization is measured without a calibrated follow-up threshold.

### shape / d2

Wasserstein distance between independently sampled RMS\-normalized within\-cloud pair distances.
Status: finite; official: 0.00151, replay 0.00151106, reconciliation 8.67362e\-19.
- Independent sampled pairs give a nonzero identical\-cloud Monte Carlo floor; row\-order changes can change the sampled realization.
- Distance\-band contributions are not anatomical regions, homologous cell pairs, or causal attributions.
- D2 is reflection\-blind; finite distance distributions do not uniquely identify a shape.

### shape / sliced\_wasserstein

Mean projected sorted\-rank distance at the one globally selected upstream sign flip.
Status: provisional\_frame; official: 0, replay 1.15246e\-16, reconciliation 0.
- Pinned sliced Wasserstein and occupancy Dice have a PCA\-frame handedness defect: proper rotation or independent row permutation can change their scores, even for identical clouds. Treat these as provisional frame\-dependent diagnostics, not standalone failure priorities or laterality tests.
- Projection rank matching is not biological cell correspondence.
- Flip spread is not a calibrated uncertainty interval or sufficient test for eigenvalue degeneracy.

### shape / occupancy

Binary occupied\-voxel Dice on the pinned canonical grid; unmatched voxels allocate 1−Dice.
Status: provisional\_frame; official: 1, replay 1, reconciliation 0.
- Pinned sliced Wasserstein and occupancy Dice have a PCA\-frame handedness defect: proper rotation or independent row permutation can change their scores, even for identical clouds. Treat these as provisional frame\-dependent diagnostics, not standalone failure priorities or laterality tests.
- Voxels are metric\-grid regions, not registered anatomical regions.
- Voxel/NN ratio below roughly 10 risks measuring sampling density; clipped outliers occupy boundary voxels.
- Point counts are descriptive, not Dice weights.

### growth / scale

Signed log RMS\-radius ratio; shell contributions sum to each side's squared RMS, not to the log ratio.
Status: finite; official: 0, replay 0, reconciliation 0.
- RMS extent is not occupied volume, biological growth, or anatomical region error.
- Coordinate units must match; all coordinates here are selected\-scope.

### growth / count

Signed log ratio of supplied selected expression row counts; original\-file counts are separate metadata diagnostics.
Status: finite; official: 0, replay 0, reconciliation 0.
- Equal cell caps can hide original\-file count differences in the selected\-scope official metric.
- Observed cell counts do not by themselves establish proliferation or cell loss.

### local / neighborhood\_mmd

Truth\-fitted PCA multi\-kernel MMD\-u of spatial\-neighborhood pseudobulks; gene means are descriptive residuals.
Status: finite; official: 0.0125, replay 0.0125015, reconciliation 5.41148e\-08.
- Gene residuals are not additive MMD contributions; no gene\-wise reconciliation is claimed.
- MMD\-u can be negative, including for identical arrays; no clipping to zero.
- Neighborhoods are constructed after cell selection and may differ from full\-file neighborhoods.
- Constant neighborhood truth can make the truth\-fitted kernel undefined even when original expression varies.

Signed finite-sample two-sided anchor kernel allocations, ordered by absolute term magnitude within this table only. These are neither causal cell errors nor diagnostic priorities; coordinates remain in each side's own unregistered frame.
| side | row_position | coordinates | allocation |
| --- | --- | --- | --- |
| prediction | 4 | \[\-2.018634080886841, \-0.46226853132247925, 0.13211746513843536\] | 4.4093e\-05 |
| prediction | 86 | \[\-1.9920319318771362, \-0.6681209802627563, 0.1053728461265564\] | 4.32728e\-05 |
| prediction | 43 | \[\-2.237082004547119, \-0.642333984375, 0.057433776557445526\] | 4.23491e\-05 |
| prediction | 0 | \[\-2.055737257003784, \-0.4860828220844269, 0.3197552263736725\] | 4.22224e\-05 |
| prediction | 70 | \[\-2.046781063079834, \-0.30012285709381104, 0.30037230253219604\] | 4.12882e\-05 |
| prediction | 93 | \[\-2.4895896911621094, \-0.6355416178703308, \-0.15327399969100952\] | 4.12473e\-05 |
| prediction | 18 | \[\-2.202876567840576, \-0.7481556534767151, 0.16611990332603455\] | 4.1214e\-05 |
| prediction | 87 | \[\-2.2339892387390137, \-0.7293202877044678, 0.0011549079790711403\] | 4.11886e\-05 |
| prediction | 75 | \[\-2.0052576065063477, \-0.626953125, 0.4057144522666931\] | 4.11677e\-05 |
| prediction | 68 | \[\-2.154982805252075, \-0.7655008435249329, \-0.00936847273260355\] | 4.1141e\-05 |

### local / moran

Pearson agreement of Moran's I on genes selected once by truth variance; centered gene products allocate correlation.
Status: finite; official: not part of the official T3 panel (extra diagnostic only), replay 0.999828, reconciliation 1.24353e\-07.
- Autocorrelation strength agreement does not establish matching expression locations or laterality.
- Constant I profiles make correlation undefined; a small gene profile can give uninformative extreme correlations.
- Moran drops the first of 16 queried neighbors; duplicate\-coordinate ties may affect self\-exclusion.
Original prediction/target counts: 256 / 256; selected: 256 / 256. Original count log ratio: 0 (D — metadata only). The official count_log_ratio uses selected counts and can be zero after both files reach the cap. Do not infer biological growth or proliferation.

## Measured figures

![Per\-gene measured mean changes. Contrast is the supplied reference (T1/T2) or matched WT (T3). Named genes have the largest literal change residuals; these are not biological importance rankings.](figures/gene_changes.png)
A — defined: Per\-gene measured mean changes. Contrast is the supplied reference (T1/T2) or matched WT (T3). Named genes have the largest literal change residuals; these are not biological importance rankings.

![Up to 12 classes with largest measured soft\-fraction gaps. Prediction labels are ignored. A low predicted fraction is not biological absence; omitted classes remain explicit report warnings.](figures/composition.png)
A — defined: Up to 12 classes with largest measured soft\-fraction gaps. Prediction labels are ignored. A low predicted fraction is not biological absence; omitted classes remain explicit report warnings.

![Largest measured unique ordered\-pair residuals. Official reconciliation uses retained occurrence weights, not the displayed unweighted bar heights. Covariance remains separate descriptive evidence, not regulatory interaction.](figures/pair_discrepancies.png)
B — defined: Largest measured unique ordered\-pair residuals. Official reconciliation uses retained occurrence weights, not the displayed unweighted bar heights. Covariance remains separate descriptive evidence, not regulatory interaction.

![Measured WT\-relative deltas on operational target\-DE genes. A line is descriptive even when the official reliability gate rejects it; the displayed gate determines official severity interpretation. No gene residual is an allocation of log severity.](figures/response.png)
A — defined: Measured WT\-relative deltas on operational target\-DE genes. A line is descriptive even when the official reliability gate rejects it; the displayed gate determines official severity interpretation. No gene residual is an allocation of log severity.

![Independent centered native XY projections. Axes are arbitrary supplied frames, unregistered and unpaired; visual differences are not pointwise errors or laterality evidence.](figures/spatial_clouds.png)
D — defined: Independent centered native XY projections. Axes are arbitrary supplied frames, unregistered and unpaired; visual differences are not pointwise errors or laterality evidence.

![Separate illustrative draw of 10,000 ordered pairs per selected cloud; self\-pairs excluded. Normalization removes size and arbitrary rigid frames. This is sampled evidence, distinct from the official D2 draw, and blind to reflection. Shared histogram bins include the complete visual draw.](figures/shape_distances.png)
B — defined: Separate illustrative draw of 10,000 ordered pairs per selected cloud; self\-pairs excluded. Normalization removes size and arbitrary rigid frames. This is sampled evidence, distinct from the official D2 draw, and blind to reflection. Shared histogram bins include the complete visual draw.

![Direct size and count summaries. The official count scalar uses selected counts; original counts are descriptive metadata only. Neither identifies biological growth or proliferation.](figures/spatial_size_counts.png)
D — defined: Direct size and count summaries. The official count scalar uses selected counts; original counts are descriptive metadata only. Neither identifies biological growth or proliferation.

![Signed finite\-sample kernel terms in independent native XY frames. Both prediction and truth terms can be negative; they are not per\-cell causal errors or anatomical correspondence. No local\-priority threshold is adopted; effective upstream MMD seed is 0.](figures/local_neighborhood.png)
A — defined: Signed finite\-sample kernel terms in independent native XY frames. Both prediction and truth terms can be negative; they are not per\-cell causal errors or anatomical correspondence. No local\-priority threshold is adopted; effective upstream MMD seed is 0.


## Evidence taxonomy and limits

A: direct or algebraic selected-sample evidence; B: sampled evidence; C: artificial noncausal sensitivity; D: descriptive companion statistics. Exact/direct does not mean exact biology. Gene and composition calculations are exact for the bounded selected cells and supplied probe within stated numerical reconciliation tolerance; pair calculations are sampled. The source-provided gene panel is used as supplied and this report is not an official leaderboard or board-validation result. The selected subset is explicitly not evidence of full-file equivalence. Probe labels are limited to its domain, so a low predicted proportion is not biological absence. Variogram/covariance diagnostics do not establish gene regulation. Any max-dense input limit controls selected input conversion only; it is not a total process peak-memory guarantee.

## Warnings

- All scores and diagnostics apply to the selected supplied\-cell subset; they do not certify unscanned original rows.
- Biological states are assessed only through global gene errors and frozen\-probe proportions; no state\-conditional expression attribution is computed.
- The target\-informed matrix edit is noncausal sensitivity and may alter normalization.
- Pinned sliced Wasserstein and occupancy Dice have a PCA\-frame handedness defect: proper rotation or independent row permutation can change their scores, even for identical clouds. Treat these as provisional frame\-dependent diagnostics, not standalone failure priorities or laterality tests.
- Spatial scores and neighborhoods apply only to selected supplied cells; full\-file count is a separate descriptive diagnostic.
- D2 uses independent random pair draws and can have a nonzero identical\-cloud floor.
- No anatomical correspondence, regional causation, or calibrated laterality is assessed.
- Runtime FutureWarning: 'n\_jobs' has no effect since 1.8 and will be removed in 1.10. You provided 'n\_jobs=\-1', please leave it unspecified.
- Null values in this report mean undefined, not applicable, or non\-finite input values; see each diagnostic definition and limitation for its reason.
