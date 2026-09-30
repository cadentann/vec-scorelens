# VEC ScoreLens

ScoreLens helps you inspect where a supplied prediction differs from a permissible local target: named gene discrepancies, probe-class proportions, sampled gene pairs, and task-specific response or spatial evidence. It supports Task 1 (T1), Task 2 (T2), and Task 3 (T3), producing private local reports alongside the pinned `veckit` scorer's unchanged scalar dictionary for a deterministic bounded subset. It does not download challenge data, contact a service, submit to a leaderboard, or obtain hidden targets.

Use only locally supplied targets that you are permitted to analyze. A pseudo-target made from released training observations is internal cross-validation, not hidden-stage performance. Keep reports private when inputs are private.

## Install

Clone the public repository, then install from its root with Python 3.11 or later:

```sh
git clone https://github.com/cadentann/vec-scorelens.git
cd vec-scorelens
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install '.[test,notebook]'
```

`matplotlib` is a normal required dependency: reports use locally generated figures and do not load remote assets. The package also pins `veckit` to commit `46d41e63f42a9aab815db20b742feeccd249cb17`; ScoreLens verifies the installed version, commit provenance, and reviewed source hashes before scoring.

## Five-minute synthetic walkthroughs

Each command writes only to a new directory. The examples are expression-like synthetic fixtures, not Challenge data or biological simulations.

```sh
# T1: reference-relative expression diagnostics
vec-scorelens-demo --task T1 --out local_outputs/t1_inputs --case direction
vec-scorelens --task T1 \
  --prediction local_outputs/t1_inputs/prediction.h5ad \
  --target local_outputs/t1_inputs/target.h5ad \
  --reference local_outputs/t1_inputs/reference.h5ad \
  --out local_outputs/t1_report

# T2: spatial diagnostics with the required supplied reference
vec-scorelens-demo --task T2 --out local_outputs/t2_inputs --case mixed
vec-scorelens --task T2 \
  --prediction local_outputs/t2_inputs/prediction.h5ad \
  --target local_outputs/t2_inputs/target.h5ad \
  --reference local_outputs/t2_inputs/reference.h5ad \
  --out local_outputs/t2_report

# T3: matched-WT response diagnostics
vec-scorelens-demo --task T3 --out local_outputs/t3_inputs --case mixed
vec-scorelens --task T3 \
  --prediction local_outputs/t3_inputs/prediction.h5ad \
  --target local_outputs/t3_inputs/target.h5ad \
  --wt local_outputs/t3_inputs/reference.h5ad \
  --out local_outputs/t3_report
```

The standard output identifies the report directory. Inspect its `summary.md`, `report.html`, `diagnostics.json`, CSV tables, `provenance.json`, and local `figures/` directory. These are generated outputs from the command, not an evaluation certificate, and a report directory is never overwritten.

Open the generated [T3 synthetic example report](examples/generated/T3/report.html) for the matched-WT response and spatial layout. The corresponding [T1](examples/generated/T1/report.html) and [T2](examples/generated/T2/report.html) synthetic reports are also included. These static artifacts contain generated synthetic values only and omit JSON/CSV ledgers; they are report-layout examples, not biological or Challenge results. Running the CLI creates the full private report.

## Why ScoreLens exists

`veckit` is the authoritative pinned scorer; ScoreLens does not replace it, introduce a new metric, or calibrate its outputs. VecBench and the Local Scoring Wrapper already accept users' predictions and provide aggregate calibration, noise, or reliability workflows. Intuition Lab teaches controlled task failures, pins the same scorer, and exposes scoring helpers that accept supplied prediction paths. The Submission Viewer already diagnoses spatial coverage and coordinate-frame problems in users' predictions. ScoreLens's narrower addition is a companion report that exports named gene and fitted-probe class terms, occurrence-weighted sampled pair evidence, T3 WT-relative response residuals and severity gates, and constrained spatial metric evidence alongside the scorer's native values.

This distinction is bounded. Named rows are not universal all-gene-pair localization, a new calibration layer, or a general spatial viewer. T2/T3 evidence must be read as selected-sample score-adjacent diagnostics rather than registered biology. A dated, public-facing comparison with its source links and limits is in [NOVELTY.md](NOVELTY.md).

## Run notebooks and checks

From the installed repository root, launch the three clean local notebooks with:

```sh
python -m jupyter lab notebooks/scorelens_synthetic.ipynb
python -m jupyter lab notebooks/scorelens_t2.ipynb
python -m jupyter lab notebooks/scorelens_t3.ipynb
python -m pytest -q
```

For noninteractive execution, use `python -m jupyter execute notebooks/scorelens_t2.ipynb` (and the equivalent T1/T3 paths). The notebooks create temporary local paths and call the installed commands; they do not fetch data or use hidden filesystem locations.

## Inputs and scope

To inspect your own T1 files, replace the three input paths below with files you are permitted to analyze and choose a new report directory:

```sh
vec-scorelens --task T1 \
  --prediction /path/to/prediction.h5ad \
  --target /path/to/local_pseudo_target.h5ad \
  --reference /path/to/reference.h5ad \
  --out local_outputs/my_t1_report
```

For T2, change `--task T1` to `--task T2` and supply spatial files. For T3, use `--task T3` and replace `--reference` with `--wt` pointing to the matched wild-type file. A hidden organizer target is neither required nor obtainable through ScoreLens. On a rerun, use a different output directory; existing reports and inputs are never overwritten.

An interrupted process may leave a hidden `.<report-name>.stage-*` folder beside the intended output, containing partial private results. It is not a completed report. Inspect it and remove it only after confirming that no run is active and its contents are no longer needed. Ordinary Python failures clean up their staging folder; abrupt termination cannot guarantee cleanup.

T1 and T2 require `--reference`; T3 requires `--wt`. The roles are validated explicitly: ScoreLens never substitutes a contrast file. T2/T3 inputs must include finite three-dimensional coordinates in `obsm['spatial_3D']`; first three columns are used and remain row-aligned with the selected expression rows. T1 has no coordinate requirement.

Supported selected input matrices are finite, nonnegative dense, CSR, or CSC AnnData `.h5ad` data with the same unique ordered genes. Target `obs['celltype']` labels train the frozen probe; prediction labels are ignored. Inputs need at least 32 cells and two genes. These are companion input policies, not the Challenge submission schema.

By default ScoreLens samples at most 256 rows from each file; `--max-cells` permits 32–1024. All supplied genes are retained in T1 up to 40,000. T2/T3 reject files with more than 2,000 genes and retain every supplied gene within that limit. This companion resource policy for spatial neighborhood calculations is not a Challenge rule. Reports show original file counts and selected counts separately. Official scalars describe selected arrays only; original counts are descriptive metadata, never a substituted full-file score.

The pair ledger is a bounded random sample of ordered pairs, not an exhaustive all-gene-pair map; ScoreLens does not promise all-pair localization across a 32,285- or 40,000-gene T1 panel.

## What each task adds

All tasks retain per-gene pseudobulk/DE ledgers, frozen-probe composition terms, sampled variogram evidence, descriptive covariance, and the target-aware MMD sensitivity. T2 adds selected-cloud form, count/scale, distance-distribution, and local spatial-expression evidence. T3 adds those spatial views plus matched-WT response deltas and guarded severity-fit evidence. T3 changes are mean transformed-expression deltas (`prediction − WT` and `target − WT`); they are not raw-count log fold changes.

Priority findings cover only the published selected-sample conventions: literal response reversal, mean residual, eligible T3 response magnitude, frozen-probe composition, variance ratio, selected spatial size, and original-count metadata. The report explicitly leaves raw official-metric ranking, SW/Dice shape, sampled covariance, local-neighborhood organization, laterality, biological health, and sensitivity-based repair advice outside priority evaluation. A result within listed thresholds is never an all-clear for those unsupported domains or for full files.

See [METHODS.md](METHODS.md) for definitions, [LIMITATIONS.md](LIMITATIONS.md) and [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md) for boundaries, [TESTING.md](TESTING.md) for commands and acceptance status, and [ATTRIBUTION.md](ATTRIBUTION.md) for upstream and Challenge notices.

## License and credit

Copyright 2026 Caden Tan. Original ScoreLens source, tests, documentation, notebooks and synthetic examples are licensed under [Apache-2.0](LICENSE). Third-party materials retain their own licenses; the upstream `veckit` MIT notice is preserved in [NOTICE.md](NOTICE.md). This license does not cover Challenge materials or datasets, official panels, dependencies, external code, user inputs, or reports derived from data whose rights belong elsewhere.
