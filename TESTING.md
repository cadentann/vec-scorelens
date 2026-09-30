# Testing and acceptance status

Run from a fresh environment installed with `python -m pip install '.[test,notebook]'`:

```sh
python -m pytest -q
python -m jupyter execute notebooks/scorelens_synthetic.ipynb
python -m jupyter execute notebooks/scorelens_t2.ipynb
python -m jupyter execute notebooks/scorelens_t3.ipynb
python benchmarks/memory_probe.py
```

The notebooks use temporary directories and local installed APIs only. They do not use Colab, hidden paths, or network fetches.

## Coverage and interpretation

The suite covers legacy T1 compatibility; task-specific T1/T2/T3 public-scorer parity; reference/WT and coordinate validation; response guard branches; selected/original count disclosure; source-pin verification; selective HDF5 reads; synthetic planted failures and healthy controls; priority-rule edge cases; report consistency; and all three notebooks. Notebook execution inside pytest overlaps the explicit notebook commands above: do not add those counts as independent experiments.

Development checks used macOS 26.5.1 arm64 and Python 3.11.0, with NumPy 2.4.6, SciPy 1.17.1, scikit-learn 1.9.1, AnnData 0.12.19 and h5py 3.16.0. Installation allows version ranges, so record actual versions with `python -m pip freeze` and check them with `python -m pip check`. Other operating systems, all combinations of allowed versions, and Colab have not been assessed.

A prior full-width synthetic T1 workflow (256 selected rows, 32,285 genes) measured approximately 1.32 GB peak process RSS. The final-verification baseline at 128 selected rows measured 1.126 GB in one run and about 0.838 GB on a repeat. The default benchmark uses a 2 GiB full-workflow RSS test budget to allow this observed variability; it is not a hard whole-process memory guarantee for other data or settings. The benchmark separately rejects unselected sparse-matrix densification and verifies source SHA-256, size, and mtime before and after scoring. The loader's `--max-dense-mib` limit bounds an input-array estimate, not downstream sklearn/scorer/report memory. Large unused sidecars were checked with guarded HDF5 reads; a large logical fixture is not a benchmark of downloading or reading a large real dataset.

On restricted hosts, point caches to new writable local folders before running notebooks or reports:

```sh
export MPLCONFIGDIR="$PWD/local_outputs/matplotlib_cache"
export IPYTHONDIR="$PWD/local_outputs/ipython_cache"
```

Upstream scikit-learn may warn that its `n_jobs` option is deprecated; this does not change the pinned scoring source. Keep warning logs and record the commit and environment for your own run. Synthetic acceptance does not establish a Challenge score, biological validation, or a calibrated false-alarm guarantee.
