# Real-data analysis for PAC-calibrated regression tolerance intervals

This folder contains the reproducible first-pass analysis for:

1. CPS 2012 hourly wages
2. Superconductivity critical temperature
3. MEPS 2016 total medical expenditure

MEPS health-care utilization, the outcome used by the public CQR preprocessing
code, is included as a sensitivity analysis.

## Main settings

- Target content: `C = 0.90`
- Outer confidence: `1 - alpha = 0.95`
- Split: approximately one third training, calibration, and evaluation
- Repetitions: 20
- Learner: gradient boosting with 300 trees
- Methods: Parametric-TI, SR-TI, ASR-TI, CQR-TI
- ASR tail-scale quantile: `tau = 0.10`
- Calibration: empirical `C + sqrt(log(2/alpha)/(2*n_cal))` score quantile

The parametric comparator uses the same boosted mean estimator, a global honest
residual standard deviation, and a normal/chi-square tolerance factor. It is an
explicit reproducible benchmark, not a reconstruction of unavailable manuscript
source code.

## Files

- `R/00_download_data.R`: downloads/bootstrap public data.
- `R/01_prepare_data.R`: constructs model matrices and transformed responses.
- `R/tolerance_methods.R`: implements the four intervals and diagnostics.
- `R/guo2024_pointwise_ti.R`: exact two-sided pointwise Guo--Young
  nonparametric regression TI using equations (11) and (14) with a
  GCV-selected one-dimensional cubic smoothing spline.
- `R/02_run_experiments.R`: runs repeated train/calibration/evaluation splits.
- `R/03_summarize_results.R`: produces tables, plots, and LaTeX.
- `tables/`: manuscript-ready numerical summaries.
- `figures/`: coverage, width, conditional diagnostic, and score-pivotality plots.
- `real_data_section_compact.tex`: drop-in replacement matching the original
  section's concise dataset--table--takeaway structure.
- `real_data_section_revised.tex`: recommended expanded section with
  split-success and conditional score-pivotality diagnostics.
- `report_ko.md`: Korean interpretation and manuscript recommendations.

## Re-run

From the Codex workspace root:

```bash
Rscript outputs/real_data_analysis/R/00_download_data.R --accept-meps-terms
Rscript outputs/real_data_analysis/R/01_prepare_data.R
Rscript outputs/real_data_analysis/R/02_run_experiments.R --reps=20 --cores=4 --trees=300
Rscript outputs/real_data_analysis/R/03_summarize_results.R
```

Before using `--accept-meps-terms`, review the
[AHRQ MEPS HC-192 data-use terms](https://meps.ahrq.gov/data_stats/download_data/pufs/h192/h192doc.shtml#DataA).

## Data sources

- CPS: R package `hdm::cps2012` and the
  [Distributional Conformal Prediction replication repository](https://github.com/kwuthrich/Replication_DCP).
- Superconductivity:
  [UCI Superconductivity Data](https://archive.ics.uci.edu/dataset/464/superconductivty%2Bdata).
- MEPS:
  [AHRQ HC-192](https://meps.ahrq.gov/data_stats/download_data/pufs/h192/h192doc.shtml)
  and the [CQR MEPS preprocessing repository](https://github.com/yromano/cqr/tree/master/get_meps_data).

## Interpretation boundary

The held-out results are empirical diagnostics. They do not directly verify the
unknown conditional PAC probability. MEPS has a point mass at zero and a complex
survey design; the current results target the observed sample distribution and
do not incorporate survey weights. Superconductivity may contain dependent
material families. These limitations should be stated in the manuscript.
