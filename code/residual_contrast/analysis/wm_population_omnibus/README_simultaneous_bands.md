# Simultaneous bands for WM residual difference curves

This adds time-resolved inference alongside the existing temporal omnibus
test. The omnibus statistic `sum(T.^2)` and its Holm-adjusted p-values are
unchanged. Bands are optional and calibrated over the entire planned family,
before choosing any regions to display. They do not use the omnibus S statistic
or depend on omnibus significance.

For the four-category driver the family is **489 ROIs x 4 WM categories x
41 sampled lags = 80,196 cells**, jointly in both directions, within one HRF
model. For the interaction driver the family is all its ROIs x sampled lags.
There is no joint correction across HRF models or across separate driver runs.

## Calculation

Use the original subject curves with no extra smoothing. For each test r,
retain the same complete-case source subjects as the omnibus engine, and
form `e_i(r,t) = D_i(r,t) - mean_D(r,t)`. The default is a **centered
whole-subject bootstrap-t**. Sample N subject indices with replacement,
**shared across all ROIs, categories and times**. At each test, use only
sampled subjects eligible under that test's original mask; recompute its
bootstrap sample size and standard error. For each replicate calculate

```
Z_b(r,t) = mean_resampled_centered_curve(r,t) / SE_resampled(r,t)
Q_b      = max_{r,t} abs(Z_b(r,t))
SE(r,t)  = sample_sd(r,t) / sqrt(n_r)
band     = mean_D(r,t) +/- q * SE(r,t)
```

Here q is the upper empirical (1-alpha) quantile of the bootstrap maxima,
using order `ceil((B+1)*(1-alpha))` without interpolation. If this order exceeds
B, the bands stay unbounded. Defaults: alpha inherited from the omnibus
config (0.05), B=9,999, independent deterministic RNG seed, minSubjects=10.
The minimum sample guard is practical, not evidence that ten subjects suffice
for a large simultaneous family. Monte Carlo calibration itself is approximate;
repeat with more replicates when crossings are borderline.

Joint subject draws preserve temporal and regional dependence, including
overlapping subject masks. Resampled zero variance or fewer than two eligible
draws count as an infinite maximum; if these occur often enough, all bands
become unbounded. This is approximate bootstrap inference, not an exact
sign-flip confidence set. It estimates centered estimation errors, rather
than flipping the uncentered curves used in the zero-null omnibus test.

An optional `cfg.scb.method='gaussian-multiplier'` implementation uses one
shared N(0,1) subject weight, with
`Z_b=sum_i g_i*e_i/sqrt(sum_i e_i^2)` and fixed observed SE. Its conditional
covariance reproduces the sample correlation of the mean estimators. It is
retained for sensitivity comparison; the bootstrap-t version is the default
because the initial asymmetric-noise coverage pilot favored studentization.

Bands target approximate simultaneous coverage under independent subject
vectors, suitable finite moments, nondegenerate variances and an adequate
sample size for the chosen joint bootstrap approximation. They do not require null
sign symmetry. Running both analyses together still requires the existing
symmetry gate for the omnibus test. Family-related subjects need a justified
dependence-aware method; changing that gate does not implement one.

Zero-variance cells, cells with nonfinite uncertainty, and tests below
minSubjects receive `[-Inf, Inf]` with explicit status. They never produce
zero-width confidence intervals. Every planned cell remains in the export;
only eligible cells contribute to the maximum, while ineligible cells have
unbounded bands. Effective sample sizes and exclusions are retained in MAT.
Interpretation after missing-data exclusions requires a valid population
assumption. Audit source extraction Flags/missing runs: unavailable conditions
may have been initialized to zero, which is not automatically identifiable
from the Results array. This extension does not repair first-stage estimation
bias, source missing-data coding, or raw-scan autocorrelation modeling.

## Run with the four-category analysis

```matlab
addpath('code/residual_contrast/analysis/wm_population_omnibus');
cfg = struct('hrfModelName','cHRF','runSimultaneousBands',true);
cfg.scb = struct('nBootstrap',9999);
% Set these ONLY using documented evidence / accepted assumptions:
cfg.resampling.independentSubjectsVerified = true;
cfg.resampling.nullSymmetryAssumed = true; % for the omnibus, not the bands
cfg.resampling.verificationNote = 'Replace with actual independence evidence and symmetry assumption';
results = run_wm_load_omnibus(cfg);
```

Use `run_wm_population_omnibus(cfg)` instead for the existing category-load
interaction. `cfg.scb.seed`, `batchSize`, `testBatchSize`, and `minSubjects`
are optional. Batch sizes change memory use, not the resampling design.
`cfg.scb.method` defaults to `'subject-bootstrap-t'`.
`cfg.scb.nBootstrap` is separate from `cfg.nBootstrap`, which controls the
older centered-bootstrap omnibus sensitivity analysis.

For a non-inferential pipeline check, use `cfg.smokeTest=true`,
`cfg.nPerm=99`, and `cfg.scb.nBootstrap=99`, without declaring real-data
assumptions verified. Smoke-test band crossings must not be reported as
inferential findings. The new code does not automatically submit any jobs.

On JHPCE, add `RUN_SIMULTANEOUS_BANDS=1 SCB_NBOOTSTRAP=9999` to the existing
load-difference or interaction launcher, preserving the required subject
independence/symmetry environment settings. Use a fresh OUTPUT_ROOT/job
directory; existing result files are protected from overwriting.

## Outputs and interpretation

`temporal_omnibus_results.mat` now includes `results.simultaneousBands`:
mean, SE, lower/upper bounds, positive/negative indicators, cell eligibility
and status, bootstrap maxima and degeneracy count, critical value, configuration, correction
family and conditional inference status. `results.categoryBandCounts` in
the four-category driver counts regions with at least one supported lag.

`simultaneous_confidence_bands.csv` has one row per planned test and lag,
including effect, SE, bounds, eligibility, direction indicators and inference
status. The four-category export also includes atlas_roi_index and category.
These are jointly corrected band crossings, not pointwise p-values.

- Lower bound > 0: supported positive mean difference at that sampled lag.
- Upper bound < 0: supported negative mean difference at that sampled lag.
- Band contains zero: insufficient evidence of a signed effect at that lag;
  not evidence of equality or adequate residual model fit.
- A significant omnibus test can coexist with no band excluding zero,
  because the omnibus accumulates evidence over lags.

Coverage is at the **41 sampled lags only**. Connected graphical lines and
shading do not assert coverage between samples or establish exact biological
onset, offset or duration. Bands are approximate frequentist intervals, not
a posterior probability that the observed curve lies within them.

Plot any 1-12 tests from the result without recalibrating the family:

```matlab
% ROI 1 in body/faces/places/tools:
plot_simultaneous_confidence_bands(results,[1 490 979 1468], ...
    fullfile(results.config.outputDir,'roi_001_bands.png'));
```

Generic single-contrast omnibus curve plots use simultaneous bands when
available; the full four-category driver suppresses those generic plots.
The explicit plotting helper supports the four-category mapping.

## Verification

```matlab
addpath('tests');
test_simultaneous_confidence_bands;
test_wm_load_omnibus;
test_wm_population_omnibus;
test_omnibus_sensitivity;
report = validate_simultaneous_band_coverage('',200,999);
```

Checks cover a direct numerical reference, shared multipliers, exact duplicate
curves, fixed missingness masks, batching reproducibility, translation/sign/
scale behavior, degenerate cells, inference gates, CSV/MAT mapping, overwrite
protection and unchanged omnibus decisions. The coverage experiment uses
Gaussian and asymmetric temporally correlated curves with shared variation
across four tests and 41 lags, n=200. It is a finite-sample check in a smaller
family, not proof of 95% coverage for the full 80,196-cell family or actual
HCP data. Repeat calibration experiments with representative subject counts,
covariance and missingness before making final study claims.

Methods background:
[Degras (2011)](https://arxiv.org/abs/0908.1980),
[Telschow & Schwartzman](https://pmc.ncbi.nlm.nih.gov/articles/PMC9268949/).
