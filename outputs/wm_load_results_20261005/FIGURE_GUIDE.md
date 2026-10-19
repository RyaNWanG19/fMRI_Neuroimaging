# Figure guide

Start with `results_report.md` for counts, interpretation, methods, and uncertainty caveats. All figures use the original saved significance decisions, with Holm correction across 489 x 4 hypotheses separately within each HRF model. They are PNG files under `figures/`.

## Suggested short results sequence

1. **Number of significant regions** (`significant_region_counts.png`). Each bar counts regions with evidence of a population-average 2-back minus 0-back residual contrast somewhere across the 41 sampled times. sHRF has fewer discoveries in each category, but discovery counts alone do not establish better fit.
2. **Overlap across models** (`model_overlap.png`). Each category has eight mutually exclusive patterns. The figure includes the none-significant group. “deriv” means cHRFderiv. Shared significance is descriptive agreement between models fitted to the same data.
3. **Mean curves and standardized evidence** (`body_examples_mean.png`, `body_examples_T.png`, and corresponding faces/places/tools files). Pick one category for a short presentation. Mean curves show direction and amplitude; T curves show the mean relative to its standard error. Panel labels refer to whole-curve significance, not individual lags. Each example is the first atlas index in its significance group, selected after testing. Examples are not group averages or representative samples.
4. **Anatomical distribution** (`body_shared_orthviews.png` and corresponding faces/places/tools files). Green regions pass in all three models for that category. Individual-model maps are also provided as `<category>_<model>_orthviews.png`. Every image uses MNI [0 -20 20] mm. The three slices do not show every region; see `shared_region_lists.md` for complete per-category shared lists.
5. **Residual magnitude** (`residual_magnitude_distributions.png`, optional). These are cumulative distributions across all 489 regions. At an x-axis value, the y-axis gives the fraction with RMS mean contrast at or below that magnitude. A curve farther left generally indicates smaller RMS contrasts. Axis ranges vary across categories; within each category the model curves share axes. This display is descriptive and does not test model differences.

## Essential interpretation text

Systematic load-dependent residual contrasts remain in many regions across all four stimulus categories. The number and magnitude of these contrasts vary across HRF models. Significance depends on contrast magnitude relative to uncertainty, rather than mean amplitude alone. These tests can miss residual mismatch shared equally by 2-back and 0-back.

## Robustness qualification

These full runs contain 99,999 sign flips and no centered-bootstrap analysis. Every significant result is flagged by the existing two-standard-error Monte Carlo screen, largely because its upper bound at the minimum p-value propagates through the 1,956-test Holm calculation. The flags are a heuristic sensitivity screen, not simultaneous confidence intervals, and do not establish that every decision is unstable. Treat precise threshold classifications cautiously pending a higher-precision assessment.

## Reproduce

With the source residuals, downloaded results, CANlab and SPM available, run from the repository root:

```matlab
addpath('code/residual_contrast/analysis/wm_population_omnibus');
summary = summarize_wm_load_results( ...
    fullfile(pwd,'outputs','wm_load_results_new'), true);
```

Use a new output directory. The function does not overwrite previous analysis folders. `false` as the second argument skips anatomical rendering. `summary_data.mat` contains numeric summaries, example indices and p-values, atlas labels, per-ROI significance arrays and RMS magnitudes, source paths, and observed-statistic audit results.
