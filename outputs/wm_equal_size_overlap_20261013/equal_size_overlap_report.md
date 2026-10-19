# Equal-sized region-set overlap across WM categories

For each category and model, rank all 489 regions and select exactly K. Compare the intersection across body, faces, places and tools. Selections use no significance filter. RMS ranks the magnitude of each population-average residual contrast curve. S=sum(T^2) ranks standardized departures from zero across 41 times. S is not a p-value, and differing null distributions can change its relationship to p-value rankings.

Three figures show four-category common counts, common count divided by K, and mean Jaccard similarity across the six category pairs. Displayed K ranges from 50 to 300; the MAT file contains every K from 1 to 489. All models share the same K and eligible ROI universe.

## RMS magnitude

| K | cHRF common count | cHRFderiv common count | sHRF common count | cHRF fraction | deriv fraction | sHRF fraction |
|---:|---:|---:|---:|---:|---:|---:|
| 50 | 12 | 12 | 13 | 0.240 | 0.240 | 0.260 |
| 75 | 25 | 22 | 25 | 0.333 | 0.293 | 0.333 |
| 100 | 38 | 31 | 37 | 0.380 | 0.310 | 0.370 |
| 150 | 64 | 65 | 67 | 0.427 | 0.433 | 0.447 |
| 200 | 115 | 105 | 107 | 0.575 | 0.525 | 0.535 |
| 250 | 167 | 149 | 148 | 0.668 | 0.596 | 0.592 |
| 300 | 213 | 212 | 186 | 0.710 | 0.707 | 0.620 |

| K | Mean pairwise Jaccard cHRF | cHRFderiv | sHRF |
|---:|---:|---:|---:|
| 50 | 0.361 | 0.342 | 0.348 |
| 75 | 0.418 | 0.401 | 0.421 |
| 100 | 0.461 | 0.419 | 0.451 |
| 150 | 0.517 | 0.524 | 0.517 |
| 200 | 0.624 | 0.598 | 0.575 |
| 250 | 0.701 | 0.672 | 0.620 |
| 300 | 0.730 | 0.744 | 0.668 |

## Standardized statistic S

| K | cHRF common count | cHRFderiv common count | sHRF common count | cHRF fraction | deriv fraction | sHRF fraction |
|---:|---:|---:|---:|---:|---:|---:|
| 50 | 9 | 6 | 12 | 0.180 | 0.120 | 0.240 |
| 75 | 25 | 19 | 20 | 0.333 | 0.253 | 0.267 |
| 100 | 40 | 33 | 33 | 0.400 | 0.330 | 0.330 |
| 150 | 82 | 75 | 72 | 0.547 | 0.500 | 0.480 |
| 200 | 128 | 114 | 114 | 0.640 | 0.570 | 0.570 |
| 250 | 165 | 149 | 155 | 0.660 | 0.596 | 0.620 |
| 300 | 209 | 206 | 193 | 0.697 | 0.687 | 0.643 |

| K | Mean pairwise Jaccard cHRF | cHRFderiv | sHRF |
|---:|---:|---:|---:|
| 50 | 0.317 | 0.319 | 0.340 |
| 75 | 0.422 | 0.390 | 0.369 |
| 100 | 0.484 | 0.450 | 0.450 |
| 150 | 0.601 | 0.576 | 0.529 |
| 200 | 0.651 | 0.626 | 0.603 |
| 250 | 0.684 | 0.662 | 0.641 |
| 300 | 0.725 | 0.724 | 0.673 |

## Interpretation limits

- Equal K controls set size, while preserving each model/category spatial ranking. It does not equalize sampling uncertainty or power.
- Uniform positive scaling of every ROI RMS in a category preserves its RMS ranking. A lower overlap can reflect relative spatial ordering or noisier ranks, beyond uniform magnitude reduction.
- Common RMS-ranked membership does not establish matching curve signs or shapes. RMS ignores sign, and S squares the standardized contrast.
- The comparison is descriptive. No overlap p-values, model-difference tests, subject bootstrap or spatial-null analysis were run.
- Exact score ties use ascending atlas index. There are 0 tied cutoffs at the reported K values across all rankings/categories/models.
- Permutation p-value ties and their Monte Carlo precision do not enter these rankings. Sampling uncertainty in observed means and standard errors remains.
- Region identities and rank score order were checked against mapping in the saved source results. Explicit set intersections agree with the cumulative-rank method.

Under independently uniform random K-region sets, expected common count is K^4/489^3 and expected common fraction is (K/489)^3. This is an illustrative size benchmark, not an anatomical significance test. At K=100 the expected count is 0.855. Spatial dependence and shared subjects invalidate treating it as the study null.

## Sources

- cHRF: `C:\Users\Ryan_\Documents\my_Git\fMRI\data\task_residual\wm_load_omnibus\AllFourCategory_LoadDiff\cHRF\job_36166327\temporal_omnibus_results.mat`
- cHRFderiv: `C:\Users\Ryan_\Documents\my_Git\fMRI\data\task_residual\wm_load_omnibus\AllFourCategory_LoadDiff\cHRFderiv\job_36166328\temporal_omnibus_results.mat`
- sHRF: `C:\Users\Ryan_\Documents\my_Git\fMRI\data\task_residual\wm_load_omnibus\AllFourCategory_LoadDiff\sHRF\job_36166326\temporal_omnibus_results.mat`
