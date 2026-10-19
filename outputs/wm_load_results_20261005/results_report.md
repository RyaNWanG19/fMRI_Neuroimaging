# WM load-difference residual results

Each test assesses whether the population-average 2-back minus 0-back residual curve departs from zero somewhere across 41 sampled lags (0-28.8 s). The four categories are body, faces, places, and tools. Full runs used 99,999 subject sign flips, alpha 0.05, and Holm correction across 1,956 ROI-category tests separately within each HRF model.

## Significant regions

| Category | cHRF | cHRFderiv | sHRF | Shared by all three | At least one model | None |
|---|---:|---:|---:|---:|---:|---:|
| body | 328 (67.1%) | 339 (69.3%) | 132 (27.0%) | 126 | 362 | 127 |
| faces | 323 (66.1%) | 298 (60.9%) | 132 (27.0%) | 126 | 345 | 144 |
| places | 345 (70.6%) | 349 (71.4%) | 161 (32.9%) | 156 | 374 | 115 |
| tools | 305 (62.4%) | 340 (69.5%) | 112 (22.9%) | 107 | 356 | 133 |

The saved decisions show systematic load-dependent residual contrasts in many regions in every stimulus category. sHRF has fewer significant regions than either canonical model in all four categories. This is a descriptive comparison: counts alone do not establish that sHRF fits better or that the models differ statistically.

51 regions are significant in every category under every model. Per-category overlaps describe agreement between decisions, not independent replication.

## Exclusive model overlap

| Category | All three | cHRF + deriv only | cHRF + sHRF only | deriv + sHRF only | cHRF only | deriv only | sHRF only | None |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| body | 126 | 181 | 1 | 3 | 20 | 29 | 2 | 127 |
| faces | 126 | 153 | 2 | 1 | 42 | 18 | 3 | 144 |
| places | 156 | 166 | 0 | 3 | 23 | 24 | 2 | 115 |
| tools | 107 | 185 | 2 | 0 | 11 | 48 | 3 | 133 |

## Monte Carlo resolution and robustness

The minimum raw p-value is 0.00001. At the minimum, 1,956-test Holm adjustment yields 0.01956. With zero exceedances, the code estimates a raw-p Monte Carlo SE near 0.00001. Its upper two-SE input is then near 0.00003, which can produce a Holm upper sensitivity value near 0.05868. Because Holm uses cumulative maxima, this early sorted-p bound can propagate to many later hypotheses. The broad flags therefore do not mean each reported p-value is individually close to 0.05. This is a heuristic two-SE sensitivity screen, not a simultaneous confidence interval. It flags possible finite-permutation sensitivity, rather than demonstrating that all flagged decisions would change.

| Category | Flagged significant cHRF | Flagged significant deriv | Flagged significant sHRF |
|---|---:|---:|---:|
| body | 328 | 339 | 132 |
| faces | 323 | 298 | 132 |
| places | 345 | 349 | 161 |
| tools | 305 | 340 | 112 |

These runs do not contain centered-bootstrap sensitivity results. The earlier interaction bootstrap findings cannot be transferred to these tests. Before treating borderline classifications as stable, increase the permutation count or use a statistically justified Monte Carlo precision assessment. No additional inference was run for this report.

## Residual magnitude

Median RMS of each ROI population-average curve, over all 489 ROIs (native residual units):

| Category | cHRF | cHRFderiv | sHRF |
|---|---:|---:|---:|
| body | 0.0128842 | 0.0114313 | 0.00406276 |
| faces | 0.0119448 | 0.0096904 | 0.00384705 |
| places | 0.013486 | 0.0120045 | 0.00405174 |
| tools | 0.0110587 | 0.0110114 | 0.00394063 |

RMS describes curve magnitude and uses no uncertainty weighting. Significance instead depends on the time-specific mean-to-SE ratios and the permutation distribution. Magnitude distributions are descriptive and have no between-model p-value.

## Example trajectories

Each category has mean and T-statistic figures for five patterns: all three models, each model only, and none. Examples are the first atlas index in each pattern, chosen after testing. Empty patterns remain explicitly empty. The example curves are not representative estimates for their groups. Each category figure uses common y limits across its panels.

### body

| Pattern | ROI | Label | cHRF adjusted p | deriv adjusted p | sHRF adjusted p |
|---|---:|---|---:|---:|---:|
| All three | 1 | Ctx_V1_L | 0.01956 | 0.01956 | 0.01956 |
| cHRF only | 83 | Ctx_7AL_L | 0.02828 | 0.51092 | 1 |
| deriv only | 60 | Ctx_7m_R | 0.0524 | 0.01956 | 1 |
| sHRF only | 72 | Ctx_5m_R | 1 | 0.7328 | 0.01956 |
| None | 4 | Ctx_MST_R | 1 | 1 | 1 |

### faces

| Pattern | ROI | Label | cHRF adjusted p | deriv adjusted p | sHRF adjusted p |
|---|---:|---|---:|---:|---:|
| All three | 1 | Ctx_V1_L | 0.01956 | 0.01956 | 0.01956 |
| cHRF only | 19 | Ctx_FEF_L | 0.01956 | 0.12518 | 0.73101 |
| deriv only | 28 | Ctx_RSC_R | 0.67728 | 0.03305 | 1 |
| sHRF only | 68 | Ctx_d23ab_R | 0.24508 | 0.22428 | 0.02972 |
| None | 27 | Ctx_RSC_L | 1 | 0.13018 | 1 |

### places

| Pattern | ROI | Label | cHRF adjusted p | deriv adjusted p | sHRF adjusted p |
|---|---:|---|---:|---:|---:|
| All three | 1 | Ctx_V1_L | 0.01956 | 0.01956 | 0.01956 |
| cHRF only | 53 | Ctx_PCV_L | 0.01956 | 1 | 0.51865 |
| deriv only | 71 | Ctx_5m_L | 0.0976 | 0.0448 | 1 |
| sHRF only | 199 | Ctx_OP4_L | 0.44892 | 1 | 0.01956 |
| None | 47 | Ctx_A1_L | 1 | 1 | 1 |

### tools

| Pattern | ROI | Label | cHRF adjusted p | deriv adjusted p | sHRF adjusted p |
|---|---:|---|---:|---:|---:|
| All three | 1 | Ctx_V1_L | 0.01956 | 0.01956 | 0.01956 |
| cHRF only | 79 | Ctx_24dd_L | 0.02828 | 0.14985 | 1 |
| deriv only | 47 | Ctx_A1_L | 1 | 0.01956 | 1 |
| sHRF only | 17 | Ctx_3b_L | 0.40872 | 1 | 0.02972 |
| None | 15 | Ctx_4_L | 0.0756 | 0.39078 | 0.36064 |

## Anatomy and interpretation boundaries

- Orthviews use deterministic canlab2018 parcel boundaries and the project atlas-index convention. Historical atlas-version identity cannot be independently established from ordinal IDs alone.
- All anatomical maps use the same MNI coordinates [0 -20 20] mm. Three slices do not display every significant region. Color denotes ROI membership, not effect sign, magnitude, or voxelwise significance.
- A significant omnibus test provides evidence of a nonzero population-average load contrast somewhere in the sampled window. It does not locate significant lags.
- Shared mismatch under 2-back and 0-back can still cancel. These are not absolute condition-specific residual tests.
- The sign-flip inference is conditional on the documented independence and joint null-symmetry assumptions. Numerical agreement does not verify those assumptions.
- Correction covers four categories within each model, not all models or the earlier interaction analysis together. Different correction families prevent direct interpretation of count changes as effect changes.

## Checks and provenance

All 5,868 observed ROI-category/model means, SDs, T curves, S statistics and subject masks were recomputed from the local source residual arrays and matched within numerical tolerance. Saved adjusted p-values and labels match the 1,956-test Holm calculation. This check does not independently regenerate the 99,999 permutations or audit upstream HRF fitting.

- cHRF: N 402-410, source `C:\Users\Ryan_\Documents\my_Git\fMRI\data\task_residual\wm_load_omnibus\AllFourCategory_LoadDiff\cHRF\job_36166327\temporal_omnibus_results.mat`
- cHRFderiv: N 402-410, source `C:\Users\Ryan_\Documents\my_Git\fMRI\data\task_residual\wm_load_omnibus\AllFourCategory_LoadDiff\cHRFderiv\job_36166328\temporal_omnibus_results.mat`
- sHRF: N 402-410, source `C:\Users\Ryan_\Documents\my_Git\fMRI\data\task_residual\wm_load_omnibus\AllFourCategory_LoadDiff\sHRF\job_36166326\temporal_omnibus_results.mat`

## Slide-ready interpretation

Across all four stimulus categories, the population-average 2-back minus 0-back residual contrast remained significant in many regions. The sHRF model yielded fewer significant regions than cHRF or cHRFderiv, while a subset of regions passed under all three models. These results identify load-dependent residual structure. They do not establish superior fit from counts alone or exclude mismatch shared across loads. Finite-permutation uncertainty warrants caution when interpreting exact significance boundaries.
