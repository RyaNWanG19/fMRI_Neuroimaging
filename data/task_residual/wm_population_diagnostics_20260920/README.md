# Local diagnostics, 2026-09-20

All three Body_LoadDiff_vs_Face_LoadDiff models completed local diagnostics.
These runs used 31 sign flips and 31 bootstrap replicates ONLY to validate
the workflow. Their p-values and decision comparisons are not final inference.
The descriptive diagnostics use all retained subjects, not a subject subsample.

| Model | Median maximum absolute lag skewness across ROIs | Largest absolute lag skewness | Largest leave-one-out mean displacement (SE norm) |
|---|---:|---:|---:|
| cHRF | 0.518517 | 4.66538 | 1.65282 |
| cHRFderiv | 0.500867 | 4.73078 | 1.80536 |
| sHRF | 0.465358 | 3.70811 | 2.03596 |

Skewness is the centered third moment divided by the second moment to power
1.5. The maximum spans the 41 lags in each ROI. These are descriptive values;
they neither verify whole-curve symmetry nor automatically trigger exclusions.
Influence is the Euclidean norm across nonconstant lags of the change in the
mean after omitting one subject, divided by full-sample standard errors.

The owner confirmed distinct unrelated participants. Real IDs were unavailable;
row IDs describe each file's order and do not establish cross-model matching.
Symmetry was accepted as a modeling assumption, not marked as verified.

The full JHPCE helper is `submit_wm_population_full_jhpce.sh`. It uses separate
output directories and 99,999 replicates for each resampling method. It has not
been submitted remotely by this local diagnostic run. Bootstrap inference is
approximate and assumes iid retained curves with suitable moments; sign-flip
inference remains conditional on symmetry. Holm families are separate for
each model. The sensitivity analysis must not be used to pick the smaller p.
