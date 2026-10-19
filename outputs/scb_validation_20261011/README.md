# Simultaneous band validation, 2026-10-11

Synthetic validation only; no HCP residual data were analyzed here.
The current implementation and its use are described in
[the analysis README](../../code/residual_contrast/analysis/wm_population_omnibus/README_simultaneous_bands.md).

Each method was evaluated on the same 200 independently generated datasets
per scenario, n=200 subjects, four tests x 41 lags (164 cells), and 999
bootstrap replicates. Temporal innovations have AR(1) correlation rho=0.85,
and shared subject variation correlates tests. The asymmetric scenario uses
centered, standardized lognormal innovations. The population mean is a known
nonzero sine curve. Coverage means ALL 164 true means lie within their bands.

| Method | Gaussian curve family | Asymmetric curve family |
|---|---:|---:|
| Gaussian multiplier, fixed SE | 94.5% (MC SE 1.61 percentage points) | 91.5% (MC SE 1.97 points) |
| Whole-subject bootstrap-t (default) | 94.5% (MC SE 1.61 points) | 95.5% (MC SE 1.47 points) |
| Pointwise +/-1.96 SE, evaluated jointly | 6.0% | 9.0% |

The initial Gaussian-multiplier report and example are in this directory.
The bootstrap-t report and example are in `bootstrap_t/`. The CSV files
pre-date addition of the explicit method column to the reusable validation
runner; directory names and this note identify the methods used.

The finite pilot supports using bootstrap-t as the default, but 200
experiments cannot precisely determine nominal coverage. These results do
not establish coverage of the actual 80,196-cell WM family, arbitrary sample
sizes, related subjects, missingness or other noise distributions.

Reproduce using the current runner (choose new output directories):

```matlab
addpath('tests');
validate_simultaneous_band_coverage('',200,999,'gaussian-multiplier');
validate_simultaneous_band_coverage('',200,999,'subject-bootstrap-t');
```

`synthetic_band_example.png` in each directory is the first Gaussian dataset,
with corrected positive/negative crossings marked at observed lags. It is an
illustration, not a real-data finding.
