# Regions significant in all four stimulus categories

The maps show the same atlas regions passing the 2-back minus 0-back residual omnibus test for body, faces, places and tools within each HRF model:

| Model | Regions |
|---|---:|
| cHRF | 220 |
| cHRFderiv | 231 |
| sHRF | 52 |

Files:

- `figures/cHRF_all_four_orthviews.png`
- `figures/cHRFderiv_all_four_orthviews.png`
- `figures/sHRF_all_four_orthviews.png`

Every map uses the same MNI coordinates [0 -20 20] mm. The overlays color deterministic canlab2018 atlas parcels. Colors denote significant-region membership; they do not indicate effect magnitude or direction. Three slices show only part of each set. Complete ROI lists and labels are in `category_overlap_report.md`.

Suggested caption:

> Colored regions show significant population-average memory-load residual contrasts in all four stimulus categories under the indicated HRF model. Shared significance indicates recurring load-dependent residual structure, without establishing identical directions or temporal profiles across categories.

These maps retain the original correction across 1,956 ROI-category hypotheses per model and the original Monte Carlo precision qualification. They add no voxelwise tests or new overlap significance tests. Atlas mapping uses the original extraction index convention and verifies the current label order against the saved summary.

To display interactively from the repository root, with CANlab and SPM on the MATLAB path:

```matlab
addpath('code/residual_contrast/analysis/wm_population_omnibus');
plot_wm_all_category_orthviews('cHRF');
plot_wm_all_category_orthviews('cHRFderiv');
plot_wm_all_category_orthviews('sHRF');
```

To export at different slice coordinates, supply a new filename:

```matlab
plot_wm_all_category_orthviews('cHRF', [0 -60 20], 'cHRF_all_four_posterior.png');
```
