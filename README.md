This repository contains the analysis pipeline and code used for the thesis investigating systematic, task-dependent information within fMRI residuals. Using datasets from the Human Connectome Project (HCP), this project tests whether residuals contain structured information that can reliably predict task states.

## Paired-Contrast Tensor Decomposition (Option B)

The separate [paired tensor module](code/residual_contrast/analysis/paired_tensor/README.md)
fits signed CP decompositions to region x time x subject contrast tensors.
It includes multiple starts, a rank sweep, an uncentered SVD benchmark,
split-half component reproducibility, explicit subject-exclusion auditing,
and CSV/MAT/figure outputs. It is descriptive and does not replace population
inference. The default uses the WM body-versus-face load interaction.

```matlab
addpath('code/residual_contrast/analysis/paired_tensor');
results = run_paired_tensor(struct('hrfModelName','cHRF'));
```

## WM Population Temporal-Omnibus Test

`code/residual_contrast/analysis/wm_population_omnibus` contains the ROI-level
population analysis. For each parcel it tests the complete mean paired-difference
curve with `S = sum(T(t).^2)`, calibrates that statistic by subject-level sign
flipping, and applies Holm correction across all planned parcels for one contrast
and one HRF model. It operates directly on subject differences and is separate
from classifier importance or lag-specific inference.

The HRF-model files contain arrays in `[ROI x time x subject x condition]` order.
The driver verifies the repository's WM slots (0-back categories in 1--4 and
2-back categories in 5--8), constructs a category load-difference interaction,
and records the permutation to `[subject x ROI x time]`.

The files do not contain real HCP IDs or family/dependence metadata. A full run
requires documented confirmation that the rows are independent. Real IDs are
optional: absent IDs, row_001 etc. track observations within each model file.
These labels do not establish matching subjects across models.

```matlab
addpath('code/residual_contrast/analysis/wm_population_omnibus');
cfg = struct();
cfg.hrfModelName = 'cHRF';
cfg.contrastName = 'Body_LoadDiff_vs_Face_LoadDiff';
% cfg.subjectIDs = subjectIDs; % optional aligned real IDs
cfg.resampling.independentSubjectsVerified = true;
cfg.resampling.nullSymmetryAssumed = true;
cfg.resampling.verificationNote = 'Owner confirms distinct unrelated participants; null symmetry assumed';
cfg.runDiagnostics = true;
cfg.nBootstrap = 99999;
results = run_wm_population_omnibus(cfg); % default B = 99,999
```

For a non-inferential smoke run, set `cfg.smokeTest = true` (default B = 999).
The minimum attainable raw permutation p-value is `1/(B+1)`; smoke-run p-values
must not be reported as final inference.

On JHPCE, submit
`code/residual_contrast/analysis/wm_population_omnibus/run_wm_population_omnibus_jhpce.sbatch`.
It reads the HRF files from
`/users/rwang/fMRI_Neuroimaging/data/task_residual` by default. Environment
variables can override the repository, data, output, model, contrast, and
permutation settings. The full-run helper below uses the actual data files in
`/users/rwang`, records the owner's independence confirmation and accepted
symmetry assumption, and submits all three models:

```bash
bash code/residual_contrast/analysis/wm_population_omnibus/submit_wm_population_full_jhpce.sh
```

Each job runs 99,999 sign flips and 99,999 centered whole-curve bootstrap
replicates per ROI. The latter is an approximate sensitivity analysis allowing
asymmetry, with iid retained subjects and suitable finite moments required.
It is not an exact test or a substitute for the primary sign-flip result.
Both corrections include all 489 planned ROIs, separately per model.
`omnibus_sensitivity.csv` reports centered skewness, extreme standardized
residuals, leave-one-out mean influence, bootstrap p-values, and disagreements.
Diagnostics do not automatically exclude subjects or prove symmetry. Bootstrap
results for nonzero constant cells are marked unreliable and assigned p=1.
Methods background: [bootstrap principle and null calibration](https://www.stat.cmu.edu/~cshalizi/uADA/19/lectures/ch06.pdf).
The configuration records assumptions; `inferenceValid` means eligible under
those assumptions, not that their truth was established by the program.

## WM ROI Influence Elastic-Net Module

The separate module in `code/residual_contrast/analysis/wm_roi_influence` identifies stable, high-confidence influential candidate parcels for WM residual interaction contrasts without modifying the existing linear-SVM pipeline.

The driver `run_roi_influence_elasticnet.m` uses subject-level paired differences:

`D = (2-back category A - 0-back category A) - (2-back category B - 0-back category B)`

and fits nested elastic-net logistic models with grouped subject folds using the symmetric paired representation `[D/2; -D/2]`. Main outputs are kept distinct:

- Haufe pattern: task-linked residual signal map.
- Elastic-net selection frequency: stability of sparse multivariable selection.
- Held-out model reliance: dependence of the fitted decoder on each parcel under held-out paired half-swaps.

The population-level sign-flip maxT regional test is separate from the decoder and uses only unstandardized subject-level paired differences. The `high_confidence_influential_candidate` flag is a reproducibility-based candidate label, not a causal, mechanistic, exclusive, or definitive biological claim.

Example MATLAB command:

```matlab
addpath(genpath('code/residual_contrast/analysis/wm_roi_influence'));
cfg = struct();
cfg.hrfModelName = 'cHRF';
cfg.contrastName = 'Body_LoadDiff_vs_Face_LoadDiff';
cfg.windowSec = [4.32 8.64];
run_roi_influence_elasticnet(cfg);
```

### Cluster/Local Split

Run heavy nested model fitting and sign-flip inference on JHPCE with plotting disabled:

```bash
cd /users/rwang/fMRI_Neuroimaging
export REPO_DIR=/users/rwang/fMRI_Neuroimaging
export DATA_DIR=/users/rwang
export HRF_MODEL=cHRF
export CONTRAST_NAME=Body_LoadDiff_vs_Face_LoadDiff
export WINDOW_SEC="4.32 8.64"
export MAKE_PLOTS=0
sbatch --export=ALL code/residual_contrast/analysis/wm_roi_influence/run_roi_influence_elasticnet_jhpce.sbatch
```

After the `.mat` file exists, generate light local outputs:

```matlab
addpath(genpath('code/residual_contrast/analysis/wm_roi_influence'));
resultFile = 'data/task_residual/roi_influence_elasticnet/Body_LoadDiff_vs_Face_LoadDiff/cHRF/roi_influence_cHRF_Body_LoadDiff_vs_Face_LoadDiff.mat';
run_roi_influence_local_postprocess(resultFile);
```

Plots are written to a `plots/` subfolder. If all three HRF model `.mat` files are provided, local postprocess also writes `hrf_model_comparison.csv` with descriptive deltas only.

### Downstream ROI Evidence Tables

After cluster outputs have been copied locally, build aggregate downstream tables and figures across completed ROI influence runs:

```matlab
addpath(genpath('code/residual_contrast/analysis/wm_roi_influence'));
outputs = run_roi_influence_downstream_analysis();
```

The downstream driver scans `data/task_residual/roi_influence` first, then `data/task_residual/roi_influence_elasticnet`, and writes:

- `roi_influence_all_runs_downstream.csv`: one row per parcel per run with the five evidence flags.
- `roi_influence_ranked_candidates.csv`: all rows sorted by convergent evidence.
- `roi_influence_hrf_consensus.csv`: one row per parcel/contrast summarizing support across HRF models.
- `roi_influence_network_summary.csv`: descriptive network-level summaries when network labels are available.
- `roi_influence_run_inventory.csv`: run metadata and thresholds used by the downstream pass.

For case-by-case inspection by HRF model, including five-evidence visualizations and Canlab orthviews ROI distribution maps:

```matlab
addpath(genpath('code/residual_contrast/analysis/wm_roi_influence'));
outputs = run_roi_influence_hrf_casewise_analysis();
```

This writes one folder per model under `data/task_residual/roi_influence_elasticnet/casewise_by_hrf/`:

- `cHRF/`
- `cHRFderiv/`
- `sHRF/`

Each folder contains the all-ROI evidence table, top candidates, evidence-specific top/support tables, a five-evidence summary, per-HRF figures, and orthviews images for candidate ROI distributions.
