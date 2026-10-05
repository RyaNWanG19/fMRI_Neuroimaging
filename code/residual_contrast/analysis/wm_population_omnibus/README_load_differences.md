# Four WM load-difference omnibus tests

`run_wm_load_omnibus` tests each category's population-average **2-back minus
0-back residual curve** against zero over 41 times. Categories and source
condition pairs are body (5-1), faces (6-2), places (7-3), tools (8-4).
The source `Results` must be ROI x time x subject x condition, with dimensions
489 x 41 x N x 8 and the established project condition order.

These are not the six category-by-category interactions and not eight
absolute condition-specific residual tests. Equal mismatch under both loads
still cancels. Nonsignificance does not demonstrate adequate model fit.

## Inference and output

- One sign per subject, shared across every ROI, time and category.
- Fixed complete-case subjects per ROI-category; sample SD is recomputed
  after each sign flip. Statistic is sum of squared one-sample t statistics.
- **Primary Holm family: 489 x 4 = 1956 hypotheses per HRF model.** No
  study-wide correction across the three HRF models is claimed.
- Defaults: 99,999 permutations; alpha 0.05. Independence and null-symmetry
  gates from the existing engine remain required for inference. The symmetry
  assumption concerns load-difference curves, not absolute residuals.
- Optional `runDiagnostics=true` runs centered whole-curve bootstrap with the
  same 1956-test correction family. It is approximate, not an exact test.
- Arrays have 1956 rows, ordered body ROIs 1-489, then faces, places, tools.
  `results.testMapping` and summary columns `category`/`atlas_roi_index`
  identify the rows. `results.categoryCounts` provides the four counts.
- Files: `temporal_omnibus_results.mat`, `temporal_omnibus_roi_summary.csv`,
  and optional `omnibus_sensitivity.csv`. Existing outputs are not overwritten.
- Generic plots are disabled. Do not pass this result directly to plotting
  code that assumes exactly 489 rows; select the category using testMapping.

## Local smoke test (not inferential)

From the repository root, with MATLAB:

```matlab
addpath('code/residual_contrast/analysis/wm_population_omnibus');
for model = ["cHRF","cHRFderiv","sHRF"]
    cfg = struct('hrfModelName',char(model),'smokeTest',true,'nPerm',99);
    results = run_wm_load_omnibus(cfg);
end
```

For a full local run use `smokeTest=false`, `nPerm=99999`, and explicitly set
`cfg.resampling.independentSubjectsVerified=true`,
`cfg.resampling.nullSymmetryAssumed=true`, and a nonempty
`cfg.resampling.verificationNote` only when those statements are warranted.
Real subject IDs can be passed in `cfg.subjectIDs`; otherwise IDs are ordinal
row labels and cannot establish independence or cross-file subject identity.

## JHPCE (three jobs, one HRF model each)

The array launcher uses the existing tested environment/executable setup.
Default data location is `/users/rwang` (the valid MAT files), and MATLAB is
`/jhpce/shared/jhpce/core/matlab/R2025a/bin/matlab`. Override `DATA_DIR`,
`REPO_DIR` or `MATLAB_BIN` if needed. Create logs before submission.

```bash
cd /users/rwang/fMRI_Neuroimaging
mkdir -p /users/rwang/logs
SMOKE_TEST=1 NPERM=99 \
  OUTPUT_ROOT="$PWD/data/task_residual/wm_load_smoke" \
  sbatch code/residual_contrast/analysis/wm_population_omnibus/run_wm_load_omnibus_jhpce.sbatch
```

Full run, only after confirming independent subjects and accepting symmetry:

```bash
INDEPENDENT_SUBJECTS_VERIFIED=1 NULL_SYMMETRY_ASSUMED=1 \
  RESAMPLING_VERIFICATION_NOTE="Independent subjects confirmed by data owner; joint null sign symmetry of load-difference curves assumed" \
  SMOKE_TEST=0 NPERM=99999 \
  sbatch code/residual_contrast/analysis/wm_population_omnibus/run_wm_load_omnibus_jhpce.sbatch
```

Set `RUN_DIAGNOSTICS=1 NBOOTSTRAP=99999` to add the optional bootstrap;
this increases runtime. Request additional Slurm time if needed. Results go
under `data/task_residual/wm_load_omnibus/AllFourCategory_LoadDiff/<model>/job_<id>`.
Full cluster jobs are not submitted automatically by the local implementation.

## Verification

```matlab
addpath('tests');
test_wm_load_omnibus;
test_wm_population_omnibus;
```

Tests check all four condition differences, row mapping, missing-data masks,
observed means/t statistics, pooled Holm adjustment for both resampling
methods, saved metadata, overwrite protection, and inference assumption gates.
