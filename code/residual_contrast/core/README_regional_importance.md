# Regional importance for the existing linear SVM

Stage A fits a linear SVM with C=1 at each lag and contrast. It already saves
standardized feature weights and covariance-based Haufe signal maps. Neither
the absolute weight nor the Haufe map alone is a test of regional importance.
Optional extensions measure held-out decoder reliance and perform nested
SVM recursive feature elimination (SVM-RFE). The original all-feature decoder
is retained alongside the selected-feature decoder for a paired comparison.

```matlab
addpath('code/residual_contrast/core');
cfg = struct('runImportance',true,'importanceRepeats',20,'importanceSeed',2718);
PredictContrast_WM('data/task_residual/WMcHRF.mat', ...
    'outputs/wm_svm_importance_cHRF', cfg);
```

For another task, pass cfg as the sixth argument to
`run_observed_stage_flexible`. Existing calls default to importance disabled.
`NBOOT` controls the existing outer subject bootstrap (default 100). Full runs
can be expensive; this extension has synthetic numerical/integration tests,
not a completed full-data analysis. Choose repetitions and bootstrap count
before interpreting rankings, and check Monte Carlo stability.

## Calculation and outputs

Each training bootstrap includes both conditions of each sampled subject;
both observations of each OOB subject remain held out together. Means,
imputation and feature scaling come from training data only. At each fitted
lag and contrast, for each retained ROI independently, swap the two held-out
condition values within each subject with probability 0.5. Repeat this
perturbation and compute baseline accuracy minus mean perturbed accuracy.
The perturbation uses an independent reproducible RNG stream and does not
change splits or baseline outputs. Linear score updates give exactly the
same predictions as evaluating the perturbed feature matrices, apart from
possible floating-point differences at zero-margin ties.

This paired half-swap measure is tailored to our within-subject task decoder;
it is not ordinary across-subject permutation importance. It preserves each
subject's two regional values but disrupts their condition alignment and
cross-region relationships. Correlated regions can substitute for one another,
or create unrealistic feature combinations under perturbation. Interpret it
as reliance of this fitted model, not unique biological relevance.

- `Reliance_obs`: ROI x lag x contrast x bootstrap accuracy loss.
- `Reliance_mean`: equal-weight mean across eligible bootstrap fits.
- `Reliance_resample_sd`: descriptive variation across bootstrap fits;
  combines fit/split variation and perturbation Monte Carlo variation.
- `Reliance_valid_boot`: number of eligible fits for each cell.
- `regional_reliance_<task>.csv`: long table with those summaries and indices.

Accuracy loss is a fraction: 0.02 means two percentage points. Negative values
are retained (the perturbation improved accuracy). Dropped features and failed
fits are NaN rather than zero; valid-fit counts distinguish missing results
from measured zero reliance. Reused OOB subjects mean bootstrap fits are not
independent replicates. SDs and positive means are not confidence intervals,
p-values, or error-controlled regional discoveries. Importance estimation
does not itself select features; the separate option below does.

## Nested SVM feature selection

```matlab
cfg = struct('runFeatureSelection',true,'runImportance',true);
cfg.featureSelection = struct('innerFolds',5, ...
    'featureCounts',[10 25 50 100 200 Inf],'boxConstraint',1,'seed',31415);
PredictContrast_WM('data/task_residual/WMcHRF.mat', ...
    'outputs/wm_svm_selection_cHRF',cfg);
```

Both extensions are disabled by default. Feature selection is performed
independently at each lag, contrast and outer bootstrap, using only the
training subjects. Inner folds keep both conditions and all bootstrap copies
of a subject together. Training retains bootstrap multiplicity, while inner
validation scores each distinct subject once. These are subject groups, not
HCP family groups; related participants require a family-aware design before
interpreting generalization to independent families.

Each inner training fold learns imputation means, standard deviations and
constant-feature exclusions from scratch. Start with all eligible features,
fit the C=1 linear SVM, rank squared standardized coefficients, discard the
weakest batch to reach the next smaller count, and refit. This is **batched
RFE on a predeclared count grid**, not single-feature-at-a-time RFE or a single
weight ranking. Tied weights favor the lower ROI index for reproducibility.
Counts above availability clip to the eligible count; `Inf` means all eligible
features. No omnibus significance or confidence-band result enters selection.

Choose the count with maximum mean inner-fold accuracy; accuracy ties within
1e-12 favor the smallest requested count. C is fixed (`boxConstraint`, default
1) for every candidate, not tuned. Refit the same elimination path on all outer
training rows, retain the chosen count and evaluate only on OOB subjects.
Unusable inner folds invalidate candidate evaluation; they are not silently
omitted. Fewer distinct training subjects reduce the inner-fold count, recorded
in the audit. A count grid omitting `Inf` forces a subset choice; retain it when
comparing whether elimination helps.

The MAT file adds `FeatureSelection` and `split_list`. Main selection outputs:

- `FeatureSelection.bootstrap{b}`: selected/eligible ROI masks, selected count,
  requested count, weights, held-out accuracy/score difference, optional
  selected-model reliance, inner candidate scores and subject-fold audit.
- `accuracy_obs` / `accuracy_mean`: selected-model OOB accuracy.
- `accuracy_difference_obs` / `accuracy_difference_mean`: selected minus
  baseline accuracy on the same OOB subjects and bootstrap fits.
- `selection_frequency`: selected-fit count divided by successful outer fits,
  including nonselection due to training preprocessing. No successful fits
  produces NaN; an unavailable region gets frequency zero with eligibility
  count zero. `selection_frequency_when_eligible` uses eligible fits instead.
- `regional_selection_<task>.csv`: ROI/lag/contrast selection frequency and
  selected, eligible and valid-fit counts.
- `selection_accuracy_<task>.csv`: baseline and selected accuracy, paired
  accuracy difference and valid-fit counts at each lag and contrast.

Selected-model reliance is stored separately in each bootstrap record.
Nonselected regions are NaN there, not zero measured importance. Selection
frequency measures predictive selection stability; no frequency threshold is
a calibrated regional significance test. Correlated regions may substitute
for one another. Inner tuning accuracy is not held-out performance. Outer
bootstrap results reuse subjects and do not provide independent replicates
for a t-test of the accuracy difference. Report the paired comparison
descriptively unless a suitable inference procedure is added.

Background: [recursive feature elimination](https://scikit-learn.org/stable/modules/feature_selection.html#recursive-feature-elimination)
and [MATLAB linear SVM fitting](https://www.mathworks.com/help/stats/fitcsvm.html).

### Cluster execution

The new launcher runs three HRF models as separate array tasks. The resource
requests (8 workers, 32 GB, 24 hours) are starting estimates, not measured full
run requirements. Create the Slurm log directory before submission. Run a
two-bootstrap pilot on one model first:

```bash
mkdir -p /users/rwang/logs
export REPO_DIR=/users/rwang/fMRI_Neuroimaging
export DATA_DIR=/users/rwang
NBOOT=2 sbatch --array=0 --export=ALL code/residual_contrast/core/run_wm_svm_selection_jhpce.sbatch
```

Then use the full defaults after reviewing pilot runtime and memory:

```bash
NBOOT=100 sbatch --export=ALL code/residual_contrast/core/run_wm_svm_selection_jhpce.sbatch
```

`INNER_FOLDS`, `FEATURE_COUNTS` (space-separated numbers, allowing `Inf`),
`RUN_IMPORTANCE`, `IMPORTANCE_REPEATS`, `OUTPUT_ROOT` and `MATLAB_BIN` can override
defaults. Outputs are isolated by HRF model, job ID and array index. A pilot
with two fits is a smoke check, not a stable regional-selection analysis.
No cluster jobs have been submitted by this implementation.

## Comparing prediction with population inference

Keep these examples in mind:

> A region can have a consistently positive residual across all subjects yet
> contribute little to predicting differences between subjects. A region can
> have an average residual near zero yet predict group membership when one
> group has positive residuals and another has negative residuals.

Those examples concern a between-subject outcome. Our actual SVM instead
decodes paired task conditions. For WM it distinguishes category-specific
load differences (e.g., Body_LoadDiff versus Face_LoadDiff). A positive offset
shared by both conditions may be uninformative, while a consistent positive
paired difference can help condition decoding. Thus nonzero population mean
paired differences and successful decoding are related questions, but regional
significance and multivariable decoder reliance remain different quantities.

Compare only matching HRF models, atlas/ROI order, contrasts and lag grids.
The category 2-back minus 0-back bands are not the same contrast as the
Body_LoadDiff minus Face_LoadDiff interaction decoder. For the interaction,
use the corresponding interaction omnibus/band run. Verify subject cohorts
and exclusions as well; file-local row numbers do not establish identity.

An eventual comparison table should retain ROI identity, lag, selection frequency, mean paired
difference, simultaneous band limits, ROI omnibus adjusted p-value, held-out
reliance and its resample variation/eligible-fit count. The omnibus is an
ROI-level whole-curve result, while bands and current reliance are lag-specific.
Report overlap and disagreement descriptively, without treating agreement as
independent replication on the same data. Avoid selecting significant ROIs or
the best lag on the full dataset before evaluating predictive performance.

## Verification

```matlab
addpath('tests');
test_svm_regional_reliance;
test_svm_feature_selection;
```

The checks compare analytic score updates with actual perturbed predictions,
confirm RNG isolation and reproduce baseline Stage A outputs with importance
enabled/disabled. They also verify constant/all-missing feature handling and
CSV persistence.
The feature-selection checks verify a training-only inner-CV reference,
subject-copy grouping and OOB separation, deterministic sparse-signal recovery,
unavailable features/count clipping, and paired baseline/selected persistence.
