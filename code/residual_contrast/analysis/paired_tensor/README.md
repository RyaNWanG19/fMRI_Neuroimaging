# Option B: paired residual contrast tensors

`run_paired_tensor` fits signed CP models to a `[ROI x time x subject]`
contrast tensor. The default contrast matches the WM population analysis:

`(2back body - 0back body) - (2back faces - 0back faces)`

The repository WM condition order is body/faces/places/tools at 0-back,
then the same order at 2-back (see `functional_clustering/build_wm_contrast.m`
and `core/PredictContrast_WM.m`). Defaults use `WMcHRF.mat`, all 41 lags,
all retained subjects, ranks 1:6, ten random starts per fit and two random
split-half repetitions per rank. No Tensor Toolbox is needed.

```matlab
addpath('code/residual_contrast/analysis/paired_tensor');
cfg = struct('hrfModelName','cHRF'); % also cHRFderiv or sHRF, independently
results = run_paired_tensor(cfg);
```

A smaller execution check (not a completed rank/stability assessment):

```matlab
cfg = struct('ranks',1:3,'nStarts',3,'maxIter',100,'nSplitHalf',1);
results = run_paired_tensor(cfg);
```

Other paired contrasts, including LANGUAGE phase-specific Story minus Math,
use an explicit source file and weights:

```matlab
cfg = struct('dataFile','path/to/language.mat', ...
    'contrastName','Present_Story_minus_Math', ...
    'contrastWeights',[-1 0 0 1 0 0]);
results = run_paired_tensor(cfg);
```

Custom files must contain `Results` in ROI-time-subject-condition order.
Set `hrfModelName` accurately for custom files. `timeSec` can override the
default `(0:T-1)*0.72`; `roiNames` and `subjectIDs` can supply real labels.
Without IDs, ordinal labels track rows only and cannot establish alignment
across HRF files. No cross-model pooling or paired model comparison is done.

## Data handling

- No across-subject centering or within-parcel standardization: these would
  change the mean-contrast and amplitude questions. Units are preserved.
- A subject with any nonfinite value in any used source condition is excluded
  globally. This is stricter than the ROI-specific omnibus exclusion scheme.
- Optional `Flags` must be condition-by-subject; a used flag other than 1
  excludes the subject.
- Without Flags, an entirely zero ROI-by-time source condition is ambiguous
  because extraction initializes missing responses with zeros. These subjects
  are excluded by default and explicitly recorded. Set
  `excludeAmbiguousZero=false` only when provenance establishes valid zeros.
  A zero *contrast* between valid nonzero responses is always retained.
- These checks cannot detect every missing-response artifact when metadata
  are absent. Source Flags/IDs should be preserved during upstream extraction.

## Outputs and interpretation

Timestamped outputs go under `data/task_residual/paired_tensor` (or
`cfg.outputDir`). The MAT file contains every rank's factors, seed/config,
per-start objective histories and convergence, subject audit, mean contrast,
and split memberships and matching diagnostics. CSVs export rank diagnostics,
subject audit, spatial and temporal weights and subject scores. PNGs show
the rank diagnostics and each rank's components. `saveOutputs=false` disables
all writes; `makePlots=false` skips PNGs.

`A` and `B` columns have unit Euclidean norm; signed amplitudes reside in `C`.
The largest absolute entry of each spatial and temporal column is positive;
the corresponding sign is transferred to C. Components are ordered by C norm.
Their energies are not additive because CP components need not be orthogonal.
The group mean reconstruction is `A * diag(mean(C,1)) * B'`.
Each component assumes a shared spatial map and temporal shape across subjects;
only its amplitude/sign varies by subject. Subject-specific latency shifts can
require multiple components. A component with large subject-score variation
can have a near-zero mean score, so large component norm is not evidence of a
population mean effect.

Energy fit means `1 - SSE/sum(D(:).^2)`, not centered variance explained.
The uncentered SVD benchmark flattens ROI x time into features and shares the
same reconstruction target; at equal rank it has more flexible spatial-temporal
patterns and should fit at least as well as CP.

Split-half checks fit the same rank independently to disjoint random subject
halves. Exact one-to-one component matching maximizes the sum of
`abs(A_left' * A_right) .* abs(B_left' * B_right)` (ranks <= 8).
Reported congruence assesses shared spatial-temporal shapes, not subject-score
reliability or held-out prediction. Split component indices are local to each
half, not matched to the full-data components. Repetitions are descriptive;
related participants would need a grouped splitting scheme before claiming
independent replication. `nSplitHalf=0` skips this check and reports NaN.

Inspect per-start convergence and split convergence before interpretation.
Large cancellation ratio (sum of component norms / reconstruction norm > 10)
is a heuristic warning, not a definitive degeneracy diagnosis. Random starts
do not guarantee a global optimum. The smallest final reconstruction error
chooses the reported start, even if it has not converged; this is flagged.

There is no automatic rank selection, p-value, biological network assignment,
or cross-validated decoder. Choose rank using reproducibility, parsimony and
fit, then validate any scientific claim independently. This module supplements
the population omnibus analysis; it does not replace it.

Verification: `addpath('tests'); test_paired_tensor` tests synthetic recovery,
signed contrasts, missingness, reproducibility, the SVD bound and output IO.

## Local interpretation of completed results

```matlab
resultFile = 'path/to/paired_tensor_results.mat';
report = report_paired_tensor(resultFile, struct('ranks',[1 2]));
```

This does not refit CP or reload the original residuals. The `interpretation`
subfolder contains rank/convergence diagnostics, group-mean reconstruction
figures, per-component parcel/temporal/subject plots, top-parcel and extreme
subject-score CSVs, and factor similarity matrices. Group-mean energy fit is
descriptive and in-sample; it is not cross-validated explanatory power.

For anatomical labels, supply a verified extraction-order mapping:

```matlab
cfg = struct('ranks',[1 2], 'roiTable','verified_roi_mapping.csv');
report = report_paired_tensor(resultFile,cfg);
```

The CSV requires `roiIndex` (1 through number of parcels) and `roiName`;
optional `network` labels are retained in the parcel tables. Without this
mapping, stored labels are used, which may be generic Parcel identifiers.
No anatomical identity is inferred solely from matching atlas size.
