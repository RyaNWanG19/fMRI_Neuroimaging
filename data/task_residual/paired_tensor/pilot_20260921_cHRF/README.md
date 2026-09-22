# Initial Option B cHRF pilot, 2026-09-21

Input: `data/task_residual/WMcHRF.mat`, Results [489 x 41 x 410 x 8].
Contrast: `(2back body - 0back body) - (2back faces - 0back faces)`.
Settings: ranks 1:3, 3 random starts, maximum 500 iterations, relative
objective-change tolerance 1e-7, seed 20260921, one split-half repetition.
No centering or parcel standardization; all 41 time points included.

402 of 410 subjects were retained. Source rows 103, 171, 200, 243, 263,
399, 403, and 407 were excluded for nonfinite source values. The source has
no Flags or real subject IDs. No entirely zero source conditions were found
among the conditions used in the contrast. Ordinal IDs only track source rows.

| Rank | CP energy fit | Uncentered SVD energy fit | Median matched split-half congruence |
|---|---:|---:|---:|
| 1 | 0.021602 | 0.025205 | 0.535322 |
| 2 | 0.040140 | 0.046001 | 0.711675 |
| 3 | 0.055038 | 0.063806 | 0.575316 |

All selected full-data starts converged. No selected full-data solution
triggered the heuristic cancellation-ratio warning. Detailed per-start and
split-fit convergence diagnostics are in `paired_tensor_results.mat`.

The initial 100-iteration run was superseded by this 500-iteration run;
rank 2 and 3 had not converged at the earlier limit. The recorded numbers
and CSV/MAT/PNG files refer to the extended run.

These are descriptive execution/pilot results, not a completed rank-selection
or replication study. Rank 2 had the highest matching score in this single
split, but no rank was selected. Three components reconstruct only about
5.5% of total uncentered contrast energy. This does not directly measure
how much population mean effect or decodable information they capture.

Reproduce from the repository root:

```matlab
addpath('code/residual_contrast/analysis/paired_tensor');
cfg = struct('ranks',1:3,'nStarts',3,'maxIter',500,'nSplitHalf',1);
results = run_paired_tensor(cfg);
```
