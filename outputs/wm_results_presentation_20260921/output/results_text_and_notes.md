# Results text

The population-average body-versus-face load-difference residual contrast
showed significant temporal omnibus effects in 332 of 489 ROIs (67.9%) for
cHRF, 304 (62.2%) for cHRFderiv, and 114 (23.3%) for sHRF. Each model used
99,999 subject-level sign flips with Holm correction across its 489 ROIs at
alpha = 0.05. The ROI-specific complete-case sample sizes ranged from 402 to
410 subjects. Across models, 104 ROIs passed in all three, and 343 passed in
at least one. These overlap counts describe consistency across model choices.

The centered whole-curve bootstrap sensitivity analysis, using 99,999
replicates, identified 329, 301, and 115 ROIs for cHRF, cHRFderiv, and sHRF,
respectively. Decisions agreed with sign flipping for 99.4%, 99.4%, and 99.8%
of ROIs. A total of 103 ROIs passed both methods in every model. This close
agreement supports stability to the calibration method, while inference
remains conditional on the documented assumptions.

The findings indicate systematic population-average residual contrasts over
the analyzed time course. They do not identify significant individual lags
or establish classifier importance, causal influence, or a significant
difference between HRF models. The correction family is separate within
each model, so the presentation does not claim study-wide control across
all model analyses. Independence was confirmed by the data owner and null
symmetry was accepted as a modeling assumption. The bootstrap provides an
approximate comparison and does not prove those assumptions.

# Speaker notes

1. **Discovery counts:** Each bar counts ROIs with evidence of a nonzero mean
   curve somewhere among the 41 lags. The larger canonical-model counts are
   descriptive and do not establish that those models are better.
2. **Model overlap:** The seven bars are mutually exclusive. The group of 191
   contains only cHRF and cHRFderiv discoveries; it excludes the 104 shared
   by all three. There are 146 ROIs with no discovery in any model.
3. **Mean trajectories:** Left and right V1 pass in all models. Right a24 passes
   only in cHRF. These are illustrative examples selected after testing.
   The curves communicate effect direction and shape, not lag-specific tests.
   Different significant/not-significant outcomes are not a test of model
   differences. The y-axis uses input residual units, not established percent BOLD.
4. **Bootstrap comparison:** The two methods differ in only 3, 3, and 1 decisions.
   Agreement includes both discoveries and nondiscoveries. The sign-flip test
   remains primary; never choose whichever method gives the smaller p-value.
5. **ROI evidence:** Tied adjusted p-values of 0.00489 arise when raw p-values
   reach the 0.00001 Monte Carlo floor. They do not imply identical effects.
   The 12, 9, and 28 near-threshold Monte Carlo flags in the respective models
   warrant caution about fine boundary decisions; this heuristic is not a
   formal confidence interval.

# Small corrections to the supplied methods slides

- Replace **Identify High Influence Regions** with **Identify Regions with
  Reliable Population-Average Residual Contrasts**.
- Use D_i(t) for the vector across ROIs, and D_i,r(t) for one ROI's scalar
  paired difference. The current equation equates those two objects.
- State explicitly that A and B are the body and face **2-back minus 0-back**
  contrasts, rather than individual task conditions.
- Add the Holm correction across 489 ROIs separately within each HRF model.

# Traceability

Results: cHRF job 35821286, cHRFderiv job 35821287, sHRF job 35821288 in
data/task_residual/jhpce_outputs_wm_population_omnibus_full/Body_LoadDiff_vs_Face_LoadDiff.
Atlas labels: local CANlab2018 atlas, matched by ROI index.
Example selection: ROI 1 (left V1), remaining shared ROI with the highest
minimum omnibus S across models (ROI 2, right V1), and cHRF-only ROI with the
highest cHRF S (ROI 122, right a24). This is an illustrative selection.
Time values follow the project's saved vector 0:0.72:28.8 s.
The JSON in build/results.json retains the full extracted numeric precision.
Editable slide chart workbooks round plotted means to nine decimal places.
