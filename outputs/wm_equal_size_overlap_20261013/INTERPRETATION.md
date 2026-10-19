# Interpretation of equal-sized overlap

Matching the number of selected regions substantially narrows the overlap difference seen with significance-based sets. It does not establish that all of the original difference is caused by set size.

For the top 100 regions in each category:

| Ranking | cHRF common count | cHRFderiv common count | sHRF common count |
|---|---:|---:|---:|
| RMS magnitude | 38 | 31 | 37 |
| Standardized statistic S | 40 | 33 | 33 |

Because each selected set contains 100 regions, these numbers also give the percentage of each set shared by all four categories. Under RMS ranking, sHRF is close to cHRF and above cHRFderiv at K=100. Its overlap is therefore not consistently lower when examining the most strongly ranked regions at equal set sizes.

The comparison varies with K. Under RMS ranking, sHRF has the largest four-category intersection at K=50 and K=150 (13 and 67, respectively). At K=300, cHRF/cHRFderiv/sHRF common fractions are 71.0%/70.7%/62.0%. Under S ranking, they are 69.7%/68.7%/64.3%. This indicates a remaining descriptive ranking difference in broader sets, considerably smaller than the original significance-set contrast.

The pairwise Jaccard analysis tells a similar story. At K=100, mean pairwise Jaccard values under RMS are 0.461/0.419/0.451 and under S are 0.484/0.450/0.450. At K=300, sHRF is lower than both canonical models under both rankings.

Suggested presentation text:

> The large overlap difference between significance-based region sets narrows substantially after selecting the same number of regions per category. Among the top 100 regions, four-category overlap is similar across models, particularly when ranking by residual-contrast magnitude. Smaller significant-set sizes therefore contribute to the original pattern. Some overlap differences remain at larger set sizes, but their statistical significance has not been assessed.

Suggested figure:

Use `figures/equal_size_common_fraction.png` for the main comparison. Each panel uses a different ranking, and each colored curve represents one HRF model. At a chosen K on the x-axis, the y-axis is the number of regions common to all four top-K category sets divided by K. Compare models vertically at the same K. These are top-ranked sets without a significance requirement.

Follow-up:

Subject resampling could assess the stability of rankings and the remaining model differences. That analysis must keep the same subjects together across models and categories and recompute mean/SE-based scores in each sample. The current comparison has no bootstrap confidence bands or formal model-comparison p-values.

These results concern spatial ranking agreement, without establishing similar signs, temporal profiles, biological specificity, or overall model fit.
