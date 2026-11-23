# Post hoc follow-up after morphology-only species-out transfer: neighbor intrinsic growth rate

**Date: 2026-11-24.** This is explicitly post hoc and motivated by the source sensitivity results. The original morphology-only strict species-out protocol was frozen and run first; it did **not** show robust transfer of unseen early-neighbor six morphological traits. No new ecological mechanism or independent experiment is claimed.

## Observation motivating the test
On the authentic van Steijn et al. (2025) source, DOI https://doi.org/10.1002/ecy.4528, our strict 15-species leave-out evaluation for age-matched `deltaC=C_pre−C_sim` gave RMSE: focal six morphology traits 0.2285; both six morphology vectors 0.2512; both plus 6 absolute trait gaps 0.2411. A sensitivity adding `GR` to the six-trait vectors reduced the gap model RMSE to 0.2178. This does not identify where the predictive gain originates.

**GR measurement:** Original `traits.csv` reports growth rate (`GR`) estimated as solo-grown 10-week total biomass divided by growing days. This measurement is related to the solo control that forms the focal \(C\) denominator. Crucial decomposition: the earlier neighbor's `GR_j` is from a **different** species than focal \(i\), so GR_j is not identical to the solo control \(A_i\) in the outcome formula. However both derive from the same source experiment and may encode intrinsic performance, original trait gradients, and common experimental conditions.

## Frozen follow-up comparisons

Use the identical 202 matched directed species pairs, 15 taxon-guild source records and strict heldout-species split: when species `j` is held out, *all* responses with j in either role are excluded from training; test only rows with `neighbor=j`. No tuning from test outcome.

Primary fixed ridge-penalty **10**, intercept unpenalized. Per fold, fit the z-score of the six morphology traits and separately z-score `log1p(GR)` on **14 training species** only; guild references fixed. Compare:

1. `focal_morph6` (previous strongest primary model)
2. `focal_morph6_plus_focal_GR` (possible solo-denominator coupling control)
3. `focal_morph6_plus_neighbor_GR` (key predictor, no direct focal denominator coupling)
4. `focal_morph6_plus_both_GR` (joint growth strategy)
5. `focal_morph6_plus_neighbor_morph6`
6. `focal_morph6_plus_neighbor_morph6_plus_neighbor_GR`
7. `focal_morph6_plus_neighbor_morph6_plus_abs_gap6_plus_neighbor_GR`
8. `guild_focal_morph6_plus_neighbor_GR` (guild-adjusted sensitivity).

For each: pooled-pair RMSE and species-equal mean MSE; per withheld neighbor MSE. Predefined main contrast is `focal_morph6_plus_neighbor_GR` vs `focal_morph6`. Calculate 20,000 resampled-species bootstrap interval for MSE gain (positive = predictive benefit), 15 species equal weighted and fraction improved. A species-bootstrap interval excludes zero only if the same held-out neighbor species generally benefit; it is not a new-experiment CI.

**Mapping specificity**: shuffle the **15 early-neighbor GR values within the three five-species functional guilds** while leaving each focal species' six morphology traits and its GR unchanged (499 shuffles, seed 20261124). Refit strict species-out `focal_morph6_plus_neighbor_GR`. Report observed species-equal MSE and null distribution, fraction of null trials as or more predictive. Mapping p is not a causal randomized treatment p.

**Caveat checks**: report MSE gain after excluding heldout Antennaria test species (previously an extreme observed neighbor effect), and fixed penalties 1 and 100 as explicitly labeled post hoc sensitivity; do not select the best penalty using test MSE. If possible record neighbor-GR source rank/mean and model coefficient orientation, but do not infer niche-preemption without a resource-access manipulation.

**Interpretation gates:**
- If GR_j does not reliably improve over focal morphology-only or cannot beat shuffled GR mapping, the prior 7-trait gain must not be represented as transferable neighbor growth mechanism.
- If GR_j improves consistently across unseen neighbor taxa, the correct conclusion is a **predictive role of independent solo-growth potential for priority-modified net neighbor impact**, under the same 15-species experiment; it is not identified as a direct competitor suppression coefficient and does not establish a general ecological law.
- Original van Steijn et al. already linked growth/acquisition traits and priority effects in sample. Incremental KNK contribution, if any, is the stricter excluded-species test, not trait–priority-effect association itself.
- Do not change KNK frozen v28–v38, Ecology Letters manuscript, Dawson-Glass / POEM Draft PRs.
