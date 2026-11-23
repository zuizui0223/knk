# Post hoc endpoint triangulation: is neighbor solo-GR predictive specifically of priority effects?

**2026-11-24.** This is a post hoc robustness analysis, specified before rescoring any further endpoint. It follows the strict 15-species-out morphology and GR_j tests on the same van Steijn et al. source and must not be presented as independent replication.

## Motivation and source
On the age-matched `Cdelta_total` scale, adding the earlier-neighbor's SOLO growth rate `GR_j` to the later focal's six morphological traits reduced held-out new-neighbor RMSE from 0.2285 to 0.2165. Species-level bootstrap of MSE improvement crossed zero, and a within-guild permutation mapping fraction was 0.04. This is suggestive, not conclusive, and the net neighbor C definition uses solo-focal biomass A_i in two denominators.

The source paper van Steijn et al. (2025), *Ecology* e4528 DOI 10.1002/ecy.4528, itself defined `RII_late=(B_late−B_sim)/(B_late+B_sim)`, which **does not use solo controls**. Its `RII_sim=(B_sim−B_alone)/(B_sim+B_alone)` reflects simultaneous competition. Original scientific finding of trait-dependent competition is not novel to KNK.

## Frozen comparisons

Use the same original 202 matched directed focal×neighbor pots, 15 fully heldout species-role folds, six focal morphological traits, same fold-only standardized `GR_j`, and fixed ridge penalty 10.

For **each outcome**, fit only two prespecified models:
1. focal morphology six traits, no neighbor-specific predictors.
2. the exact same focal traits **plus earlier-neighbor solo GR_j**.

Evaluate:
- `RII_late_total`: priority effect over same partner simultaneously planted, **no solo control in the outcome**.
- `RII_sim_total`: ordinary simultaneous neighbor-vs-alone effect, the negative-control-like comparison for priority specificity.

Report per-pair RMSE, species-equal MSE difference, 20,000-bootstrap of **15 heldout neighbor species** with fixed seed 20261124 and number of species improved. Additionally shuffle GR_j across five species within each guild, 199 times (seed 20261126 and 20261127) to diagnose species-specific mapping; this is **not a randomized intervention p-value**.

Interpretation:
- If GR_j helps both RII_late and RII_sim similarly, the capacity may be a general competitive effect, not formed only by early arrival.
- If only RII_late benefits, stronger *predictive* support for a history-gated effect, still not proof of mechanism; stronger at one endpoint alone is not a direct test of relative improvements.
- If neither benefits, the Cdelta gain may rely on normalization or post hoc choice. Report faithfully.
- GR_j is from post-growth solo plant controls, not an ex-ante field measurement; no external species pool validation.

This endpoint triangulation is kept separate from frozen KNK v28–v38, the submitted-candidate Ecology Letters paper, and Draft PRs #3 and #4.
