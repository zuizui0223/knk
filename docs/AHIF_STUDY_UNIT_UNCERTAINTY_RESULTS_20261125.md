# Study-level directional effects: physical-unit sensitivity results

Verified run https://github.com/zuizui0223/knk/actions/runs/45586032532, job 139107756848, scorer commit 620b6929b731a59112192e24cd2e2d2fb8739e7b. Artifact 15009155441 SHA256 ccc7ac50f40f00be0c71f112504e8147d6700bce63685278feff42c493cfdb84. Original source and prior score replays, synthetic checks and all-ones weighted replay passed. Five existing sources, 27 panels, 4622 dependent fold-date cells. All 999 physical-unit weight draws valid in each source and every single-cluster deletion valid. Exact aggregate and paired covariance: validation/ahif_v1/study_unit_uncertainty_summary_20261125.json. Separate 15-row exploratory study-level ledger: database/ahif_v1/study_unit_effects_exploratory_20261125.csv. All individual derived draws retained in artifact; code and seed reproduce them. No changes to prior effect ledger or manuscript.

Primary effect: TV(signblind)-TV(signed), pp. Values below are percentile sensitivity ranges, NOT confidence intervals or calibrated SE.

| Source | Clusters | Point pp | Weight perturbation range (2.5–97.5 percentiles), pp | Positive draws | Leave-one-cluster range pp |
|---|---:|---:|---:|---:|---:|
| Narwani | 20 | 8.6480 | [0.6456,11.2848] | 986/999 | [3.5774,10.5652] |
| Kansas | 4 | -0.6578 | [-7.2322,4.4909] | 397/999 | [-6.2677,5.0803] |
| POEM | 40 | -4.8550 | [-7.3051,-2.3740] | 0/999 | [-6.1552,-3.7971] |
| van Steijn | 10 | 18.9509 | [16.0756,18.3751] | 999/999 | [17.9588,18.5510] |
| Mancuso E1 | 44 | 16.1746 | [13.2944,17.3862] | 999/999 | [15.3513,18.1361] |

Three positive sources have positive 2.5th percentile and all single-cluster deletions positive; POEM retains negative sign throughout; Kansas unresolved. The van Steijn original point lies outside the perturbation interval, showing reweighting changes center in this nonlinear statistic. Do not interpret the intervals as sampling CIs or positive-draw fractions as probabilities/p-values. No inferential significance. Keep study, lab, and endpoint dependence. POEM/van Steijn potentially related, Kansas original linkage conditional, Mancuso read composition separate. 

Secondary recent+older signed vs older signed forecast: Narwani +8.38 pp [3.96,9.59], Kansas +1.18 [-.03,2.21], POEM +9.46 [7.32,11.69], field +14.98 [13.53,15.56], Mancuso +11.25 [9.42,11.97]. In POEM recency improves a signed historical forecast but remains worse than the equal-mass signblind comparator. This is a useful counterexample to universal sign persistence. Same-direction gross share ranges separately archived. No causal proof of directional reorganization or ecological memory. The first endpoint's comparators differ in sign info only; the second differs in historical recency and is not evidence for sign-specific superiority.

Meta-analysis gate: 0/15 rows eligible for formal inverse-variance synthesis. Conditional perturbation SD not sampling SE; fixed sites/years and small arms, conditional data quality, limited publication/lab independence, differing measurements. Do not compute pooled random-effects estimate. Next external target: Saine 2025 fungal log introduction (446 complete logs, 5 sites, 2019 baseline and only two post-introduction dates; separate one-lag sign-forecast lane, not original 3-timepoint mixture). Harmonize one-lag estimand with previous sources before any comparison, preserve design strata and shared controls; distinguish relative-read composition. No new independent publication or outcome claimed from this turn's eligibility review.
