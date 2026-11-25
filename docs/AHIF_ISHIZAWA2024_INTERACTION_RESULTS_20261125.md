# Ishizawa 2024: low-order interaction assays predict identity and rank structure differently

2026-11-25 JST. Original paper PNAS DOI https://doi.org/10.1073/pnas.2312396121, source version https://doi.org/10.5281/zenodo.8433832. Original seven file MD5s verified; source run https://github.com/zuizui0223/knk/actions/runs/45525669923, job 138908229857, scorer commit 9518c77e883ddae2e8a435666cbc29ec2fb1a003. Four synthetic tests passed and artifact uploaded. Pre-outcome contract docs/AHIF_ISHIZAWA2024_INTERACTION_PRESCORE_20261125.md. Raw 2377 rows, fields strain/system/nspecies/rep/cfu, 448 unique strain/subset cells, complete 127 nonempty subsets of 7 strains. Source replicate labels exist; this first analysis uses cell geometric means and does not yet verify physical flask-level replicate alignment. Observed per-strain replicate counts: one singleton n=6 and six n=9, pairs n=6, trios/four/five n=5, six-member groups six n=5 and one n=12, full seven n=9. No physical-flask-level covariance or CI claimed. No raw source redistribution.

Prediction from only 1–2 or 1–3 member subsets reproduces original author point log CFU predictions with maximum absolute differences 7.993605777301127e-15 and 7.105427357601002e-15. Original published higher-richness outcomes never used to fit. 64 held-target composition types, 35 four-, 21 five-, 7 six-, 1 seven-member. All same strain pool and overlapping small-community assay information. These are a fixed finite collection of targets, not independent experiments or random samples of communities. No p-values or extrapolation beyond this experiment.

| Endpoint | Mono mean error | Pair mean error | Trio mean error | Pair-to-trio improvement | Improved targets |
|---|---:|---:|---:|---:|---:|
| named_TV | 0.291917 | 0.213103 | 0.116690 | 45.242% | 55/64 |
| rank_TV | 0.165258 | 0.179880 | 0.108842 | 39.491% | 52/64 |
| identity_assignment_TV | 0.126659 | 0.033223 | 0.007848 | 76.379% | 29/64 |
| simpson_abs_error | 0.097423 | 0.152278 | 0.077680 | 48.988% | 50/64 |
| top1_abs_error | 0.117439 | 0.157436 | 0.093634 | 40.525% | 47/64 |
| log_sum_geometric_cfu_abs_error | 0.206935 | 0.206937 | 0.087324 | 57.802% | 48/64 |

## Matched error decomposition

For each target named TV error = sorted rank TV error + identity-assignment residual. Pair-to-trio mean named improvement .096413 = .071037 rank error improvement + .025375 identity residual improvement. Descriptive rank component share 73.680%; relative identity-assignment error reduction 76.38% because its pair baseline is already small. This is an exact mathematical decomposition of prediction error, NOT a direct measure of the mechanism or causality. Pair vs mono: named error decreases .291917 to .213103, while rank error increases .165258 to .179880, Simpson error increases .097423 to .152278, top1 error increases .117439 to .157436. Pair model can identify relative winners but mispredict dominance magnitude. Trio improves both; rank improvement occurs 52/64. Trio does not outperform mono on rank at every richness (7-member single target rank .244949 vs mono .163899; likewise total target n=1 not independent evidence). Report all null/worsened outcomes.

| Target richness | n | Mono rank | Pair rank | Trio rank | Pair→trio named improvement | Pair→trio rank improvement |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 35 | 0.172401 | 0.163445 | 0.092471 | 0.096023 | 0.070975 |
| 5 | 21 | 0.156261 | 0.181785 | 0.114210 | 0.097110 | 0.067575 |
| 6 | 7 | 0.156731 | 0.241168 | 0.155156 | 0.103287 | 0.086011 |
| 7 | 1 | 0.163899 | 0.286047 | 0.244949 | 0.047280 | 0.041098 |

Equal-richness-group weighting preserves improvement .085925 named and .066415 rank (absolute mean TV). Original paper already demonstrated higher-order interaction information improves composition prediction and discusses attenuation of pairwise effects; novel narrow descriptive diagnostic is named vs rank error partition, not first interaction-to-community link. No randomized arrival order in main factorial subset design, so this is not an identified history-specific pathway and cannot open AHIF's causal bridge. CFU geometric-mean profiles, not arithmetic per-flask relative abundances or total biomass. Source replicate labels require exact alignment audit for any future physical-flask uncertainty. Next prescore: compare original pair/trio models with equally informed mean-small-subset, additive log-response and unlabelled template controls to distinguish contextual interaction terms from more data and direct smaller-community abundance cues; independent temporal-order study still needed. All targets and six endpoints retained in validation/ahif_v1/ishizawa2024_low_order_interaction_v1.json. DB v0.6 remains 34 studies,161 effect records,6 prediction records pending separate low-order lane integration; frozen main/manuscript unchanged.