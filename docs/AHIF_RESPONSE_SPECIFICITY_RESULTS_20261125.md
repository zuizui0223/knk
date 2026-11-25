# AHIF matched response specificity — verified results

Run https://github.com/zuizui0223/knk/actions/runs/45573316614 succeeded, job 139065018320, scorer commit 57eefe22285d0787593c8843737f70892162ea07. Source rechecks and synthetic checks passed. Artifact 15003192215 sha256:3ae2e6b7b99310cb87521fb8dc30df10fdfd33024ed80dbb087a7016aab727ca. Exact aggregate JSON archived at validation/ahif_v1/response_specificity_summary_20261125.json. Same four experiments, all 19 panels, 3,470 dependent split-time cells: pond 1,872, Kansas 78, POEM 1,200, field 320; all eligible on common absolute and signed denominators. No new study.

Gain pp (positive = lower TV with specified recent contrast instead of comparator). B older response base; S both-arm state abundance, SQ sqrt thereof, V contrast variance proxy, R recent absolute response; mixtures equal weight.

| Source | BR vs B | BR vs BS | BR vs BSQ | BR vs BV | R vs SQ alone | Signed vs signblind | Recent+old signed vs old signed | Same-direction gross fraction |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Pond | 8.7043 | -0.0317 | 4.1698 | -0.6430 | 6.0821 | 8.6480 | 8.3811 | .678845 |
| Kansas | 1.1802 | -1.2964 | -1.0655 | .4437 | -5.2390 | -.6578 | 1.1815 | .552362 |
| POEM | 17.4666 | .9807 | -3.7548 | -2.9973 | -9.1298 | -4.8550 | 9.4552 | .506058 |
| Field | 14.9963 | 7.1907 | 5.1187 | 5.6604 | 9.0578 | 18.9509 | 14.9819 | .841092 |
| Equal-source mean | 10.5868 | 1.7108 | 1.1171 | .6160 | .1928 | 5.5215 | 8.4999 | .644589 |
| Exclude Kansas | 13.7224 | 2.7132 | 1.8446 | .6734 | 2.0034 | 7.5813 | 10.9394 | .675332 |

Matched state/variance controls sharply attenuate contrast-specific advantage; only field has consistent superiority across all three matched state/variance mixture comparators. R alone is almost tied with SQ on average; R alone has lower TV than BR overall, not a tuned model choice. No broad response-specific superiority supported. Directional signed forecast improves vs same predicted absolute mass with equal +/- split in 2/4 sources; signed mixture improves over older signed alone in 4/4 (different comparison, not sign-specific recency). Opposite and abstain equal-source fractions .269542 and .085869; direction shares weight held gross response, not species counts. Signed forecasts preserve previous estimated treatment-control direction; not proof of a causal species interaction, dynamical memory, or change in interaction coefficients. Relative composition not tested; contrasts use raw measurement units normalized by gross. Dependent folds/time/reused arms, no p-values or formal meta-analysis. Kansas original source linkage conditional; POEM three states, eligible only in this t>=2 protocol. Portal target is temporal difference, unlike future treatment-control contrast. Preserve nulls, don't optimize these sources further. Next priority independent repeated treatment-control multispecies data; preregister exact sign estimand and replicate/block/measurement rules. Main, PR6 and frozen manuscript unchanged.
