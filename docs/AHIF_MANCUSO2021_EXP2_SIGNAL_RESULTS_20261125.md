# Mancuso second experiment and signal/noise diagnostic — results

Successful run https://github.com/zuizui0223/knk/actions/runs/45583849255, job 139100303164, commit 58e9ef359fdf2b5528a3b1063493f8895522c9d6. Source hashes and synthetic checks passed; prior Experiment1 exact recency and sign gain reproduced within 1e-10. Artifact 15007327273 SHA256 71618993b486ee77fa35caa713f781bdb62480960e59566c64e993d19a13d85f. Exact aggregate archive validation/ahif_v1/mancuso2021_exp2_signal_summary_20261125.json. Same publication and original inoculum, distinct physical reactors. No new independent publication, no changes to PR8 or original effects ledger/manuscript.

Experiment2: 252 author flag-free day1–6 profiles, 40 unique complete reactors (240 profiles) in six eligible matched contrasts, 648 dependent recency fold-date cells, 864 specificity and signed cells. .3:4 missing, .6:4 only 2 complete reactors; STOP those contrasts. Full 103 genus-taxonomic bins, not absolute cell abundance. All 42 mismatched SequencingType labels from the two source tables fall in selected Experiment2 scope. Strict exclusion of affected reactors leaves zero eligible panels, so the result is conditional on trusting the published composition file. Named bacterial genus subset (74 bins) gives recency +5.3651 pp, sign gain +23.2041 pp, same-direction .928832. This is not resolution of sequencing discrepancy. Original author quantitative ASV diversity analysis used 698 samples rarefied to 6840 reads; 824 flag-free source IDs are not identical to that population.

| Metric | Exp1 | Exp2 |
|---|---:|---:|
| Matched recent vs older gain pp | 8.402604 | 5.159878 |
| BR over old B pp | 10.199666 | 5.917829 |
| BR over BS pp | 1.232000 | -.362094 |
| BR over BSQ pp | 3.245467 | 1.465036 |
| BR over BV pp | 6.025659 | 2.968952 |
| Signed vs signblind gain pp | 16.174600 | 21.656793 |
| Recent+old signed vs old signed pp | 11.252141 | 8.543501 |
| Same-direction weighted gross | .799933 | .910634 |
| Dominant/lag response top1 overlap | .722222/.598380 | .580247/.498457 |

Sign gain positive all six Exp2 comparisons; matched abundance state BR-over-BS negative in mean. The common .3 h^-1 × 1/day contrast gives .948049/.930176 same-direction gross and +29.7130/+23.0216 pp signed gain in Exp1/Exp2. Descriptive replication in same system, not meta-analysis or independent publication. High-frequency 16/day not present in Exp2, cannot replicate its decline. Extreme rates may reflect biofilm/mixing as noted by authors; do not infer extinction or resistance.

## Signal/noise diagnostic (post hoc to Exp1)

At matched .3,.4 rates, target days3–6:

| frequency/day | Crossfit squared signal | Half-contrast variation proxy | signal/variation | full-arm TV | sign gain pp |
|---|---:|---:|---:|---:|---:|
| 1 | .689456 | .042008 | 16.412413 | .774264 | 28.327470 |
| 4 | .556359 | .076016 | 7.318998 | .658193 | 24.249192 |
| 16 | .022487 | .044145 | .509382 | .212361 | .640813 |

Predictor day2–5 ratios 14.602119/8.562454/.465962. Cross-product is not direct ecological signal identification; between-reactor heterogeneity, compositional closure, sequencing variation, and shared reference arms remain. Low sign gain at 16/day is compatible with weaker contrast and ordinary replicate variation, not proof of intrinsically faster sign reversal. In Exp2 .9,1.5 matched rates, 1/4 per day target ratios 7.216685/4.565370, sign gain +20.8425/+18.5754 pp, same-direction .922927/.882491. Predictor ratios 6.968261/8.601122, not monotonically decreasing. No claim of causal moderation or mediation. All results are relative-composition contrasts normalized by gross, not population growth directions. Need independent publication with repeated comparable interventions and calibrated uncertainty to formalize synthesis.
