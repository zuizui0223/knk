# Mancuso 2021 DDR64 independent response diagnostic

Successful run https://github.com/zuizui0223/knk/actions/runs/45574399668, job 139068620650, scorer commit bea55e84613b606aab1c6d8c0f3993e11af082aa. Artifact 15003443652 sha256:373ad0762d2ffbe90210e771f2b94f7012a7b4e4d13f42fa029e3679bd5d1880. Both source hashes and exact approved sample join, prior scorer synthetic tests passed. Exact aggregate JSON validation/ahif_v1/mancuso2021_response_summary_20261125.json. Source eLife DOI 10.7554/eLife.67175.

880 genus-bin count profiles, 824 author-approved IDs, 56 technical exclusions. 42 joined sequencing-type mismatches across full sources, zero on selected Experiment1 days1–6; strict sensitivity identical, not independent evidence. 377 valid source profiles in fixed date scope; 44 unique reactors and 264 profiles in eight eligible comparisons, 103 original bins (including unresolved/non-target bins). Original row totals 141–142083: normalized to relative read composition, not cell abundance. Named-genus subset 74 bins, reclosed compositions, retained share .7764–1 across approved scope. Four contrast STOPs due missing complete reactors: rate .1 all three fluctuating arms lack four complete constant reference reactors; rate .2 frequency1 incomplete. Experiment2 excluded before scoring.

Primary: 864 recency fold-date cells, 1152 specificity/signed cells, dependent/repeated. Matched recent vs old TV .345359 vs .429385, gain 8.4026 pp. Historical response H vs sqrt state +8.9830 pp, vs variance proxy +11.1372 pp. Dominant top1 overlap .722222 > lag response overlap .598380, a counterexample to universal dominance/response separation. Same-time matched mixture BR gains over B +10.1997 pp, BS +1.2320 pp, BSQ +3.2455 pp, BV +6.0257 pp. R alone beats sqrt state +5.7774 pp. Signed vs signblind TV .472657 vs .634403 (+16.1746 pp), recent+older signed vs older signed +11.2521 pp. Same-direction gross .799933, opposite .179101, abstain .020966. Direction gain positive 6/8 panels; named-genus sensitivity +16.2338 pp, same-direction gross .805373, recency +8.5992 pp. These are relative read composition differences, so no absolute abundance, total compensation, causal species interaction, or independent fold-level inference.

## Descriptive post-score frequency observation

Only rates .3 and .4 have complete comparisons at all three fluctuation frequencies. Equal-rate means:

| Pulses/day | Same-direction gross | Opposite gross | Signed gain pp |
|---:|---:|---:|---:|
| 1 | .938679 | .060277 | +28.3275 |
| 4 | .929312 | .067950 | +24.2492 |
| 16 | .635022 | .354809 | +.6408 |

Rate .2 frequency16 same .474799, signed gain -1.9488 pp; rate .3 frequency16 signed gain -6.9490 pp; rate .4 frequency16 +8.2306 pp. Preserve all panels. Frequency is experimentally assigned but the prediction-score moderation contrast is exploratory, no formal randomization or independent uncertainty, shared constant arms. High-frequency dilution may yield smaller treatment-control signals, not necessarily biologically faster reversal; signal-to-noise, compositional closure, and taxonomic grain can contribute. Need holdout second experiment or separate publication, response magnitude/noise-aware comparison and compositional sensitivity before interpreting frequency dependence. The original article already showed dilution-frequency-dependent composition/diversity, not a new discovery of the frequency effect itself. New estimand is directional predictability in disjoint physical reactors. Separate compositional lane; no pooling with absolute measurements. No modifications to existing AHIF effect counts, PR6/8 or KNK main/manuscript.
