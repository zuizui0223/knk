# POEM absolute biomass compensation — observed role opposition, not established causal compensation

2026-11-25 JST. Source article https://doi.org/10.1111/jvs.70026, original ZIP https://doi.org/10.5281/zenodo.10119982. Run https://github.com/zuizui0223/knk/actions/runs/45536266086, job 138943495447, scorer commit 551b2034dbed049e8120764c672654f1bf7a91e5; four synthetic unit tests, exact original ZIP MD5 and SHA256 checks, 5210 source rows, six harvest reconstructions, original persistent plot identities and actual scoring/artifact upload passed. Original parser reused by exact Git blobs from frozen PR4. 49 taxa, 14 sown, 35 nonsown, 40 seeded physical plots (20 per cohort, 5 per each S/F/G/L treatment), 50 original plots including excluded B. Unknown removed according to authors' code. Absolute identified dry shoot mass, g per sampled equal-area quadrat, not total production or root biomass. Cohorts fixed, no climate inference. Same source experiment already registered in AHIF. 
## All role-aligned estimates
| Season | Role | Difference g | CI lower | CI upper | Interval family |
|---:|---|---:|---:|---:|---:|
| 1 | early_sown_group | 0.11025 | -0.18783 | 0.40834 | 1 |
| 1 | other_sown | -1.11439 | -1.55965 | -0.66912 | 1 |
| 1 | nonsown | 1.08910 | -1.27045 | 3.44866 | 1 |
| 1 | total_identified_shoot | 0.08497 | -2.48345 | 2.65339 | 1 |
| 2 | early_sown_group | 2.35951 | -2.09771 | 6.81672 | 4 |
| 2 | other_sown | -10.23399 | -15.93755 | -4.53044 | 4 |
| 2 | nonsown | 3.74302 | 0.74853 | 6.73751 | 4 |
| 2 | total_identified_shoot | -4.13147 | -11.86720 | 3.60427 | 4 |
| 3 | early_sown_group | 2.53610 | 0.16560 | 4.90661 | 1 |
| 3 | other_sown | -0.99691 | -4.45836 | 2.46455 | 1 |
| 3 | nonsown | 0.76694 | -4.47733 | 6.01121 | 1 |
| 3 | total_identified_shoot | 2.30613 | -1.91518 | 6.52744 | 1 |

Season 2 primary was selected because prior rank/name retention differed between cohorts, explicitly exploratory. Strict prespecified directional criterion (early group increases AND another group decreases) FAILED; early-group estimate +2.3595 g with simultaneous interval crossing zero. Do not relabel this criterion as a success. A distinct, simultaneous-interval-supported observed opposition is other-sown −10.2340 g (interval [−15.9375,−4.5304]) versus nonsown +3.7430 g ([.7485,6.7375]). This is a possible recruitment-mediated offset of sown biomass losses, not direct competition/compensation causality or evidence of equivalent total mass. Total −4.1315 g, CI [−11.8672,+3.6043]. Nonsown gain is 36.57% of other-sown estimated loss, descriptive ratio, no causal share interpretation. No independent direct growth or interaction measure. Equal weight across cohorts and three early groups. Shared synchronous control means role-mapped and averaged at plot level before variance calculation; not 3 independent controls. Approximate independent-group Welch-Satterthwaite CI, small n=5 per group. Seasons 1 and 3 intervals unadjusted descriptive sensitivity only.

## Cancellation and held-out signed gap
| Season | Mean within-contrast taxon G | Mean N | Mean C | Ratio C/G | Cross-fit taxon gap | Role pooled gap |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 7.29025 | 1.83939 | 5.45086 | 0.74769 | 1.40552 | 2.22524 |
| 2 | 19.61541 | 4.70321 | 14.91220 | 0.76023 | 8.63424 | 12.81392 |
| 3 | 14.79183 | 2.57071 | 12.22112 | 0.82621 | 4.37566 | -0.61040 |
Each G=sum abs taxon mean difference, N=abs sum difference, C=G−N nonnegative by arithmetic. The ratio of averaged C to averaged G is not the average of per-contrast ratios. In season 2 all six source cohort-by-treatment comparisons have positive taxon cross-fit signed gaps, mean 8.6342 g, versus observed mean C 14.9122 g. Role-aligned pooled signed gap +12.8139 g. Cross-fit learned signs on disjoint plot sets, 256 seeded twofold splits, same persistent-plot partitions across ages and S control reuse; repeated splits are not independent, no p-value or CI. Raw cancellation upward biased by noise and mathematical inequality. Positive cross-fit supports descriptive sign reproducibility, not an unbiased compensation ratio or formal causal inference. Season 3 pooled role gap is negative, despite positive mean taxon gap, hence do not generalize over ages or roles. 

## Source contrast ledger
| Cohort | Season | First group | Taxon G | N | C/G | Cross-fit gap |
|---|---:|---|---:|---:|---:|---:|
| 2020 | 1 | F | 8.07946 | 4.68288 | 0.42040 | -5.02680 |
| 2020 | 1 | G | 9.54485 | 0.01751 | 0.99817 | 6.69236 |
| 2020 | 1 | L | 8.00887 | 1.07268 | 0.86606 | 3.03378 |
| 2021 | 1 | F | 6.76479 | 2.32844 | 0.65580 | 1.26620 |
| 2021 | 1 | G | 5.39476 | 2.37608 | 0.55956 | 0.16331 |
| 2021 | 1 | L | 5.94876 | 0.55873 | 0.90608 | 2.30430 |
| 2020 | 2 | F | 34.47621 | 10.79426 | 0.68691 | 18.03520 |
| 2020 | 2 | G | 14.06243 | 1.71523 | 0.87803 | 3.25238 |
| 2020 | 2 | L | 23.30972 | 8.54274 | 0.63351 | 11.49246 |
| 2021 | 2 | F | 21.63200 | 3.94562 | 0.81760 | 12.35358 |
| 2021 | 2 | G | 11.63251 | 2.60577 | 0.77599 | 2.77777 |
| 2021 | 2 | L | 12.57956 | 0.61564 | 0.95106 | 3.89405 |
| 2020 | 3 | F | 23.58868 | 3.07246 | 0.86975 | 11.04414 |
| 2020 | 3 | G | 14.06298 | 0.74912 | 0.94673 | -0.76071 |
| 2020 | 3 | L | 27.77306 | 5.89636 | 0.78770 | 10.84189 |
| 2021 | 3 | F | 8.04928 | 0.59362 | 0.92625 | 5.73481 |
| 2021 | 3 | G | 9.44155 | 4.91259 | 0.47968 | -1.78507 |
| 2021 | 3 | L | 5.83541 | 0.20011 | 0.96571 | 1.17892 |

Lear 2024 resident-count contrast does not show comparable pooled cancellation; the two experiments use different absolute units, taxonomic grains, interventions and aggregation. Neither a mechanism nor a formal cross-system contrast is established. The new hypothesis is that non-sown recruitment can partly offset loss of intended community members in field assembly, tested by explicit independent interventions on recruitment or a source-matched opportunity covariate in future. Do not infer direct competition, resource release, establishment rates, causal mediation, total-biomass equivalence or first demonstration of ecological compensation. All 18 species vectors, role endpoints and split summaries archived in validation/ahif_v1/poem_field_absolute_compensation_v1.json. Main/PR3–5 unchanged; existing effects unchanged.