# Saine: baseline biological summaries add modest rank-response prediction beyond read depth

Source run https://github.com/zuizui0223/knk/actions/runs/45389204461, job 138428134687, commit c4acdf4a9f596f09366503f0b2a3ee67e616cbfc. Original file checks, four source-scorer tests and full restricted-wild scoring succeeded. 446 complete logs, five fixed sites, 999 wild replicates. Same source experiment; no new independent evidence.

| Training year | Response | Incremental held-site R-loss gain over depth | Gain of depth over zero | raw approximate p | Holm across four | improved sites |
|---|---|---:|---:|---:|---:|---:|
| 2020 | named_TV | 0.0225% | -0.3277% | 0.118 | 0.236 | 1/5 |
| 2020 | rank_TV | 1.9792% | -0.8155% | 0.01 | 0.036 | 5/5 |
| 2021 | named_TV | -0.0658% | 0.0650% | 0.371 | 0.371 | 0/5 |
| 2021 | rank_TV | 2.3685% | 0.5170% | 0.009 | 0.036 | 3/5 |

| Year | Response | HC3 joint Wald statistic | Approximate p | Holm across four |
|---|---|---:|---:|---:|
| 2020 | named_TV | 1.52140 | 0.653 | 0.9630000000000001 |
| 2020 | rank_TV | 3.15279 | 0.321 | 0.9630000000000001 |
| 2021 | named_TV | 2.41075 | 0.446 | 0.9630000000000001 |
| 2021 | rank_TV | 8.65284 | 0.017 | 0.068 |

Biological predictors improve rank TV treatment-response R-loss by 1.98% / 2.37% compared with depth-plus-log-context treatment model; approximate Holm .036 in both lanes. Cross-year improvement 5/5 sites; same-year 3/5. Total gain vs zero is 1.180% and 2.873%, respectively. Direct 2021 rank moderation diagnostic does NOT pass Holm (.068). This supports modest prediction within this source, not causal moderation, proof of biology-specific mechanism, or a statistically established rank-vs-named effect difference. Restricted null includes nonzero depth/context effects and assumes independent log errors and conditional linear mean; read count is imperfect detection proxy; initial richness/Simpson/top1 may still reflect measurement. Post-QC conditioning; same heldout 2021 observations in two lanes; no independent replication. Standardized richness/concentration/top1 coefficients are correlated and not separately tested. Do not claim individual modifier or structural buffering.

## Stop rule and independent evidence target

Freeze Saine follow-up as an exploratory prediction lead. No more endpoint fishing. Pascual-García et al. 2025 (Nature Communications 16:3002, DOI 10.1038/s41467-025-57591-2, https://github.com/apascualgarcia/ReplayEcology) is a separate data family: 275 parent communities x four revived endpoints. Their original study already established composition-conditioned reproducibility, so initial-state dependence per se is not novel. No introduced-species contrast in that experiment; separate spontaneous-succession prediction lane. Need directly test identity-free initial state prediction against depth and parent identity. Lear et al. 2024 DOI 10.1111/ele.14493, data DOI 10.5281/zenodo.12528023 has a distinct invasion/disturbance history and final resident profile; original study already tests diversity and timing. Li et al. 2019 DOI 10.1038/s41396-018-0283-x, Dryad 10.5061/dryad.b5rg5cs measures pairwise niche/fitness differences and community impact, but the deposited file description lists derived impact rather than the full underlying resident abundance trajectories. Hu et al. 2024/2025 DOI 10.1038/s41559-024-02618-y, Dryad 10.5061/dryad.8gtht76xz, already links resident dynamics to response; data archive 1.93 GB, not acquired. None independently scored yet. No unsupported broad novelty.

General PR CI on commit c4acdf4a9f596f09366503f0b2a3ee67e616cbfc failed during collection due to missing SciPy despite dedicated workflow success. Add SciPy to dev dependency (not change source scoring) and rerun whole CI. DB remains v0.5 with 161 records, no new causal-bridge or pooled effect. Aggregate validation/ahif_v1/saine2025_depth_biology_v1.json.