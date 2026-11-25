# AHIF routing shape and fixed-site contrast — verified results

Date 2026-11-24. Source van Steijn 2026 Ecology 10.1002/ecy.70291, original Zenodo 10.5281/zenodo.15363824. Run https://github.com/zuizui0223/knk/actions/runs/45382068138, job 138403349241, prescored commit f7a32f29f9645755c0df7625a88737078eb161ba. All seven source MD5s verified, five new unit tests passed, original 40 physical plots from 10 randomized blocks retained. Same publication as previous field analysis. No new biological experiment; exploratory follow-up specified after earlier outcomes but before these computations.

## Rank shape at matched richness and Simpson concentration

| Site | Original rank E (new seed) | Geometric residual E | Holm p | One-dominant residual E | Holm p |
|---|---:|---:|---:|---:|---:|
| MAL | 0.3224 | 0.2084 | 0.0356 | 0.1111 | 0.0708 |
| UME | 0.0324 | -0.0429 | 0.7545 | 0.1203 | 0.0708 |

MAL residual rank-change treatment signal remains under the prespecified geometric baseline (Holm .0356), but not under the alternative one-dominant baseline (Holm .0708); UME neither. Geometric mean squared L2 fit residuals .00535–.00854 across the four site-year sets vs one-dominant .02464–.05438. Geometric is descriptively a closer match; cannot dismiss the alternative on outcome-dependent grounds. No claim that all functions of richness/Simpson have been excluded; these are deterministic derived outcomes, not adjustment for a causal mediator. Negative E indicates below-median null fit, not negative interaction.

## Fixed-site heterogeneity

| Endpoint | Difference norm (original endpoint units) | Approximate wild p | Holm p |
|---|---:|---:|---:|
| nontarget_rank | 0.3196 | 0.04000 | 0.04098 |
| nontarget_richness | 4.4989 | 0.01073 | 0.03220 |
| nontarget_simpson_concentration | 0.2289 | 0.02049 | 0.04098 |

All three fixed-site omnibus contrasts yield approximate Holm p<.05 under centered independent-block symmetric-error wild bootstrap (1024 signs). This is NOT an exact treatment-randomization p-value, and only five blocks/site, so site differences remain exploratory. The two locations are fixed, not randomly sampled regions or randomized climate. No cross-endpoint comparisons of norms in incompatible units.

For scalar FtS–StF change, Simpson difference MAL−UME = -0.3168 (unadjusted 95% Welch interval [-0.5660, -0.0677], two-endpoint Holm 0.0400). The within-site contrast means are -0.3254 vs -0.0085. Richness exchange is not detected (Holm 0.6558). Joint four-schedule richness heterogeneity can therefore reflect calendar timing rather than the reciprocal order contrast. Not direct competition measurement.

## Prior-art constraints and independent source targets

- Locey & White (2013) Ecology Letters, https://doi.org/10.1111/ele.12154: feasible set conditioning on richness and total abundance already explains substantial SAD shape; geometry cannot be treated as novel or automatically mechanistic. Future stronger challenge: conditional integer partitions matching richness, total recorded hits and preferably concentration; distinguish hit-count sampling from individual abundance.
- Izsak & Pavoine (2012), https://doi.org/10.1016/j.ecolind.2011.06.030: rank-curve shape and curvature already established as ecological descriptors. Original SAD review McGill et al. (2007), https://doi.org/10.1111/j.1461-0248.2007.01094.x. Do not claim first identity-free rank-abundance representation.
- Stuble & Young 2020, https://doi.org/10.3390/d12020071, original methods: 75 plots (3 California sites x 5 treatments x 5 replicates) planted November 2011, May 2015 survey; source discussion explicitly cites the same plots' 2014 cover data in Young et al. 2017 https://doi.org/10.1111/rec.12384. Mark those as overlapping evidence for future curation; no second independent experiment by counting two papers. Species-by-plot matrices not acquired. Werner 2016 https://doi.org/10.1890/15-1918.1 from same research team is a separate candidate, precise plot lineage must be verified before declaring independence.
- Torres et al. 2023 Dryad https://doi.org/10.5061/dryad.80gb5mkws, source metadata explicitly includes field non-target cover and mesocosm biomass after invader-removal history. Distinct intervention family from arrival; candidate direct independent validation, source bytes not yet acquired this turn.

New result source JSON: validation/ahif_v1/vansteijn2026_routing_shape_fixed_site_v1.json. Archived as a separate exploratory diagnostic, NOT appended to 141 existing v0.4 effects in this step. No meta-analytic pooling, no new independent publication, no direct interaction-to-structure causal bridge. No frozen main or manuscript changes.
