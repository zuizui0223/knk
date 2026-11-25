# AHIF: honeybee non-target dropout and cage-grain audit

2026-11-24. Prasad et al. (2026), The ISME Journal, https://pmc.ncbi.nlm.nih.gov/articles/PMC13099264/; original repo https://github.com/Aiswarya-prasad/2025_aprasad_PriorityEffects. Source published derived abundance matrices; not a new biological intervention. Original paper already reports nonfocal strain responses.

## Actual source analysis

22 strains / 12 species, inoculation day 0 and 3, destructive samples day 10. One cage per treatment. Experiment 3 has two shared full reciprocal order control cages and six first-community strain dropout cages. Original Visualize_results.Rmd (blob 54e77aa1be81bb63794762d65da58442d5592, lines 1376–1485) exports absolute counts; dropout branch omits low-detection filter, so explicitly apply authors' 15-sample exclusion file. All shared control IDs and values match. Source blob digests in JSON. After QC: 69 distinct bees across 8 cages. No treatment replication within cage; no bee-level p-values or CIs.

Both strains of the focal species excluded from both arms, leaving 20 strains / 11 species. Source values <=1 copies recoded nondetected. Keep zero-roster bees as separate detection outcomes; for positive rosters normalize within bee, sort separately before averaging for rank profiles, aggregate by species before species-level normalization/sorting. TV=half L1 distance between two cage mean profiles. Conditional positive-roster distances are descriptive, susceptible to treatment-dependent selection; not causal estimates. Median nonfocal copy ratio includes zeros. Estimated gene copies are not biomass.

| Focal | Dropout order | Bees ctrl/drop | Positive ctrl/drop | Strain named TV | Strain rank TV | Species named TV | Species rank TV | Median copy ratio |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Bifidobacterium apousia | Am1-B1 | 9/9 | 9/9 | 0.187 | 0.113 | 0.170 | 0.117 | 0.978 |
| Bifidobacterium apousia | Bm1-A1 | 10/9 | 10/9 | 0.112 | 0.060 | 0.108 | 0.063 | 1.045 |
| Bifidobacterium sp1. | Am4-B1 | 9/9 | 9/9 | 0.171 | 0.105 | 0.135 | 0.110 | 1.002 |
| Bifidobacterium sp1. | Bm4-A1 | 10/8 | 10/8 | 0.126 | 0.111 | 0.126 | 0.121 | 1.640 |
| Lactobacillus apis | Am7-B1 | 9/10 | 9/10 | 0.247 | 0.195 | 0.207 | 0.193 | 1.816 |
| Lactobacillus apis | Bm7-A1 | 10/5 | 10/5 | 0.177 | 0.106 | 0.111 | 0.098 | 1.957 |

All 6 contrasts and unfiltered sensitivity stored in validation/ahif_v1/prasad2026_nontarget_dropout_cage_descriptive_v1.json. Both focal strains omitted to avoid mechanical structural zeros. Small distance is not equivalence; large distance is not proof of a causal effect under one cage per treatment. No pooled meta-analysis or promotion to 118 normalized prior effects.

## Novelty correction and contrary evidence

Stuble & Young 2020, Diversity 12:71, https://doi.org/10.3390/d12020071 already documented grass/forb arrival order affects nonsown nonnative cover after four seasons. Torres et al. 2023 https://doi.org/10.5061/dryad.80gb5mkws manipulated removal order/timing and measured nonremoved taxa, a different intervention lane. Prasad et al. already showed nonfocal latecomer responses after dropout. Mere non-target history propagation is NOT novel. Narrow test: fixed-roster changes in non-target dominance architecture independently of focal inclusion, temporal persistence, and replicated independently measured interaction effects predicting architecture.

POEM from KNK Draft #4 (Alonso-Crespo et al. 2025) age 2, 2020 cohort unsown35 rank E=.136875445, p=.1128; named E=.06319444, p=.1471. 2021 cohort stopped because one zero-unsown plot. Thus no independent affirmative replication of delayed non-target rank propagation there yet. Prior source-authenticated field time-change run 45376516869 passed 4 tests, Zenodo metadata HTTP403, no new temporal outcomes. Existing field data unchanged. Keep draft, no frozen main or manuscript edit.
