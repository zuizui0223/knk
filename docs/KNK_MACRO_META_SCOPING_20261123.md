# KNK macroecology / meta-analysis scoping programme (2026-11-23)

**Status: four distinct source-study experiments now have scores (three original source endpoints plus POEM repeated field cohorts), all NONPOOLED; NOT a completed effect-size meta-analysis.** This is a NEW branch separate from Dawson-Glass Draft PR #3 and the frozen v28–v38 Ecology Letters submission.

## Question
Across experimental ecological communities, **where in the present does assembly history remain distinguishable**? Compare named taxon or functional-group composition, the same values sorted without identity, within-guild hierarchy, and ecosystem/host-state outcomes. Test whether these modes vary with organism group, actual number of generations, timing treatment, species richness and experimental timescale once enough independent studies are available.

## Published macrosynthesis as a screening baseline, not a KNK discovery
Zou & Rudolf (2023, *TREE*, DOI https://doi.org/10.1016/j.tree.2023.08.001) reviewed **268 experimental studies**. Their published Box 1 / Figure I reported:
- No explicit arrival-time manipulation (initial state proxy): **46**
- Timing manipulated among **1–2 species: 107**
- Timing manipulated among **3–10 species: 72**
- Timing manipulated among **more than 10 species: 43**
- **222** explicitly manipulating timing; among them **58** used graded timing, not just categorical early/late
- Generation lengths: of 257 classifiable studies, **109** lasted <=1 generation, **27** between 2 and 5, and **121** >5 generations.

The source review's code and underlying 45,448-byte SupportingData.xlsx are publicly archived at https://github.com/hengxingzou/Zou2023TrendsEcolEvol. The individual spreadsheet rows have **not** been independently screened in this phase. These 268 reports are an older *review universe*, not 268 meta-effect sizes or 268 new KNK studies.


**Complementary restoration-focused review, not effect-size evidence:** Weidlich et al. (2021), *Restoration Ecology*, DOI https://doi.org/10.1111/rec.13317, systematically screened restoration assembly-history studies (43 eligible, 30 temperate grassland). Its frequently cited positive priority-effect counts address whether *any* treatment effect was reported, **not** whether history is retained in named versus identity-free abundance architecture; they are not new independent KNK outcomes or meta-effect estimates. The review universe also differs from Zou & Rudolf (2023), so the two denominators must not be added together.

## Three independent datasets with matched identity-stripping calculations

| Source | n experimental units | Endpoint | E(labeled) | E(sorted) | Difference ΔE | rank p |
|---|---:|---|---:|---:|---:|---:|
| Leopold & Busby 2020 | 233 | 5 fungi | +.1467 | +.1486 | −.0019 | .0001 |
| Alonso-Crespo et al. 2023 | 19 | 3 grassland functional groups | +.8034 | +.1949 | +.6085 | .0769 |
| Dawson-Glass et al. 2025 | 24 | 15 plant species | +.395885 | −.001846 | +.397731 | .5256 |

The first two are **frozen KNK** v28/v30/v31 results; the third is **unmerged post-submission Draft PR #3**. None of v32–v37 is a new independent study.

Statistic: null-calibrated Bray–Curtis distance partial R², E=R²observed−median(design-respecting permutation R²null). The E values depend on factor design, sample size, taxon dimensions, response grain and null calibration. ΔE is a *descriptive difference*, not a fraction of true information retained or a causal interaction effect.

For Dawson-Glass, the independently recalculated 24-pot n=360-row original CSV (SHA256 `1514d4501a0bdf7a61ace6c250566c5a0b8e4a24a0f4737646d182ebad7b8230`) gives R²labeled=0.4289594, R²sorted=0.0328995, ΔR²=**+0.3960599**. A fresh **5,000 whole-pot stratified bootstrap** (seed 20261123) gives a **percentile 95% interval +0.1889 to +0.5415**. This is within-experiment sampling uncertainty **not a meta CI**, and source-study intervention assignment still needs paper-design verification. Raw CSV is not redistributed.

## Fourth source, repeated-cohort field experiment: POEM (NEW real-data scoring)

**Source archive retrieved & hash-verified 2026-11-23.** The Zenodo v1.0 ZIP https://doi.org/10.5281/zenodo.10119982 was read in ephemeral GitHub Actions runs. Six original Excel shoot-biomass workbooks, **49 identified species** and 5 source treatment levels (S,F,G,L,B) were audited. 50 unique assigned plots across 2 start-year cohorts, each observed for 3 growing seasons. The whole-plot F/G/L/S **primary** history test uses **20 unique plots per cohort and age**; B is excluded as an unseeded natural-succession control. All snapshots reuse the *same plot identities and assignment* and **belong to one publication cluster**.

| Cohort age 2 | Named-species E (p) | Identity-free sorted E (p) |
|---|---:|---:|
| POEM2020 (harvest 2021-06) | +0.1872 (.0038) | **+0.4270 (.0022)** |
| POEM2021 (harvest 2022-06) | **+0.3226 (.0003)** | +0.0230 (.3955) |

Both highlighted results survive Holm adjustment across the 12 primary study-snapshot tests, but **differences in significance alone do not prove an interaction**. A treatment-stratified plot bootstrap of the cohort contrast in \(R²_{named}-R²_{sorted}\) gave **+.5308, 95% +.1331 to +.7336** (2021−2020). This interval is conditional on the **two** tested establishment years, not a climate-year replication CI. Hellinger–Euclidean reanalysis corroborated the opposite directions; F/G/L-only sensitivity agreed. The 2020 rank result weakens if restricted to the 14 intentionally sown taxa. Neither this contrast nor the paper's original data determines whether weather, unsown species recruitment, initial size or interactions were causal.

A no-pseudoreplication sensitivity concatenated all **three annual 49-taxon compositions per plot** into one 147-coordinate vector. Cohort2020 named E=+.1286 (p=.0011) vs sorted E=+.2372 (p=.0004). Cohort2021 named E=+.1928 (p=.0002) vs sorted E=+.0888 (p=.048). The latter is not evidence of complete rank-memory absence.

Full prescored protocol: `docs/POEM_FIELD_PRE_SCORE_CONTRACT_20261123.md`. Real-data results and limitation review: `docs/POEM_FIELD_RESULTS_AND_BOUNDARIES_20261123.md`. Result JSONs: `validation/macro_meta/poem_field_cohort_age_results_v1.json` and `validation/macro_meta/poem_temporal_channel_sensitivity_v1.json`. Source Actions passed: https://github.com/zuizui0223/knk/actions/runs/45085174472 and https://github.com/zuizui0223/knk/actions/runs/45085506739. **No raw ZIP or workbooks are committed.**

**Do not append both start-year cohorts as two independent effect sizes in the future meta-analysis.** They share one experimental publication/research group and geography; the within-source temporal and endpoint contrasts are meaningful, but not independent source clusters.


## Remaining independent candidate datasets (not yet scored)

1. **POEM is now scored** in its own section above; next source to retrieve is Toju/Tucker nectar microbiome or another direct-order randomized grassland.
2. **Toju et al. 2017 floral nectar microbes**. Dryad https://doi.org/10.5061/dryad.0j4b2, sample × OTU metadata XLSX; manipulates initial dominance and tests persistence across floral generations; **separate initial-dominance lane**.
3. **Carlström et al. 2019 synthetic 62-strain Arabidopsis phyllosphere**. DOI https://doi.org/10.1038/s41559-019-0994-z; raw ENA PRJEB32997, code https://github.com/cmfield/carlstrom2019; drop-outs/late introductions are **not** simple full order permutations.
4. **Tucker & Fukami 2014 nectar microbial communities**. Dryad https://doi.org/10.5061/dryad.1cg26. Temperature × arrival timing; source XLSX available but uninspected.
5. **Torres et al. 2023 inverse arrival/removal**. Dryad https://doi.org/10.5061/dryad.80gb5mkws; **field 51 plots and mesocosm 66 pots**, same publication (cluster them), 41/40 field-cover variable ambiguity to audit. Removal of invasive species is **not pooled with initial colonization arrival**.
6. **Weidlich et al. 2017 grassland**. Published experimental report https://doi.org/10.3389/fpls.2016.02008; 64 plots, sowing first PFG × soil/diversity; full raw species×plot source still unverified.

The **16-record** targeted (non-exhaustive) candidate map is `database/macro_meta/candidate_studies_v1.csv`. Pairwise 2-species mechanisms and initial-condition proxies are flagged as *separate lanes* rather than misleading mainline replications.

## Meta-analysis protocol and gates

**Primary population:** studies with manipulated arrival of taxa or functional groups and genuine unit × taxon abundance. **Unit:** whole plant/pot/plot/patch, not species measurements or multiple timepoints on the same unit. Block (plot, genotype, site or year), repeated measures, selection history, taxa, zeros, and intervention target must be source audited.

**Scores per experiment:** same input abundance multiset and dimensionality, original name-indexed and name-erased sorted response; same experimental units and label-permutation stream within a study. Store R²obs, null R² distribution, E, paired ΔE and design-aware study-level uncertainty. A separate predictive sensitivity uses held-out classification **relative log-loss skill**; fold-level preprocessing and model selection must be entirely within training folds. Guild-within-guild rank is a different representation that retains guild identity.

**One publication ≠ multiple independent effect sizes:** build cluster `publication_cluster`; code exposes only one independent study per publication unless justified by separate source designs. For POEM years, use year nested within experiment/study; for repeated plot-year measurements do not count each time point as a separate experiment. Explicitly mark intervention lanes: direct arrival, initial dominance, late introduction, inverse removal, pairwise mechanisms and initial-condition proxies.

**Conservative pooling gate:** need >=5 independent publication clusters in the **same intervention × endpoint-granularity lane**, design-verified original data and paired study SE before a random-effects meta-estimate is even considered. This minimum is NOT enough by itself for stable moderator regression: reserve temperature/taxonomy/time moderators for a larger audited pool (prefer ~10+ independent clusters per moderator). This rule is implemented in `benchmarks/knk_macro_meta_gate.py` and synthetic tests; by design it **does not output an illegitimate pooled estimate** today.

**Current verdict: BLOCKED FOR POOLING.** Four distinct source publications have some scored real-data structural results (two frozen, Dawson-Glass and POEM both in separate drafts), but POEM contributes repeated cohort/age scores rather than one independently comparable study-level effect. The four vary in design, taxon/functional-group grain, and study-level paired uncertainty. This is a source-audited evidence map and new field experiment, **not a random-effects meta-analysis or cross-taxon universal law**. Also, KNK v38 explicitly has not identified which biological interaction generates fungal dominance, and the new Dawson-Glass result cannot identify competition separately from plant age and growth duration.

## References and audit trail

- KNK main result ledger: https://github.com/zuizui0223/knk/blob/main/docs/ECOLOGICAL_RESULTS.md
- KNK paired identity-stripping: https://github.com/zuizui0223/knk/blob/main/docs/eog_identity_storage_rank_abundance_v31_result.md
- New third-experiment (unmerged): https://github.com/zuizui0223/knk/pull/3
- Literature review frame: Zou & Rudolf (2023) https://doi.org/10.1016/j.tree.2023.08.001
- Complementary terminological synthesis: Stroud et al. (2024) https://doi.org/10.1016/j.tree.2024.02.004

**No changes to existing frozen manuscript or experiment results in main.**
