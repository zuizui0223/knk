# AHIF-DB v0.1 — Assembly History × Interaction Formation

**An original curated ecological evidence database for KNK, not a new systematic review, not an independent new biological experiment, and not a completed meta-analysis.**

## The specific gap this database will test

Most priority-effect studies focus on whether immigration or founder history changes communities or on *why* a pair of species changes competitive response. The pilot DB links **what was manipulated**, **what was measured**, **in what experimental unit**, and **which community state retained that history**.

Proposed ecological hypothesis, NOT yet established:
**History changes realized biotic interaction strengths; the resulting asymmetry determines whether history appears in species identity, unlabeled dominance geometry, spatial process or ecosystem output.**

This is more precise than saying that priority effects exist. Previous reviews and mechanistic experiments have already established priority effects and niche pre-emption/modification (Zou & Rudolf 2023, Stroud et al. 2024, van Steijn et al. 2025); novelty lies in source-audited, same-panel, causal **interaction-to-structural-storage** comparisons and robust conditional predictions.

## What is in v0.1

- `studies.csv` — **publication-level** candidates, DOI or original data URL, author/year/domain, true source-audit status, overlap group and direct source anchor. A study is NOT one statistical replicate.
- `panels.csv` — experimental panel / start-year cohort / separate within-publication experiment, actual biological unit, intervention type (arrival, removal, initial dominance, resident presence, soil solution, inoculum, drought, generation count), endpoints available and **same-panel causal identification**.
- `effects.csv` — numeric **aggregate evidence** with its exact source estimand, CI and p-value where reported, experimental unit, status (frozen/published-source-reanalysis/post hoc) and explicit interpretation limits. Both null-like and positive findings are included. Source effects do not become independent observations merely because there are many endpoints or multiple treatments.
- `hypotheses.csv` — falsifiable hypotheses, necessary intervention/endpoints, what would support/refute, and currently available evidence. A genuine scientific gap has to survive screening, not just be labelled novel.
- `novelty_queries.sql` — reproducible audit queries for the number of independent publication/source clusters, named-vs-rank structural panels, mechanism channels and the **strict same-panel causal bridge**.

The SQLite database plus source table CSVs and JSON checks are compiled from the above tracked text inputs, using `benchmarks/build_ahif_db.py`. Any builder failure blocks the artifact. The raw original study datasets (CSV/XLSX from prior researchers) are not copied or redistributed.

## Evidence classes and allowed claims

| Evidence grade | Required source conditions | What can be claimed |
|---|---|---|
| Candidate only | Paper or source URL, no source-level experimental-unit audit | Screening lead, NOT a confirmed empirical effect |
| Source-audited structure | Independent sample/plot IDs and named + sorted rank scored with design-aware permutations | History related to identity/rank architecture *in that experiment* |
| Source-audited interaction outcome | Pair design, solo/mixed or arrival comparison, focal age/size controls reviewed | History-dependent net impact on measured focal biomass, NOT dynamic alpha unless repeatedly identified |
| Source-audited soil legacy | Source-conditioned solution/inoculum and plant growth responses | Conditioning route in the study; not the identity of a specific microbial or nutrient mediator |
| **Same-panel causal interaction → storage** | Randomized history, partner absence/removal or conditioning/barrier controls, *same units* measured for multi-species named & sorted structure, and identified mediator/interaction effect | Test the hypothesized mechanism-to-storage chain |
| Meta-ready effect | Same intervention/endpoint metric, audited randomization and comparable study-level uncertainty, independent publication clusters with dependence modeled | Design-aware random-effects synthesis after the gate is passed |

**AHIF v0.1 currently contains no verified same-panel causal mechanism-to-structural-memory bridge.** That means no such bridge was *identified among the pilot records*, not that a comprehensive search has proved none exists in the literature. This absence is the reason to conduct a targeted independent systematic screen of literature and new same-panel intervention experiments.

## Source status and dependence

The foundational 16 experiment/panel candidates came from existing **Draft PR #4** (separate macro-meta track). New source-audited interaction datasets were assembled in **Draft PR #5**, e.g. Fragata 2022, van Steijn 2025, Delory 2021, and Jing-family 2023. This DB branch is intentionally **independent of those drafts**; citations point to their exact documents and do not imply the draft effects are in the frozen journal manuscript.

Important source-grain constraints:
- POEM 2020 and POEM 2021 cohorts belong to **one publication** and use repeatedly measured plots; treat neither repeated years nor 2020/2021 cohorts as independent publication replications.
- Alonso rhizobox, POEM, Delory and van Steijn have potentially related research groups: preserve lab/overlap group.
- Jing 2023 Zenodo and Jing 2024 Dryad may share source data; they are put under one provisional publication/overlap cluster **CL_JING_UNRESOLVED** until original authors' data lineage is checked.
- The 268 studies in Zou & Rudolf 2023 are an **external prior review universe**, *not* 268 AHIF-curated effect sizes or independent database studies. Our 21 pilot publication records were purposefully selected and cannot estimate literature-wide prevalence.
- Negative and null-ish measurements are not proof of equality; percentile intervals may rely on author-selected exclusions, source trays, sample/size matching, or questionable exchangeability.

## Reproduction and SQL

Run (local Python 3.11+ / sqlite3 only):

```bash
python -m benchmarks.build_ahif_db \
   --source-root database/ahif_v1 \
   --output-dir /tmp/ahif-db-v1
sqlite3 /tmp/ahif-db-v1/ahif_v1.sqlite \
   < database/ahif_v1/novelty_queries.sql
```

The builder checks exact required columns, unique keys, foreign keys, numeric CIs, p-values, source URLs, publication overlap, source-audit status, response grain and scientifically overstrong causal-path claims. It writes:
- `ahif_v1.sqlite`
- `ahif_v1_validation.json`
- `ahif_v1_gap_report.md`
- CSV copies, suitable for review before importing to other analytical systems.

CIs and p-values are not automatically pooled. No global AHIF effect estimate is computed. No publication's raw data is bundled.

## Primary prior-art references
- Zou & Rudolf (2023), *Trends in Ecology & Evolution*, DOI https://doi.org/10.1016/j.tree.2023.08.001 — 268 reviewed empirical studies, frequency-dependent versus trait-dependent PE framing.
- Stroud et al. (2024), *Trends in Ecology & Evolution*, DOI https://doi.org/10.1016/j.tree.2024.02.004 — cross-disciplinary prior-art/terminology survey.
- van Steijn et al. (2025), *Ecology*, DOI https://doi.org/10.1002/ecy.4528 — actual 15-species priority-effect/competitive-trait experimental precedent; do not claim the trait association is KNK novelty.
- Fukami (2015), *Annual Review of Ecology, Evolution, and Systematics*, DOI https://doi.org/10.1146/annurev-ecolsys-110411-160340 — foundational preemption/modification framework.
- Original KNK source evidence: GitHub Draft PR #3 (Dawson-Glass), #4 (macro-meta, POEM), #5 (interaction and soil history), existing frozen v28–v38 in main.

## Next strictly documented curation waves

1. Read all original rows of the **Zou & Rudolf 2023 SupportingData.xlsx** into a separate **review-candidate staging table**, without marking them source-audited. Deduplicate DOI, publication, experiment, experiment year and lab. Audit the original workbook field/record count before using it.
2. Independently screen 2024–2026 new source papers and their supplementary datasets with a PRISMA-style selection/exclusion ledger. Do not claim *first* or *universal* without reviewing prior studies against this exact joint-path criterion.
3. Prioritize randomized **arrival × competitor presence × preconditioned soil × root barrier**, and multivariate community abundance in the **same experimental plots**, then preregister causal interaction estimands and rank-vs-identity outcomes before scoring.
4. Archive eligible effect sizes with design-aware study-level standard errors, laboratory cluster and potential overlap flags. Only then evaluate a multi-level meta-analytic model, keeping each intervention/endpoint lane separate.

**Do not merge without source-method and author review; existing Ecology Letters v28–v38 remain frozen.**
