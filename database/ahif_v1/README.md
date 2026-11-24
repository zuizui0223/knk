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
- SQLite `review_candidates` — **268 original Zou & Rudolf 2023 citations and source-coded descriptors**, independently hash-pinned and staged from the public SupportingData.xlsx on GitHub Actions. All have the status `not_independently_screened`; they are **not** included in the curated 26 study records or the 45 scored effects, and are not counted as completed AHIF evidence. Source staging CSV is in the CI artifact; the original workbook is not redistributed.

The SQLite database plus source table CSVs and JSON checks are compiled from the above tracked text inputs, using `benchmarks/build_ahif_db.py`. Any builder failure blocks the artifact. The raw original study datasets (CSV/XLSX from prior researchers) are not copied or redistributed.

## Independent review-data audit and nearest-prior-art screening (2026-11-24)

Original author-reviewed Zou & Rudolf 2023 source `SupportingData.xlsx` (DOI https://doi.org/10.1016/j.tree.2023.08.001, Git blob `be2777d4e4bb1454db393b94a2118ada79118dfa`) contains **268 bibliographic rows, 14 fields**. AHIF audits the original author `LitReview.R` design bin and stores a **numeric-corrected version alongside**, not in place of, the original. One paper by Fukami & Morin (2003) has 18 manipulated species; the original R `case_when` classifies it as "No temporal treatment" because only literal `>10` is allowed. Original bins **46/107/72/43** become numerical **45/107/72/44**. [See audit](../../docs/AHIF_REVIEW268_SOURCE_AUDIT_AND_RICHNESS_GAP_20261124.md).

Exactly **13/268** original review rows have a specific `Fitted Competition Model` entry; **10** of those have manipulated temporal treatments, **none of the 44** numerically recoded 11+ manipulated-species studies has a model code in the review table. The remaining `-` codes **do not mean a lack of actual species interaction measurement** and source paper auditing is incomplete. Thus this is **not** a cross-literature prevalence finding or proof of novelty.

The GitHub Actions reproducible source pipeline now creates:
- `review_candidates`: **268 unreviewed bibliography records** with original source row provenance.
- `review_priority`: **21 unreviewed near-prior-art titles**, the union of 11 title-specific risks and 13 original-review model-coded publications (3 overlaps), with foreign keys to original bibliography rows.
- `review268_design_audit.json`: both original R and corrected numeric categories, raw coded-model counts and all classification discrepancy rows.
- `review268_design_corrected_bins.csv`: all 268 source-linked revised descriptors with original code retained.
- `hypotheses.csv` includes H08, a **falsifiable >=11-species interaction-to-structure generalization hypothesis**; not yet confirmed.

Strongest independent raw-data target: Leopold et al. (2017), DOI https://doi.org/10.1111/ele.12803. Source original Dryad DOI https://doi.org/10.5061/dryad.7p2cv lists exact three CSV names, but source acquisition **failed with ordinary HTTP 401/403 access controls** (verified [Actions run](https://github.com/zuizui0223/knk/actions/runs/45316446032)). No new fungal rank E was computed. A fully specified *before-source* analysis protocol is [here](../../docs/AHIF_LEOPOLD2017_PRE_SOURCE_IDENTITY_RANK_CONTRACT_20261124.md).

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
- The 268 studies in Zou & Rudolf 2023 are an **external review-sourced, unreviewed bibliography staging table**, *not* 268 AHIF-curated effect sizes or verified independent database studies. Our targeted curated publication records were purposefully selected and cannot estimate literature-wide prevalence.
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
- `ahif_v1.sqlite` (including the separate unreviewed bibliographic review candidate table when compiled by CI)
- `ahif_v1_validation.json`
- `ahif_v1_gap_report.md`
- CSV copies, suitable for review before importing to other analytical systems.

CIs and p-values are not automatically pooled. No global AHIF effect estimate is computed. No publication's raw data is bundled.

## Nearest-prior-art challenge audit

See [AHIF novelty defense](../../docs/AHIF_NOVELTY_DEFENSE_AND_PRIOR_ART_20261124.md) for explicit **claims that are already known** and close studies that could defeat a novelty claim (e.g. Leopold et al. 2017 arrival×fungivory×nitrogen; Yan et al. 2024 root barriers; 2021 empirical assembly-graph theory; 2026 multi-factor invasion resistance). These are marked **paper-only** until the original pot/plot-level data are source audited. It is not legitimate to advertise a global absence of a joint same-panel study from this selected pilot.

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
