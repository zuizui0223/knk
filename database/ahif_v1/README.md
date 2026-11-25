# AHIF-DB v0.4 — Assembly History × Interaction Formation

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
- `novelty_queries.sql` — reproducible audit queries for independent publication/source clusters, named-vs-rank structure, mechanism channels and the **strict same-panel causal bridge**.
- `prior_art_verified.csv` / SQLite `prior_art_verified` — **7 closely related original papers** checked for substantive methodological novelty risks, each linked to its curated study with a strict DOI and source/raw-status integrity gate. Five are paper/metadata-level and two have outcome-source audits. These are **not seven independent effect estimates**. 
- SQLite `review_candidates` — **268 original Zou & Rudolf 2023 citations and source-coded descriptors**, independently hash-pinned and staged from the public SupportingData.xlsx on GitHub Actions. All have the status `not_independently_screened`; they are **not** included in the curated 30 study/candidate records or the 118 numeric evidence records, and are not counted as completed AHIF evidence. Source staging CSV is in the CI artifact; the original workbook is not redistributed.

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

## 2026 original Tribolium source added to the evidence map

Zou & Rudolf (2026), *Ecology* DOI https://doi.org/10.1002/ecy.70454, original author [GitHub version](https://github.com/hengxingzou/Zou2026Ecology/tree/5f28e5da78853bb31de1e1fa86f20f4892d7f2bb) was independently source-hash verified, **414 real experimental vials** (378 mixed, 36 single) and **42 published fitted competition-alpha parameter summaries**. The AHIF extension restricts to 126 equal-initial-frequency mixed vials and computes all 9 arrival-reversal contrasts for species-identified T abundance share versus identity-free **two-species maximum dominance share**, with 400 whole-vial within-density conditional permutations and 20,000 bootstrap replicates. Results are in [the source report](../../docs/AHIF_ZOU2026_414_VIAL_RESULTS_20261124.md) and `validation/ahif_v1/zou2026_414vials_two_species_audit.json`.

The original paper already established phenology-dependent fitted competitive coefficients. The **same adult counts** used to fit those coefficients are also used to calculate the group-level adult proportions, so coefficient-to-composition correlation is **in-sample**, not independently validated prediction, and two-species top1 cannot demonstrate a >=11-species ecological law. This update makes the pilot **30 curated study/candidate records, 32 panels, 118 aggregate numeric rows, 8 falsifiable hypotheses; 10 original-source-audited studies and 7 labeled/unlabeled-rank panels (6 multispecies or functional-group plus one 2-species top1)**. No same-panel causal mediation to ranked architecture is identified.

## 2026 original 20-species Swedish field source and 40-plot repeat design

The original source of van Steijn et al. (2026), DOI https://doi.org/10.1002/ecy.70291, [Zenodo 15363824](https://doi.org/10.5281/zenodo.15363824) was **actually acquired and authenticated against all seven original file MD5 hashes**. The original source processing `Site!="GAL" & Type=="Transect"` yields 2,832 positive pinpoint species records on the same **40 physical randomized field plots** in 2021–2024. Two sites ×5 blocks×4 independent plot treatment labels, 160 plot-year states but **40 independent original plots**. The record has 76 recorded taxa, 18 observed of original 20 sown, 58 naturally recruited/non-target taxa.

[Original source run with 9,999 site-block permutations](https://github.com/zuizui0223/knk/actions/runs/45332458601) successfully computed all **eight site-year named+sorted abundance-profile contrasts** and both site-specific 4-year single-plot trajectories. All **16** primary siteyear named/sorted E outcomes survive Holm correction; the site 4-year E values are positive in both representations. These are **new AHIF representation tests of an already-published empirical priority experiment**, not independently measured competition-network mechanisms. The sown20 vs nonsown58 predeclared subcommunity analyses are also preserved, including late-year non-target rank signals (e.g. MAL 2024 non-target E_sorted=.4799, unadjusted p=.0002). A group with zero source sown Slow taxa is stopped without pseudocounts (UME 2022-24, 2/5/5 field plots absent). All 4-year repeated observations remain assigned to original physical plot IDs.

The result source report is [here](../../docs/AHIF_VANSTEIJN2026_FIELD_HISTORY_STORAGE_RESULTS_20261124.md), input/source protocol [here](../../docs/AHIF_VANSTEIJN2026_FIELD_PRE_SCORE_CONTRACT_20261124.md), authentic aggregate JSON `validation/ahif_v1/vansteijn2026_field_history_storage_summary_v1.json`. The curated DB now records **54 correlated source-derived numeric observations** from this single 40-plot published experiment (16 primary whole-76 E +4 full 4-year E +2 4-year ΔR² +16 sown20 +16 nonsown58). The fact that some non-target species exhibit history-associated abundance architecture is **not proof those species all arrived after the sown guild** or that direct interspecific interaction coefficients mediate the effect.

The separate `methodological_prior_art.csv` / SQLite `methodological_prior_art` stores **5 nearest research-method threats** (sPriority, an in-press priority-metric methods manuscript whose DOI is not yet known, Blonder et al. 2024 LOVE 10.1111/ele.14535, original Swedish 20-species field study, and a 2026 microbial founder-size experiment). These are method references, **not independent ecological effects**. Novelty cannot rest on general priority-effect meta-analysis or predicting assembly from experiments because those concepts already have substantial prior art.

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


## Carlström 2019 metadata-only novelty audit

[Source audit and transmission test](../../docs/AHIF_CARLSTROM2019_METADATA_AND_TRANSMISSION_GATE_20261124.md): seven substantive prior-art checks; unchanged 30 study records, 32 panels, 118 numeric effects, ten outcome-source audits. Original R script plus three metadata blobs checked, 921 filtered sample records (not 921 proven independent units), one R2 treatment-label mismatch identified. No abundance matrices obtained and no new ecological effect promoted. The separate group late-arrival and single-strain dropout assays motivate a resident-community transmission test, conditional on raw-source and overlap checks.


## v0.4 — historical routing (current counts supersede historical sections)

32 study/candidate records, 34 panels, 141 aggregate effect records (118 prior + 21 field tests + 2 previously archived POEM non-target sensitivities), 9 hypotheses and 10 raw-source-audited studies. Eight named/rank panels include seven raw-source-audited panels and one derived-matrix descriptive cage panel. Nine substantive prior-art checks. The 21 field tests reuse the original 40 plots; adjusted p-values are linked via routing_test_metadata.csv and SQLite v_routing_tests. Three exact order-exchange statistics are not E. Separate cage_descriptive_contrasts.csv holds 6 QC and 6 sensitivity contrasts; no inferential p or CI, no promotion to effects. Prasad DOI 10.1093/ismejo/wrag056 links to renamed public repo ID 983583082, commit c994bcb7504db27301e43f37c1c39282fb7e6329. Stuble & Young 2020 is an additional paper-only novelty threat. POEM unsown 2020 age2 sensitivity retained with both null results; 2021 zero-roster stop remains documented. H09 is predictive and unvalidated. No pooling or causal-bridge promotion.


## v0.5 — independent fungal experiment, source-verified negative results
33 study/candidate records, 35 panels, 161 numeric effect rows (141 + 20 Saine tests, including two baseline placebos), 11 source-audited study records. Nine named/rank panels (8 source-audited and 1 descriptive cage panel); 10 paper-level prior-art checks. 460 assigned logs, 446 complete positive QC triplets, 4653 non-target logical OTU columns, 4646 detected after original read-depth QC. 20 results from one new source experiment; only two-year any-introduction named change survives its own four-test Holm family (E=.002368, adjusted p=.014). All rank and among-introduced identity tests fail to reject. Two eligibility tests separate from effects in source JSON. Meta gate closed, no new causal bridge. See docs/AHIF_SAINE2025_SUCCESSION_RESULTS_20261124.md.


## v0.6 — source-linked prediction lane
34 study/candidate records, 36 panels, unchanged 161 intervention/diagnostic effect records; 12 original outcome-source-audited records (not raw sequencing reprocessing), 10 named/rank measured panels including one descriptive cage panel and one prediction-only revival panel. 11 paper-level prior-art checks, 10 hypotheses. New tables prediction_outcomes.csv (6 rows) and prediction_folds.csv (48 rows) separate Saine restricted-wild intervention R-loss from Replay spontaneous terminal SSE. Never pool gains across these estimands or count heldout folds/replicates as independent studies. Replay adds one independent source experiment, not an independent intervention-effect replication. 275 parent communities, 14 collection locations, 1100 endpoints, 1454 ASVs. Named gain .009660 (Holm .007), rank gain .115494 (Holm .002), 13/14 rank folds improve. No direct named-vs-rank contrast yet, no causal moderation, no direct competition coefficient. Baseline persistence and initial full-profile predictors are next falsifiers. See docs/AHIF_REPLAY2025_STATE_RESULTS_20261125.md.

## v0.7 — finite combinatorial interaction-assay prediction lane
35 study/candidate records, 37 panels, unchanged 161 ecological effect/diagnostic records, 13 original outcome-source audited records, 11 named/rank panels, 12 substantive paper checks, 11 hypotheses. Existing six held-origin prediction endpoints and 48 fold rows unchanged. New finite-panel tables `low_order_model_errors.csv` (54 model-endpoint means) and `low_order_model_contrasts.csv` (48 matched contrast-endpoint rows) cover one Ishizawa experiment with 64 overlapping target subset types. These are descriptive, no p-values, no meta pooling. Nine models include equally informed direct means and coefficient controls. Direct trio mean rank TV .084737, original trio .108842; context-dependent extrapolation not necessary in this panel. Source label-level 692 system-rep groups complete; physical flask pairing not independently audited. No history-specific causal bridge. See docs/AHIF_ISHIZAWA2024_CONTROLS_RESULTS_20261125.md.

## v0.8 — absolute response opposition ledger
36 study/candidate records, 38 panels, 161 existing effects unchanged, 14 original outcome-source audits, 13 paper-level checks, 12 hypotheses. New derived tables cancellation_contrasts.csv (44 rows: 18 POEM taxon contrasts, six POEM age aggregates, 20 Lear pooled or level summaries) and cancellation_intervals.csv (27 rows: 12 POEM role intervals, 15 Lear resident intervals) are generated directly from verified JSONs. These contain distinct scales/aggregation scopes, not independent study effect sizes. POEM season2 other-sown and nonsown opposite signs after simultaneous intervals, strict early-guild-increase criterion fails. Lear pooled cancellation weak; ambiguous comparison descriptive only. No causal mediation, equivalence, pooled effect or proof of direct compensation. Existing prediction and low-order tables unchanged. See docs/AHIF_POEM_COMPENSATION_RESULTS_20261125.md and docs/AHIF_LEAR2024_CANCELLATION_RESULTS_20261125.md.

## v0.9 — reference-dependent field response opposition
Same 36 study/candidate records, 38 panels, 161 effects, 14 original outcome-source audits, 13 paper checks, 12 hypotheses. Cancellation ledger extends to 92 contrast summaries and 59 interval rows from three already registered source families; Van Steijn adds 32 site/year/history/reference taxon contrasts, 16 pooled/mean summaries and 32 role intervals. Two field experiments show opposing role responses but do not support context-free recruitment buffering. Van Steijn primary 2022 vs late synchronous: early +63.25, other sown -21.05, nonsown -40.70 hits; strict early opposition passes, POEM-like nonsown gain fails. Same staggered histories vs early synchronous give nonsown +23.00 hits (secondary unadjusted). Reference-dependent difference is exact contrast algebra; no independent replication or causal recruitment mediation. Pin hits, grams and resident colony counts remain distinct units. See docs/AHIF_VANSTEIJN_FIELD_OPPOSITION_RESULTS_20261125.md. No meta pooling or main changes.

## v1.0 — separately indexed curated synthesis and temporal forecasts
No increment to 36 original registry records, 38 panels, 161 ecological effect rows, 14 original outcome-source audits or 92/59 cancellation summaries/intervals. One separate secondary synthesis source (Ladouceur 2020), 12 experiment labels (two ambiguous STOP), 3,393 plot-state and 58,426 taxon-state archived source rows, not independent study counts. New derived tables curated_seed_experiments.csv (12), curated_seed_comparisons.csv (59 including one unmatched-arm stop, 58 scored), curated_kansas_forecast.csv (81 model/time summaries including 15 concurrent sign-only rows and 6 aggregate rows). Direct strict two-direction descriptive gate 5/23 eligible cells, indirect 20/30; none of three eligible final direct cells. Kansas selected follow-up, 4 block labels and 15 source age codes: every nonzero method has same positive directional gates at 9,10,12,13,14, but zero magnitude prediction best overall; adding lag control features worsens history-ridge SSE by 9.29%. Exact equal split enumeration adds marginal code13 compared with prior random split weighting, not new independent evidence. No inference of causal temporal recurrence, seed membership, compensation, measurement-method effect, or meta pooling. See docs/AHIF_KANSAS_FORECAST_RESULTS_20261125.md.