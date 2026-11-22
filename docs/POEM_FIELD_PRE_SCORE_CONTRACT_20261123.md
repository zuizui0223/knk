# POEM 2020/2021 independent field experiment: pre-score analysis contract

**State:** schema and source-design audit completed; **this contract is recorded before any KNK history-storage outcome is computed.** Independent extension only; KNK v28–v38 and Dawson-Glass Draft PR #3 remain unchanged.

## Reproducible data source

Alonso-Crespo et al. (2025), *Journal of Vegetation Science*, DOI https://doi.org/10.1111/jvs.70026; Zenodo v1.0 https://doi.org/10.5281/zenodo.10119982. The 6.6MB ZIP is retrieved on the ephemeral GitHub runner and verified against the Zenodo-file **MD5 `3faa5f7d1619a3beae5941ba763b7a26`**, SHA256 `7541ce47b61244f34a6b3b7f1187a2ef16f79bc173d73e6abed206203939ef8f`. It is **never committed or redistributed** by KNK.

Each cohort (POEM2020 and POEM2021) has 25 experimental mixture plots: 5 each assigned S (synchronous), F (forb-first), G (grass-first), L (legume-first), or B (no seed; free succession). Repeated shoot harvests: POEM2020 in 2020–2022 and POEM2021 in 2021–2023, three **growing seasons** per cohort. Each plot has 4 equal-area quadrats in its first 2020 harvest and 2 in subsequent harvests; original dataset contains 6 Excel biomass tables (the first 2020 harvest stores explicit zeros; others largely list positive harvested entries).

## Unit, reconstruction, exclusions

- **Experimental unit** = entire independently assigned mixture plot, never a plant row, quadrat, or repeated survey. Two cohorts are **one publication/source study cluster**.
- Source columns: `Year, Harvest, Plot, Quadrat, Species, SDW_g`. `Plot` identifies numeric site, treatment (S/F/G/L/B) and within-treatment replicate number. The numbered replicate suffix is **not assumed to be a randomized blocking factor**.
- Validate same 25 unique plot IDs/treatment assignments across all three ages in a cohort; check at least 2 equal-area quadrats per plot, nonnegative finite masses, all six expected harvests.
- If the same species has multiple positive records in one plot×quadrat (observed in the archival workbook), **sum fragment masses**. Otherwise a duplicate-deletion rule can discard real biomass. Calculate mean grams per quadrat for each plot×species. Do not treat quadrats as independent.
- Source authors explicitly remove `Unknown`. Exclude it. Store all observed uniquely identified species from 6 tables in a fixed **49-species** union, explicitly zero-fill unlisted biomass species within an *observed quadrat*. Justification: these are harvested-biomass records, with 2020 first harvest fully enumerating 30 species and zeros, whereas later harvests list mainly presences. **Sensitivity must separately restrict to 14 sown species** (defined in the original source study), and report whether full panel zero coding alters interpretation.
- **Exclude B (free succession)** from primary and sensitivity models because B did not receive the same 14-species mixture. Primary uses S/F/G/L = **20 distinct plots per cohort per age**, four 5-plot groups. Complementary F/G/L-only = 15 plots to test timing-order contrasts without S.

## Scores and hypotheses

For each cohort×age cell (six **repeated** experimental states) and for 4-level S/F/G/L primary:
1. Form per-plot **49-named-species** dry-biomass proportions, preserving actual species indexes across cohorts/ages; comparator sorts the identical 49 proportions descending and erases species names. Each has exactly the same biomass multiset and dimensionality.
2. Bray–Curtis distance and Gower-centered one-way factor fit R². Use **9,999 whole-plot history-label permutations per cohort**, keeping 5 group labels of each type, seed **20261123**. **Use the same permutations for names and ranks** and reuse the same permuted assignment across all 3 ages of a cohort.
3. E = R²obs − median(permutation-null R²), with permutation upper-tail p and null central 95% interval. ΔE=E_named−E_sorted is *descriptive only*; it is not Shannon/mutual information, and a smaller p does not quantify an absolute loss of ecological memory.
4. Report 6 paired plots (2 cohorts×3 ages), but **do not pool them** as six independent effects; the independent replication unit is the cohort/plot allocation, and both cohorts belong to one study.
5. **Predeclared sensitivities:** (a) 14 sown taxa named vs sorted (fixed names, same groups); (b) F/G/L-only original allocation; (c) scalar plot total biomass and largest-species share. If an entire plot has zero 14-sown-taxon mass, **STOP** that sensitivity cell rather than removing it silently.
6. Multiple snapshots are not independent p-value replications; show unadjusted and Holm-adjusted descriptive p-values across 12 primary named/sorted tests. Do not declare no-history if a test is non-significant and do not search additional post hoc guild partitions for significance.

## Material interpretability bounds

- The field experiment distinguishes **formation schedules**, not solely competitive niche pre-emption: sowing date, age, composition of unsown recruits, weather and cohorts can all matter.
- Scales are g per harvested equal-area quadrat, *not root biomass*, whole-plot ecosystem production, or biomass of the full 3×3m plot.
- Minirhizotron features record repeated root images; their roots cannot be assigned to named plant species. Do not merge them as simultaneous species-specific root ranks.
- All 6 harvests are from 50 persistent plot allocations (25 per cohort). For S/F/G/L primary there are only 40 independent plots across cohorts, not 120. The study contributes at most **one** source cluster to a later meta-analysis.
- This candidate is not part of the frozen Ecology Letters letter. If the source aggregation cannot be reconciled with the published methods, STOP and report discrepancies rather than inventing positive results.

## Files

Audit code: `benchmarks/probe_poem_archive.py`, `benchmarks/audit_poem_design.py`. Scoring code: `benchmarks/score_poem_field_storage.py` (separate development, not yet scored at this contract freeze). Source metadata only artifacts from GitHub Actions run 45084814558, ZIP source confirmed. Primary results, if successful, will be written as aggregate-only JSON artifact from an ephemeral runner.
