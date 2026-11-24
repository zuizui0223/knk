# Non-target historical memory appears later: field follow-up and locked analysis plan

2026-11-24. Follow-up to the authenticated van Steijn 2026 experiment (DOI 10.1002/ecy.70291, original data DOI 10.5281/zenodo.15363824). Existing 40 plots across two sites and four years, not new independent replicates. No changes to frozen main or manuscripts. The present re-tabulation is explicitly post hoc relative to the previous scored outcomes.

## Actual computed results from previously archived aggregate values

Source: `validation/ahif_v1/vansteijn2026_field_history_storage_summary_v1.json`, immutable parent commit `77b38236220169fb8b6cd5dd26a9d23d54c19d2f`, blob `27c6ee5c3210bcafadf175648e5cb3a240c53efd`. Original source-verified run https://github.com/zuizui0223/knk/actions/runs/45332458601. Recomputed all 16 non-target named/rank p-values with Holm within this explicitly defined exploratory family; no arbitrary cherry-picking.

| Site | Year | Sown rank E | Non-target named E | Non-target rank E | Named Holm p | Rank Holm p |
|---|---:|---:|---:|---:|---:|---:|
| MAL | 2021 | 0.3387 | -0.0109 | 0.0464 | 1.0000 | 0.9510 |
| MAL | 2022 | 0.2772 | 0.1226 | 0.0525 | 0.2808 | 0.9892 |
| MAL | 2023 | 0.3788 | 0.3149 | 0.1835 | 0.0016 | 0.1310 |
| MAL | 2024 | 0.1618 | 0.3635 | 0.4799 | 0.0016 | 0.0026 |
| UME | 2021 | 0.6733 | 0.0084 | 0.0078 | 1.0000 | 1.0000 |
| UME | 2022 | 0.5729 | 0.3134 | 0.0704 | 0.0016 | 0.9204 |
| UME | 2023 | 0.4687 | 0.2630 | 0.3542 | 0.0120 | 0.0154 |
| UME | 2024 | 0.4836 | 0.0839 | 0.1426 | 0.5768 | 0.5768 |

Both sites had weak non-target rank signal in the first monitored year. The non-target rank signal peaks in year 4 at MAL (Holm p=.0026) and year 3 at UME (Holm p=.0154); it weakens in UME in year 4. First-to-last non-target rank E changes +.4335/+ .1348, whereas sown rank E changes -.1769/-.1896 (MAL/UME); the difference of these changes is +.6104/+.3244. **These contrasts are descriptive, not a statistically tested time-by-treatment effect.** E is null-centered geometry, not information flow or proportion of transmitted causal effect. Non-target species are not all demonstrably late immigrants; seed banks and propagule supply remain possible explanations. The Holm correction is for this family only, not every analysis in KNK.

## Narrow follow-up hypothesis

A manipulated sowing schedule can leave a detectable imprint on the subsequent within-roster dominance dynamics of species that were not themselves sown. Distinguish an indirect community-state response from an identified causal species-interaction pathway. Separate changes in non-target total pin hits, share of all hits, richness, and internal evenness/architecture; do not mistake compositional closure, absences or treatment-induced sampling differences for absolute biomass.

## Frozen before executing NEW endpoint calculations (not preregistered before prior outcomes)

Use original seven MD5-verified Zenodo files and existing source_panels audit. Keep all 40 physical plots, five complete randomized blocks per site, common target20 and non-target58 rosters, no zero-total imputation. Within each roster normalize separately, sort per plot/year, then subtract the 2021 rank vector from the 2024 rank vector in the SAME plot. Use Euclidean distance on signed change vectors (never Bray-Curtis for signed values).

**Primary family:** two site-specific four-schedule randomized-block tests of the non-target rank change vector, 9,999 within-block whole-plot assignment permutations; Holm for two tests. Sharp null: every plot's change outcome invariant to assignment. A significant result means treatment-associated rank dynamics, not independent proof of a shift in E or a causal mediator. An early null and late significant snapshot are not themselves evidence of difference.

**Secondary family:** non-target named change, sown named/rank changes, non-target detected richness, top-one within-roster share, Simpson concentration, total non-target pin hits and non-target share of all hits; 2 sites x 8 endpoints, same assignment stream, Holm across 16. Richness and concentration are possible ingredients of rank changes; these checks are not an identified mediation model. Rank is a relative quantity. Pin hits are not absolute biomass.

**Order-exchange family:** FtS versus StF within the same block, keeping simultaneous treatments out of this contrast. For 2024–2021 non-target rank change, exact conditional swaps of the two asynchronous treatment labels, 2^10 across both sites; a two-sided multivariate statistic sums the per-site squared mean difference normalized by fixed site sum of squared paired differences. Secondary same test on non-target named change and 2023–2021 rank change; Holm across three. Each site has five pairs only; do not interpret as 10 independent sites. Order exchange changes both guild-specific planting dates, so cannot isolate a pure order effect at equal dates. The simultaneous SFearly and SFlate groups are preserved in the omnibus analysis; no contrast is a direct interaction-coefficient manipulation.

**STOP:** source-hash, plot linkage, complete roster, zero total, missing observations or exchangeability failure. No silent deletion. All analyses exploratory in the context of the previously observed trajectory. No random-effects pooling, no new publication cluster, no causal pathway promotion. New outcomes must be separately checked after the workflow runs; summary results above are already calculated.
