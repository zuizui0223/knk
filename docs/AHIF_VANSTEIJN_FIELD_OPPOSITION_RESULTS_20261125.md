# Van Steijn 2026 field response opposition: compensation depends on reference schedule

2026-11-25 JST. Source article https://doi.org/10.1002/ecy.70291, original Zenodo https://doi.org/10.5281/zenodo.15363824. Successful run https://github.com/zuizui0223/knk/actions/runs/45537486831, job 138947551329, scoring commit 5143c3e29c60caa2fe738acfcaf8da670246a23e. All three synthetic tests, seven original file MD5s, exact original 2832 Transect-row and persistent plot reconstruction, scoring and upload passed. 40 original plots, 2 fixed sites ×5 blocks ×4 histories, 160 repeated plot-year states, 78-column roster (76 detected taxa +2 all-zero sown taxa); 20 sown, 58 nonsown. Same already registered source experiment; no additional study. Source counts are plant pin hits, not mass, physical abundance or recruitment rates. All eight site-year total-hit distributions vary across plots, range 81–327 observed total hits; variable totals do not rule out finite sampling, detection or vegetation layering constraints. Nonsown not necessarily newly arrived. 

## Primary second-season late-synchronous reference

| Role | Difference hits | Simultaneous CI lower | Upper |
|---|---:|---:|---:|
| early_sown_group | 63.250 | 50.670 | 75.830 |
| other_sown | -21.050 | -37.010 | -5.090 |
| nonsown | -40.700 | -55.233 | -26.167 |
| total_hits | 1.500 | -33.727 | 36.727 |

2022 prespecified primary after prior POEM lead: shared SFlate control, within-block mean of two staggered role-aligned contrasts, equal site weighting. Primary POEM-like other-sown negative AND nonsown positive criterion FAILS. Strict early-sown positive and another component negative criterion PASSES. Other-sown AND nonsown both negative. No causal compensation/recruitment mechanism or total-hit equivalence. Four approximate Welch-Satterthwaite simultaneous 95% intervals, five blocks per site, conditional fixed sites. Total difference +1.5 with wide CI. 

## Reference schedule dependence (secondary descriptive)

| Year | Reference | Early | Other sown | Nonsown | Total | Within-contrast taxon C/G | Taxon cross-fit gap |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 2021 | SFlate | 45.15 | -0.50 | -3.80 | 40.85 | 0.4576 | 7.2008 |
| 2021 | SFearly | 12.50 | -33.15 | 4.00 | -16.65 | 0.6649 | 18.5055 |
| 2022 | SFlate | 63.25 | -21.05 | -40.70 | 1.50 | 0.9158 | 118.9188 |
| 2022 | SFearly | 24.15 | -60.15 | 23.00 | -13.00 | 0.8804 | 94.5059 |
| 2023 | SFlate | 45.60 | -27.95 | -15.20 | 2.45 | 0.8261 | 71.6069 |
| 2023 | SFearly | 32.90 | -40.65 | 5.00 | -2.75 | 0.8167 | 71.1989 |
| 2024 | SFlate | 44.35 | -26.05 | -19.30 | -1.00 | 0.8520 | 80.8083 |
| 2024 | SFearly | 35.90 | -34.50 | -1.60 | -0.20 | 0.8296 | 61.5192 |

The **same two staggered histories** yield a positive nonsown difference +23.00 hits in 2022 vs SFearly but −40.70 vs SFlate. The contrast difference is entirely the two reference means' difference: SFearly minus SFlate nonsown = −63.70 hits, total = +14.50 hits. Early and other role components each differ by +39.10 hits because the synchronous baseline is role-averaged across the two germination groups. This is exact linear contrast algebra, not new independent evidence, and not proof of recruitment mediation. The SFearly result was a prespecified secondary sensitivity with unadjusted intervals only. Reference-sensitive opposition prevents describing the two plant studies as a simple replication/failure of a context-free recruitment-buffering law. Distinguish advancing one group from delaying the other relative to its synchronous comparator. 

## Validation of response direction

Primary mean within-contrast species G=141.30, N=11.90, C=129.40 hits, ratio C/G=.915782. Arithmetic magnitude not unbiased biological compensation rate. All 4 site-by-staggered contrasts have positive held-block signed gaps (2 sites × 2 treatments), mean 118.919 hits; pooled role gap 131.088 hits. 256 seeded twofold splits of the original five blocks per site, shared control retained and same splits across years; descriptive, not independent replicates or p-values. Same directional role signs at both fixed sites for SFlate, but larger magnitude in UME; no causal site/climate inference. SFearly sensitivity has positive nonsown 2022 in both sites but not consistent in later years. 

## Interpretation gate

Current evidence: opposing species/group responses in two field experiments; identity of the compensating group depends on comparator/schedule/context. POEM biomass opposition (other-sown loss, nonsown gain) and Van Steijn hit opposition cannot be pooled as the same effect or direct comparison of mechanisms. No treatment manipulation isolating recruitment availability. Need a crossed recruitment exclusion or independently measured pre-existing propagule opportunity, paired with explicit same-date schedule controls. All 32 species vectors, eight role/aggregate contrasts, all secondary intervals and hit coverage are archived at validation/ahif_v1/vansteijn2026_field_opposition_v1.json. Do not count years, comparisons or split draws as independent studies. No main or frozen PR3–5 modifications.