# Villarreal donor mismatch results and boundaries — 2026-11-26 JST

Parent PR18 frozen 01075120396027d571517028b4cc2b25eb7c5577. Exploratory mechanism follow-up motivated by the earlier same/different source comparison, not prospective independent validation. Prescore 6d8fdeceeeeebae7ed85b82885fd45c395e504ba; verified run https://github.com/zuizui0223/knk/actions/runs/45716017323, job 139562130689, code ccd4f3b7949aa8b67477d334620d81463c3a569c, artifact 15070086440 SHA256 9fc942311c884b846a693e046ade33438627c192cd052e4d43c757bbdff5bae3. Original source checks, 20 pair crosswalk and previous 2+2 primary factorial results reproduced exactly; all 512 cells per panel defined. Original PR18 data and source-study registry unchanged.

## Pre-transfer mismatch (T2, day 50)

| Group | Bray-Curtis on relative counts | Jaccard | Donor mass in undetected recipient taxa | Donor-only detected taxa |
|---|---:|---:|---:|---:|
| Aa | .427393 | .507088 | .307551 | 6.7 |
| Ab | .858342 | .909216 | .815505 | 15.9 |

TV ranges Aa .331725–.634335, Ab .764901–.941476: no observed overlap of this covariate between the two source groups. T1 corroborates large between-source difference (.378559 vs .845018). 'Undetected' is not absence from the whole tank and donor relative mass is not actual transferred abundance. Treatment origin and compositional mismatch are confounded.

## Within-group diagnostic
Held-group score attribution (same PR18 physical 2+2, repeated split outcome attributed to included recipient pair, NOT independent pair effects): pre-transfer TV vs mass_symmetric, pooled Spearman .715789, within Aa -.369697, Ab +.248485, group-centered Pearson -.238636. One-row deletion Spearman ranges Aa [-.533333,-.133333], Ab [.033333,.483333]; no CI or p-value. Singleton held-recipient diagnostic (2 training / 1 held, own linked donor or disjoint A0 reference as specified), pooled .434586, Aa -.721212, Ab +.187879, centered Pearson -.436723; Aa deletion range [-.816667,-.616667], Ab [-.116667,.400000]. Primary pooled trend cannot be interpreted as within-origin general positive relationship. Novelty mass vs mass_symmetric also fails consistent within-group positivity (Aa -.454545 / Ab -.006061; singleton -.769697 / -.042424). Sign results also archived. No causal moderator or biological rejection established; one site, only 10 pairs/group, shared fold/reference dependencies and detection error. Do not count attributed rows as independent study-level effects.

## Actual donor identity and detection

Risk set = taxa not detected in recipient before interval; new detection at next sample. Enrichment = difference in donor prevalence of new vs still undetected taxa in that risk set. Actual paired donor compared to nine alternative same-lake donors. All weights equal pair (defined subset) then interval. Pre-transfer T1→T2 placebo, three post-transfer pulses. New detections Aa 26/56/49/7, Ab 38/85/75/47 across four intervals; these counts are dependent repeated detections, not unique colonizations. Last Aa interval enrichment defined 6/10 pairs, all others 10/10. 

| Group | Actual donor post-pulse enrichment advantage vs random other donor | Post minus placebo advantage | Hard binary Brier reduction vs random other donor | Whole-trajectory reference upper tail (post / post-minus-placebo) |
|---|---:|---:|---:|---:|
| Aa | +.069172 | +.105466 | +.021886 | .150 / .189 |
| Ab | +.045111 | -.012059 | +.044334 | .051 / .593 |

Actual donor enrichment above alternative mean averaged across pulses; the Ab advantage relative to its own placebo is negative, and interval-level advantages vary in sign. Hard 0/1 donor predictor is compared to another random hard predictor, not a calibrated probabilistic baseline. Whole-trajectory reference uses 999 within-origin donor reassignment permutations, not exact randomization p-values (pair assignment mechanism and exchangeability not established); no significant finding claimed. Coverage and other metrics remain in summary. Detection association does not establish immigrant ancestry, survival, establishment or causal facilitation. Source group pool and shared water can contribute. Exact results and 20 diagnostic pair rows archived; full per-pair/interval metrics in artifact. All meta_eligible=no. 

## Separate candidate gate
Bresciani et al. 2025 ISME Journal wraf041 DOI 10.1093/ismejo/wraf041; PRJNA1147197; public author repo luanabresciani/MCC-Soil-experiment commit b330e9b62b0f9eedc8f2a3c95a0fd3738b8fcd05. Article methods explicitly destructive sampling on days 1,5,15,30, three soil mixtures and three unmixed controls at three dilution levels, six planned replicates per condition×date. Author metadata file Microbial/metadata.csv Git blob 92e3004d5bf6c50372e504f9f42b1439b25e88df: 428 unique profiles vs 432 planned, four missing in AxC dilution5 (T1 4, T15 5,T30 5, T5 6). ASV table Microbial/asv_table.csv blob 7da00f75ba3aacfa61ce421d6079e79968b3a2cc, 11,148,940 bytes, located not acquired. Article says 15,000 rarefied reads, verify raw object before scoring. Rep labels are NOT repeated physical microcosms; use separate destructive-sampling population-mean contrast lane, no repeated-trajectory interpretation. Donors were homogenized before experiment within soil×dilution; source-batch dependence, not 432 independent source treatments. T1 donor controls are AFTER mixing and cannot be relabeled pre-transfer baseline. Primary donor mismatch continuous hypothesis cannot be replicated by this source as currently documented. No new biological score or independent publication count for this gate.
