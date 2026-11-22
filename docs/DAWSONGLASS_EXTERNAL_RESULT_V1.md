# Independent 15-species mesocosm history-storage validation — scored results
**Status:** source-supplied and scored locally, 2026-11-22; **unmerged Draft PR only**. Current KNK v28–v38, manuscript, and submission are unchanged.

## Source and audit
- Dawson-Glass et al. 2025, *Ecology*, DOI https://doi.org/10.1002/ecy.4504. Public source https://doi.org/10.5061/dryad.sbcc2frh4.
- User-provided original CSV `Dawson-Glass_ECY24-0464_assembly_data.csv`; SHA-256 `1514d4501a0bdf7a61ace6c250566c5a0b8e4a24a0f4737646d182ebad7b8230`.
- 360 records = 24 pots × 15 species, 12 ambient order / 12 warm order; 21 recorded zeros. No pot×species missingness, duplicates, negative or nonfinite biomass, mixed treatment codes, or all-zero pots.
- Names retained verbatim from original CSV; do not silently normalize source spellings.

## Primary treatment-retention test
Each pot contributes the 15 named-species biomass **proportions** and the **same values sorted descending** to remove species identity. Bray–Curtis distances, one-way design-aware R², **9,999 same whole-pot label permutations, seed 20261122**. `E=R²obs−median(R²perm)`. Treat unblocked label exchangeability as a study-design assumption requiring author confirmation.

| Endpoint | Observed R² | Null median | E | Permutation p |
|---|---:|---:|---:|---:|
| Species-labeled composition | 0.428959 | 0.033074 | **+0.395885** | **0.0001** |
| Species-unlabeled sorted rank profile | 0.032899 | 0.034745 | **−0.001846** | 0.5256 |
| Original raw biomass matrix | 0.430140 | 0.032697 | +0.397443 | 0.0001 |
| Total biomass | 0.044018 | 0.021279 | +0.022739 | 0.3284 |
| Largest rank share | 0.015222 | 0.020284 | −0.005061 | 0.5654 |
| Shannon entropy | 0.002390 | 0.020457 | −0.018067 | 0.8314 |

`ΔE(labeled−rank)=+0.397731` is a **metric-dependent descriptive contrast**, not a quantitative fraction of ecological information.

**Source-method cross-check:** original raw-biomass Bray–Curtis pseudo-F independently reproduced at **16.60595** (df 1,22), equal to the published **16.61** at two decimals.

## Out-of-pot prediction
Fixed L2 logistic model, penalty 4, preprocessing fitted inside each training fold, each pot held out once.

| Representation | Held-out log-loss | Accuracy | Gain versus training-prior baseline |
|---|---:|---:|---:|
| Species named | **0.300853** | **21/24** | **+0.434854** |
| Species names removed | 0.890278 | 14/24 | −0.154571 |
| Training-fold prior | 0.735707 | — | 0 |

Cross-validated fits reuse highly overlapping training data; no naive independent significance claim is made.

## Biological interpretation and exploratory robustness
- **Identity changes rather than dominance intensity.** Mean relative biomass of *Rumex crispus* changes from 10.7% (ambient) to 42.4% (warm order), whereas *Dactylis glomerata* is 37.2%→20.3% and *Bromus inermis* is 35.8%→23.4%. *Rumex* has top biomass in **0/12** versus **8/12** pots; mean maximal-species share remains similar (0.486 vs 0.466).
- Hellinger–Euclidean representation: labeled E=+0.349980, p=0.0001; sorted E=+0.048753, p=0.0989. Euclidean on raw proportions: labeled p=0.0001, sorted p=0.2681. Hellinger group-dispersion permutation p=0.939.
- Omit one pot (24 repeats, 1,999 permutations): labeled E remains +0.374 to +0.428, all p≤0.001; rank-profile p≥0.2975. Omit *Rumex* and renormalize 14 species: labeled E=+0.072907, p=0.0262 — major leverage, not the sole contributor. These are **exploratory checks**, not preregistered new discoveries.
- The 2025 source paper reports that *Rumex* advances its germination schedule by 19 days under warmed germination and explicitly warns that time to grow cannot be separated from direct competition. No inference of niche pre-emption, niche modification, or absolute competition mechanism is justified.

**Interpretation:** a second independent grassland experiment strongly supports **identity-rich historical storage**, contrasted with the earlier fungal architecture-rich result. One experiment (n=24) cannot prove a universal mechanism or exact absence of rank-profile effects.

## Repository discipline
The actual raw CSV is not committed or redistributed. Keep this result on the separate post-submission analysis branch and PR #3 until scientific and author review. Frozen KNK v28–v38 remain untouched. Detailed numerical machine summary: `validation/dawsonglass_external_postsubmission/result_summary_v1.json`.
