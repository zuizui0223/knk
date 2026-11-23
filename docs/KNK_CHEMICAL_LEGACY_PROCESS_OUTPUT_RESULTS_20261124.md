# KNK: interaction formation through an environmental (soil-solution) legacy — independently audited

**2026-11-24 | separate post-submission PR #5; original KNK v28–v38 unchanged.** This is an independent *existing published experiment*, NOT a newly conducted experiment and NOT a demonstration of universal interaction-network formation.

## Why this experiment matters

The earlier KNK plant-pair result (van Steijn et al. 2025) measured the effect on a 10-week focal plant of having a neighbor already established for three extra weeks. It cannot resolve neighbor biomass/resource capture from niche modification. A separate, mechanistically complementary source is Delory et al. (2021), *Plant, Cell & Environment* 44:1215–1230, DOI https://doi.org/10.1111/pce.13999, data/code https://doi.org/10.5281/zenodo.4010162.

The Delory experiment has a **two-phase intervention**:
1. Forb communities, grass communities, or unsown soil **condition soil solutions**.
2. A **single later-arriving plant with no living resident competitor alongside it** receives one solution for 50 days. The responses are root surface density RSD (proxy for exploration), total dry plant weight PDW, and estimated N uptake, as well as other leaf/root traits.

The original authors already reported that soil chemical legacies can affect root foraging of *Dianthus deltoides* without parallel measurable reductions in N uptake, whereas *Festuca rubra* is less responsive. KNK's new-to-repository work is a source-faithful and sensitivity-graded reanalysis and a **direct statistical comparison of process vs output responses**, not discovering environmental-conditioning priority effects for the first time.

## Original source and critical filtering correction

Original source: `BenjaminDelory/Rhizometabolome_2020@v1.0.0`, `data/Data_plants.csv`, Git blob SHA `d11b4b45acb13bd9c240b004332e76b09eb6c3e3`. The released file contains **60 plants ×44 variables**, 2 species×3 treatments×10 plants. All source values and analysis were audited with pinned Git identity in [GitHub Actions](https://github.com/zuizui0223/knk/actions/runs/45240301551).

**Important erratum:** The original source's R notebook `R/Analysis_PlantData_RhizoMetabolome_2020.Rmd`, lines 156–165, explicitly excludes **SampleID 3** (manipulation error: received 50 rather than 25 mL solution on day 2) and **SampleID 22** (author-identified phenotypic outlier). Both are *Dianthus*; resulting analysis **n=58**, with `Dianthus×Forbs n=9; Dianthus×Grasses n=10; Dianthus×Soil n=9; Festuca n=10 per solution`. The earlier KNK 60-plant analysis was a valid *unfiltered sensitivity*, but **not source-faithful replication**. This was corrected and both outputs retained transparently.

## Original-author-filtered effects: 58 plants

Percentage change versus the unsown-`Soil` solution reference within the same focal plant species. Bracketed intervals are **20,000 whole-plant bootstrap percentile CIs of absolute mean differences**, not percentage-change CIs. `RSD` has units cm²/cm³; PDW is grams; Nuptake is grams.

| Later focal / solution source | RSD mean change | RSD absolute difference 95% | Plant dry mass change | N uptake change |
|---|---:|---:|---:|---:|
| *Dianthus deltoides* / Forbs | **−32.66%** | **[−0.08230, −0.03469]** | −8.02% | −3.63% |
| *Dianthus deltoides* / Grasses | −16.30% | [−0.05287, −0.00491] | +2.07% | +7.54% |
| *Festuca rubra* / Forbs | −2.65% | [−0.07364, +0.07177] | +5.27% | −7.25% |
| *Festuca rubra* / Grasses | +4.57% | [−0.04340, +0.06977] | +8.39% | +2.89% |

The conditionally exchangeable individual-plant label-permutation diagnostics for Dianthus RSD differences are p=**0.0007** (forbs) and p=**0.0375** (grasses); **after Holm correction across 12 planned contrasts, p=0.0084 and p=0.4125**, respectively. The full 60-plant sensitivity gives Dianthus Forbs RSD −20.79% with CI of absolute difference [−0.07052,+0.00881] and p=.1171, indicating that source-author filtering matters. These permutation p values **are not guaranteed causal/randomization-based p**: conditioned solutions were administered in grouped trays, and original dataset lacks tray/batch IDs.

Neither PDW nor Nuptake shows a well-resolved difference in the 58-plant simple mean-difference analysis. **This is NOT evidence of equality of mass or nutrient uptake.** Our narrower descriptive inference is that the root-foraging proxy had a clearer and larger response.

## Stronger direct test: does the root-foraging process respond *more* than biomass or N uptake?

Simply observing one statistically discernible process response and one uncertain output response is NOT a valid test of their difference. We therefore specified a **separate exploratory log ratio-of-ratios** in `docs/DELORY_PROCESS_OUTPUT_DIRECT_CONTRAST_POSTHOC_20261124.md`, with all 2 focal species×2 conditioning legacies×2 output responses checked (8 total). Bootstrap retains the dependence of `RSD`, `PDW` and `Nuptake` within the **same sampled plant**.

\[
D_{RSD,output}=\log(\overline{RSD}_{legacy}/\overline{RSD}_{soil})
-\log(\overline{output}_{legacy}/\overline{output}_{soil}).
\]
A negative value means RSD responds **more negatively in proportional terms** than the output; it does not imply that the output is constant.

| Plant / legacy | RSD vs PDW log-response difference (95% CI) | RSD vs N uptake log-response difference (95% CI) |
|---|---|---|
| **Dianthus / Forbs** | **−0.312 [−0.432,−0.204]** | **−0.358 [−0.568,−0.160]** |
| **Dianthus / Grasses** | **−0.198 [−0.313,−0.087]** | **−0.251 [−0.454,−0.054]** |
| Festuca / Forbs | −0.078 [−0.258,+0.091] | +0.048 [−0.210,+0.279] |
| Festuca / Grasses | −0.036 [−0.190,+0.118] | +0.016 [−0.206,+0.220] |

The most important source-robust point is **Dianthus/Forbs**: even **without** the original-author exclusions, the RSD-vs-PDW contrast remains −0.239 (95% [−0.387,−0.086]), and RSD-vs-N uptake −0.249 (95% [−0.485,−0.009]). In contrast, the direct RSD-only mean-difference result becomes statistically uncertain without exclusions. This underscores why comparing two responses **directly** is stronger than comparing their separate p-values.

**Interpretation limit:** the within-plant ratio-of-ratios does not identify a nutrient-use efficiency or causal buffering mechanism, because uptake here is derived from final biomass×tissue N concentration, the outcomes share component measurements, and source solution/tray batches are unobserved. Comparisons were selected after seeing the first analysis, thus represent explicitly exploratory post hoc statistics.

## Biological evidence synthesis: two distinct routes to ecological history

| Route | Source manipulation | Directly measured | Still missing |
|---|---|---|---|
| **Live biotic preemption/competition** | van Steijn et al. 2025: same-aged focal plant, neighbor established for three extra weeks | Age-matched net neighbor effect on focal biomass, dependent on earlier partner identity | Actual resource competition, early-neighbor biomass and separated size effects |
| **Environmental conditioning (niche modification)** | Delory et al. 2021: soil solution from former grass/forb community applied without a live competitor | Root foraging, root traits, biomass, N uptake | Identification of active metabolites, microbial roles, condition-mediated effects on a whole multi-species community |
| **Spatial resource occupancy and coexistence** | Fragata et al. 2022: arrival reversed in two spider mites | Published source shows altered competitive outcome and modelled coexistence; KNK exact 90-box source audit | Directly comparable multi-species architecture measurement |

**Do not add these different interventions' effect sizes into one pooled mean.** Their mechanistic interpretation, organism endpoints, and independent experimental units differ.

### What remains genuinely new and falsifiable for KNK

**Hypothesis: history-to-interaction-to-structure transmission is channel-specific.** Earlier species can store history both by (i) retaining live occupancy and changing direct neighbor effects and by (ii) leaving a conditioned environment that modifies later plants even in absence of live competitors. These pathways can be observable in root-foraging and species-composition states while being muted in aggregate standing biomass or nutrient acquisition.

None of the current studies directly links the strength of **these measured mechanism paths** to **name-labeled versus identity-free rank structural history storage** in the same multi-species community. That is the next causal experiment: randomized *live resident retained vs removed* × *conditioned vs fresh soil solution* × *arrival order* × *root contact barrier* in at least 3-species assemblages, simultaneous measures of absolute species biomass, resource chemistry, root exploration, and named/rank community profiles. Source rank-storage predictions should be formulated before measuring outcomes; otherwise post hoc transformations multiply apparent discoveries.

## Reproducibility

- Source hash audit and pre-score primary targets: `benchmarks/score_delory_chemical_legacy.py` and `docs/DELORY_CHEMICAL_LEGACY_EXTERNAL_PRE_SCORE_20261124.md`.
- Transparent original R source filter correction: `docs/DELORY_SOURCE_AUTHOR_EXCLUSION_AUDIT_20261124.md`.
- Direct RSD-vs-output ratio contrast: `benchmarks/audit_delory_process_output_ratio.py` and its post hoc contract.
- Synthetic validation: `tests/test_delory_chemical_legacy.py`, `tests/test_delory_process_output_ratio.py`.
- Derived aggregate JSON: `validation/interaction_formation/delory_chemical_legacy_filtered_versus_full_v1.json`.
- Original-source and results CI: [source-58 vs 60](https://github.com/zuizui0223/knk/actions/runs/45240301551), [paired process/output test](https://github.com/zuizui0223/knk/actions/runs/45240737274), both succeeded.
- Raw original CSV is **not committed to KNK**; frozen main analysis and manuscript untouched.
