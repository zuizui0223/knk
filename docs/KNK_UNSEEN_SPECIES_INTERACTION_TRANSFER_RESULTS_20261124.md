# KNK: do traits predict the interaction formed by an unseen earlier-arriving species?

**2026-11-24.** Post-submission exploratory reanalysis of one pre-existing published experiment. **Source experiment/trait-priority relation already published** by van Steijn et al. (2025), *Ecology* e4528, DOI https://doi.org/10.1002/ecy.4528, Zenodo source https://doi.org/10.5281/zenodo.10512662. The new question is **strict out-of-species prediction** rather than a new priority-effect discovery.

## Motivation: named partner effects do not automatically generalize

Earlier KNK PR #5 documented that after controlling focal identity (and even early-neighbor guild), **the identity of the pre-established neighbor** is strongly associated with the focal plant's age-matched biomass response. But this result reused 15 named neighbors. If it merely memorizes a 15-species lookup table, it cannot predict a previously untested early-arriving species.

**Source:** 15 plant species, 5 grasses/5 forbs/5 legumes; 74 solo controls, 206 simultaneous and 206 sequential focal data rows, 202 exact directed-pair matches. The **focal itself was grown 10 weeks** in solo, simultaneous and sequential conditions; the earlier neighbor had 13 weeks to grow. All raw files were checked against original Zenodo v2 MD5 digests and read ephemerally on GitHub Actions. No raw data redistribution.

## Strict transfer test protocol

**Entire species withheld in both roles:** for each target early neighbor `j`, remove every observation involving j as a focal or neighbor from training. Predict only outcomes with `neighbor=j`. This is substantially harder than the previous leave-one-*pair*-out prediction that left the species in training. Focal and neighbor morphology values are standardized using **only 14 training species per fold**; the held-out species' trait values (but not its response) are available. Fixed L2 ridge penalty **10**, no model selection within test folds. Every one of 202 paired directed outcomes is predicted once. This is still within *one* experiment's species pool, not independent flora sampling.

**Outcome:** age-matched `deltaC=C_pre−C_sim`, each `C=(B_mixed−B_focal_alone)/(B_mixed+B_focal_alone)`, where C compares a ten-week focal species to a ten-week-alone control.

### Main prespecified morphology comparison

Six source morphology traits: `SLA, LA, RMF, Plant height, Lateral spread, LDMC`; `GR` excluded from this first validation because solo-plant growth rate is derived from control biomass and risks outcome-related coupling.

| Strict unseen-neighbor model | Out-of-species RMSE |
|---|---:|
| Global intercept | 0.2642 |
| Focal + neighbor functional guild | 0.2483 |
| **Focal six morphology traits only** | **0.2285** |
| Focal plus early-neighbor six morphology traits | 0.2512 |
| Both six morphology vectors and absolute trait differences | 0.2411 |
| Guild plus both morphology vectors and absolute trait differences | 0.2349 |

**Negative scientific result:** Adding the unseen early-neighbor morphology and even pairwise trait distance does **not** reliably outperform the focal morphology alone. For `both_morph_abs_gap vs focal_morph`, species-equal MSE gain is **−0.00528**, bootstrap 95% **[−0.02333,+0.01049]**, 7 of 15 held-out neighbor species improved. A 199-shuffle within-guild trait-vector mapping placebo ranks actual morphology as **31.5%** upper fraction in predictive success. This is not a formal randomized-treatment p-value.

The original source authors used in-sample trait PCA and demonstrated trait dependence; our more stringent failure does **not** refute their result.

## A separate, explicitly post hoc clue: earlier-neighbor growth potential

The sensitivity analysis including `GR` produced better errors. A new protocol was therefore committed **after** the first negative result but **before** follow-up scoring: `docs/VANSTEIJN_NEIGHBOR_GR_POSTHOC_TRANSFER_CONTRACT_20261124.md`.

The focal's GR_i is derived from the focal solo control and could be mathematically coupled to its C denominator. To separate this concern, include only `GR_j` **of the earlier neighbor j**, a distinct species from focal i, with the same six focal morphology traits.

| Model on age-matched deltaC | Out-of-species RMSE |
|---|---:|
| Focal morphology 6 | 0.2285 |
| + focal's own solo GR_i | 0.2296 |
| **+ earlier neighbor's solo GR_j** | **0.2165** |
| + both GR_i and GR_j | 0.2176 |
| + early neighbor morphology 6 + GR_j | 0.2262 |
| + early neighbor morphology 6, gap 6 and GR_j | 0.2166 |

New unknown-neighbor relative RMSE improvement from GR_j: **5.3%**. The species-equal MSE gain is **+0.00552**, 15-species bootstrap 95% **[−0.00122,+0.01426]**, **10 of 15** held-out neighbor species improved. The CI includes zero; **not** a robust generalizable law. The observed GR_j mapping's MSE **0.04688** ranks ahead of 96% of 499 randomized within-guild GR mappings (mapping upper fraction **0.040**), an explicitly post hoc negative-control diagnostic, *not* experiment-randomized proof.

The transfer gain is strongly influenced by `Antennaria dioica`: after excluding its heldout test fold, the species-equal gain drops to **+0.00211**. With fixed penalty lambda=1 or 100 the sign remains positive; only lambda=100 gives a bootstrap CI just above zero, and **we do not switch to that lambda retrospectively**.

### Additional, source-original endpoint triangulation (post hoc)

To show that the GR_j signal is not an artifact of the focal-alone denominator, score the source authors' original `RII_late=(B_late−B_sim)/(B_late+B_sim)` (no solo biomass enters this denominator), and ordinary `RII_sim=(B_sim−B_alone)/(B_sim+B_alone)`, using **identical complete species-out folds**, the two fixed models, ridge lambda=10, and 199 within-guild neighbor-GR mappings each.

| Outcome | Focal morphology RMSE | Plus earlier-neighbor GR RMSE | Relative RMSE improvement |
|---|---:|---:|---:|
| **Original late-vs-simultaneous RII (priority effect)** | 0.25038 | **0.21597** | **13.7%** |
| Original simultaneous-vs-alone RII (baseline competitive effect) | 0.18483 | **0.17554** | **5.0%** |

Source-grain MSE gain for the priority-only RII: **+0.01653**, species bootstrap 95% **[−0.00088,+0.04003]**, **10/15** species improved. Under the 199-shuffle GR mapping placebo, observed species-equal MSE **0.04685**, random-null median **0.06271**, at least as good fraction **0.005**.

For simultaneous competitive RII, MSE gain **+0.00344**, species-bootstrap 95% **[−0.00116,+0.00872]**, **7/15** species improved; same mapping fraction **0.005**. Thus **growth potential predicts ordinary competition too**. A larger relative improvement for the priority-only outcome is suggestive but **not a preplanned or independently replicated quantitative test of priority-specificity**.

**What survives:** earlier-neighbor solo growth performance is a better candidate predictor of unseen partner influence than the six measured static morphological traits in this panel; this does not isolate nutrient preemption from spatial competition, root/canopy access, or habitat modification. Growth rate is itself measured from **solo control plants after ten weeks**, not a universally available ex-ante trait.

## A mechanism hypothesis that is now precisely testable

**Growth-opportunity-gated interaction hypothesis:** With the focal age/density held constant, increasing the *actual pre-arrival biomass or resource-capture rate* of neighbor j increases its age-matched net suppression of focal i, beyond focal trait tolerance and simple guild identity.

Test in a new factorial experiment with (i) focal species, (ii) neighbor identity, (iii) **pre-arrival neighbor biomass matched vs allowed to differ**, (iv) neighbor arrival timing, (v) neighbor-present vs age-matched focal-alone controls, and (vi) root-contact or resource-access barriers. Measure repeated absolute biomass and resource uptake. Only a growth-potential × early-neighbor-access intervention can distinguish *biotic interaction formation* from an intrinsic growth advantage and background resource depletion.

Finally test whether the **experimentally measured** history-dependent interaction matrix predicts the **same experimental communities'** named-vs-identity-free abundance structure. Currently the pairwise interaction study and the four KNK multispecies-history studies are separate experiments: no mediation path has been identified, no random-effects meta-analysis should combine their raw pair outcomes.

## Provenance and claim boundary

- Original source and in-sample trait claim: https://doi.org/10.1002/ecy.4528; source files https://doi.org/10.5281/zenodo.10512662
- Morphology-only strict LOSO protocol: `docs/VANSTEIJN_UNSEEN_NEIGHBOR_TRANSFER_PRE_SCORE_20261124.md`.
- Post hoc role-specific GR protocol: `docs/VANSTEIJN_NEIGHBOR_GR_POSTHOC_TRANSFER_CONTRACT_20261124.md`.
- Endpoint triangulation protocol: `docs/VANSTEIJN_GR_ENDPOINT_TRIANGULATION_PRE_SCORE_20261124.md`.
- Executables: `benchmarks/predict_vansteijn_unseen_neighbor_traits.py` and `benchmarks/decompose_vansteijn_neighbor_growth_rate.py`, with synthetic no-leakage tests in `tests/`.
- Aggregates: `validation/interaction_formation/vansteijn_strict_unseen_neighbor_transfer_v1.json`, `vansteijn_neighbor_growth_rate_transfer_v1.json`.
- Source-authentic GitHub Actions runs:
  - https://github.com/zuizui0223/knk/actions/runs/45205233811
  - https://github.com/zuizui0223/knk/actions/runs/45206403995
  - https://github.com/zuizui0223/knk/actions/runs/45207378978
- Core prior result remains the already published van Steijn trait–priority-effect association; this is a **small-panel transfer evaluation**, not a novel experimental discovery or direct competition coefficient. No changes to frozen v28–v38 or present Ecology Letters manuscript.
