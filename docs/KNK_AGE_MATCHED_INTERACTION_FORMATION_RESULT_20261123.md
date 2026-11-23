# Assembly history changes realized biotic-neighbor effects — KNK direct-source pilot

**2026-11-23 | Draft interaction-formation line only | original KNK v28–v38 frozen**

## Verified experimental age structure

**Essential correction** to earlier preliminary interpretation: van Steijn et al. (2025), *Ecology* (DOI https://doi.org/10.1002/ecy.4528) deliberately provided **10 weeks of growth for the focal plant** in each condition. Simultaneous focal+neighbor and focal-alone controls were planted in week 3 and harvested in week 13. The sequential treatment introduced the early neighbor in week 1, the *focal* in week 3, and harvested both in week 13. The early neighbor gets **13 weeks**, but the focal always gets **10 weeks**. Methods and Figure 1 independently confirmed via https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecy.4528 and https://pub.epsilon.slu.se/36436/1/van-steijn-t-l-h-et-al-20250207.pdf.

This means the assay can estimate a **history-specific net effect of neighbor presence on focal biomass at a fixed focal age**. However, the experimentally bundled mechanism of early neighbor establishment (size, root/canopy preemption, nutrient consumption, microbial soil conditioning) is not separated. See method erratum `docs/VANSTEIJN_FOCAL_AGE_METHOD_CORRECTION_AND_NET_NEIGHBOR_EFFECT_20261123.md`.

## Source identity and matched experimental units

Public source Zenodo v2 https://doi.org/10.5281/zenodo.10512662. 15 source-verified taxa (5 grasses/5 nonleguminous forbs/5 legumes), 74 age-matched solo control pots, 206 simultaneously grown directed-pair focal observations, 206 sequentially grown directed-pair observations. There are **202 exact matched focal × neighbor directed combinations**; 4 only in simultaneous and 4 only in sequential raw tables. **One pot per directed species pairing**; do not treat 202 rows as 202 independent species replicates.

R original code defines `RII_sim = (S-A)/(S+A)` and `RII_seq = (L-S)/(L+S)`, where \(S\) is simultaneous focal biomass, \(A\) is mean solo focal biomass and \(L\) is late-focal biomass after a 3-week-preestablished neighbor.

## A history-modified age-matched net neighbor response

Because solo, simultaneous and later-arriving focal plants each grew 10 weeks, the following is a meaningful total-priority intervention contrast:

\[
C_{ij}^{sim} = \frac{S_{ij}-A_i}{S_{ij}+A_i},\quad
C_{ij}^{pre} = \frac{L_{ij}-A_i}{L_{ij}+A_i},\quad
\Delta C_{ij}=C_{ij}^{pre}-C_{ij}^{sim}.
\]

Negative \(C\): mixed focal had lower biomass than the alone control. \(\Delta C<0\): a neighbor with a three-week establishment advantage imposed a stronger *net observed cost* on focal plant biomass than the same neighbor entering simultaneously. This is **not** an estimate of the time-varying per-capita competition coefficient \(\alpha_{ij}(t,H)\). The source already demonstrated priority effects with a different normalization. Our contrast is a newly derived description using source controls, not a new manipulative experiment.

### Quantitative results

| Endpoint | Value |
|---|---:|
| Mean age-matched net neighbor RII, simultaneous | **−0.2591** |
| Mean net neighbor RII, neighbor pre-established 3 weeks | **−0.7864** |
| Mean \(\Delta C\), with 15 focal species equally weighted | **−0.5289** |
| 20,000-bootstrap 95% interval over **15 focal species**, \(\Delta C\) | **[−0.6059, −0.4444]** |
| Focal species whose mean \(\Delta C\) is negative | **15/15** |
| 15-species exact sign-flip (symmetry diagnostic) | p=0.000061 |
| Directed pairs in which simultaneous net effect is positive but pre-established effect negative | **22/202** |
| Directed pairs with negative net neighbor effect when pre-established | **202/202** |

**Important:** 22 sign changes are unreplicated single-pair observations. Some simultaneously positive outcomes could arise by sampling variation, and these sign changes are **not** validated mutualism-to-competition conversions.

### Who comes early modulates the net effect: species identity versus guild

- After subtracting each focal species' mean \(\Delta C\), the early-neighbor species identity explains **30.56%** of the remaining schedule-effect variation. Within-focal neighbor-label shuffle null median **7.17%**, p=0.0001 (9,999).
- After adjusting both the focal species and early-neighbor grass/forb/legume guild, named early-neighbor identity still tracks **25.19%** of the residual variation. Conditional within-focal×guild permutation-null median **6.11%**, p=0.0001.
- Guild category itself accounts for 7.36% of focal-centered variation, but actual 5/5/5 mapping is **not exceptional** among arbitrary balanced guild allocations (placebo upper fraction 18.2%). The test is a mapping specificity check, not evidence that true functional guild is biologically irrelevant.
- In leave-one-*unordered*-species-pair-out prediction with fixed L2 ridge (both orientations withheld together), \(\Delta C\) RMSE improves **0.2223 → 0.1997** when adding early-neighbor identity to focal identity (relative 10.2% reduction). This is within-source pair prediction; no external taxon validation.

### Extreme named partners: robust and nonrobust conclusions

Across 15 earlier neighbor choices, with a familywise **maxT mapping** guard (9,999 within-focal label shuffles):
- *Antennaria dioica*: adjusted mean \(\Delta C\) is **+0.3093** (less harmful than focal-average history cost), **maxT p=0.0001**.
- *Briza media*: adjusted mean +0.1326, maxT p=0.2351 — **not exceptional after partner search correction**.
- *Trifolium repens*: adjusted mean −0.07545, maxT p=0.9479 — **not exceptional after partner search correction**.

Some of the more extreme \(RII_{seq}\) partner results in the earlier exploratory calculation do **not** translate into equally extreme net neighbor changes relative to age-matched solo controls. The latter endpoint is more relevant to interaction formation and should govern the interpretation.

### Above vs below

The previously audited late-arrival RII difference between above and below was −0.00476, 15-focal bootstrap interval [−0.02786,+0.02026], sign-flip p=.712, so there is no broad support for exclusively aboveground or belowground late-schedule costs. The new score independently computes organ-specific age-matched \(\Delta C\), but the primary response remains whole-plant biomass.

## What this genuinely adds to KNK

The **population outcome mechanism** is now more directly connected to experimental evidence:

```
Early-neighbor identity × establishment head start
      ↓
age-matched net impact on the focal species (ΔC)
      ↓ hypothesized, NOT directly tested in the same dataset
which species becomes dominant / rank-shape vs identity-memory of an assembled community
```

The earlier KNK experiments (Leopold fungi, Alonso rhizobox, Dawson-Glass, POEM field) measured assembly history versus group abundance architecture, but did not measure the same full directed pair \(\Delta C_{ij}\) network. Hence **the arrow from an interaction network to a community's storage signature remains untested**. Do not pool the 202 directed pair values with four independent multispecies experimental systems as if all were the same outcome or independent studies.

The existing literature preempts several novelty claims:
- Original van Steijn et al. (2025) already concluded that **focal and neighbor competitive traits** explain priority-effect strength.
- Fragata et al. (2022) already showed sequence-dependent competition and predicted coexistence through niche pre-emption in spider mites (DOI https://doi.org/10.1111/ele.14021).
- Powell-Romero, Wells & Clark (2024), *Ecology Letters* DOI https://doi.org/10.1111/ele.14509 already showed that ignoring order of arrival makes asymmetric biotic interactions difficult to infer from co-occurrence.

New KNK claim should therefore be narrow: **the age-matched realized effect of an established plant neighbor varies strongly across directed partner species, and named neighbor identity retains predictive structure beyond coarse grass/forb/legume grouping, with strong source-aware selection controls.**

## Next decisive experimental manipulation

In an independently repeated plant panel, cross **neighbor identity × history of establishment × neighbor presence × spatial root-contact barrier**, with constant **focal age and density** in all cells. At 3+ time points measure absolute plant mass, soil moisture/mineral N, and root/canopy occupancies. Fit \(g_i(t;H,j\text{ present})-g_i(t;H,j\text{ absent})\) under matched conditions and compare it with final composition and its rank-sorted counterpart in the *same experimental community*. Only then can the full chain *formation of biotic interactions → retention of history in community architecture* be causally evaluated.

## Provenance / current status

- Frozen pre-score contract: `docs/VANSTEIJN_PAIR_INTERACTION_PRE_SCORE_20261123.md`; transparent age correction `docs/VANSTEIJN_FOCAL_AGE_METHOD_CORRECTION_AND_NET_NEIGHBOR_EFFECT_20261123.md`
- Source-locked scoring: `benchmarks/score_vansteijn_partner_effects.py`, `benchmarks/audit_vansteijn_partner_specificity.py`, `benchmarks/score_vansteijn_age_matched_neighbor_effect.py`
- Exact aggregate-only result `validation/interaction_formation/vansteijn_age_matched_net_neighbor_result_v1.json`
- Authentic public-source Actions execution https://github.com/zuizui0223/knk/actions/runs/45096126914; source MD5 pass and tests success.
- The original KNK Ecology Letters manuscript, frozen versions v28–v38, and Dawson-Glass/POEM Draft PRs remain unchanged. This work is a **separate unmerged Draft PR**.
