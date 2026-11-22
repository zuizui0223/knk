# Dawson-Glass independent validation: selection-corrected competition audit (2026-11-22)

**Status: exploratory, same 24 experimental pots; no new independent replication.** This document refines, and in places weakens, the interpretation of the earlier within-treatment *Rumex*–grass association. **Do not alter KNK v28–v38 or the frozen Ecology Letters manuscript.**

## Exact source, grain and provenance
- Dawson-Glass et al. (2025), *Ecology*, DOI https://doi.org/10.1002/ecy.4504; Dryad https://doi.org/10.5061/dryad.sbcc2frh4.
- User-supplied original `Dawson-Glass_ECY24-0464_assembly_data.csv`, SHA-256 `1514d4501a0bdf7a61ace6c250566c5a0b8e4a24a0f4737646d182ebad7b8230`.
- 360 source observations = 24 pots × 15 named species, 12 in each order; 21 recorded zero-biomass entries. Raw source remains **outside GitHub**.
- Reproducible code in `benchmarks/audit_dawsonglass_competition_caveats.py`; aggregate results in `validation/dawsonglass_external_postsubmission/selection_caveat_audit_v1.json`.

## 1. Dominant species drive most, but not all, turnover
The mean between-arm species-labeled **relative biomass** profiles are separated by Bray–Curtis distance **0.3841**. We decomposed its numerator into 15 absolute species-wise mean differences.

| Species | Difference in average biomass share, warm−ambient | Fraction of summed absolute species changes |
|---|---:|---:|
| *Rumex crispus* | +0.3169 | 41.25% |
| *Dactylis glomerata* | −0.1689 | 21.99% |
| *Bromus inermis* | −0.1240 | 16.13% |
| Remaining 12 taxa | — | **20.63%** |

The three dominant-shifting species therefore account for **79.37% of this particular L1 turnover summary**. This is a descriptive accounting identity of two arm means, not a partition of mutual information, competitive impact, or historical-retention E.

## 2. A post hoc search correction weakens the evidence for specific competition
Previously, *Rumex* absolute grams versus *Bromus + Dactylis* absolute grams, after demeaning within the randomized order treatment, gave **r=−0.5306**. The same sign appears separately in ambient **r=−0.5542** and warm-order **r=−0.5703** mesocosms (12 per arm).

The two grasses were selected *after* inspecting biomass responses. We therefore treated **all 91 pairs of the other 14 species** as possible comparison groups and permuted the complete `Rumex` residual vector *within treatment arms*, holding the observed covariance among the 14 comparison taxa fixed. 49,999 permutations; seed 20261124.

- Selected-pair unadjusted **one-sided p=0.00368**.
- MaxT minimum-correlation **familywise one-sided p=0.0669** for selecting the strongest negative pair among 91.

**Caveat:** This corrects selection of the competitor pair **conditional on Rumex already being nominated**. Rumex itself was identified after seeing the outcome, so this is *not* a fully corrected discovery p-value. Evidence for a specific biologically targeted antagonistic pair is thus weaker than the former unadjusted p≈0.0058 implied. No direct competitive interaction has been established.

## 3. Why an absolute-grams negative correlation is still nonmechanistic
Set `R = Rumex grams`, `G = two grass species grams`, `O = the other 12 species grams`, `T=R+G+O`. After removing arm means, **exactly**:

```
Cov(R,G) = Cov(R,T) - Var(R) - Cov(R,O)
-32.51855 = -7.69616 - 26.46920 - (-1.64681)
```

This is an **algebraic identity**, not a biological model. It demonstrates why negative covariance in absolute abundance does not uniquely identify competition: shared total-pot constraints, growth-time differences, and unmeasured environmental factors can induce the same observation. Whether the pots were actually capacity-limited is *unmeasured*, not demonstrated by this formula.

## 4. History exists beyond the three biggest species
When the three dominant-shifting species were removed, the remaining 12 species still distinguished arrival-order regimes:

| Post hoc representation | Null-calibrated E | Permutation p |
|---|---:|---:|
| Remaining 12 taxa, raw grams | +0.2390 | 0.0001 |
| Remaining 12 taxa, reclosed relative grams | +0.3159 | 0.0001 |
| Remaining 12 taxa, relative values sorted without species labels | +0.0572 | 0.1233 |

Both raw grams and relative composition of the *remaining* 12 species retain an order-associated fingerprint; this need not be driven only by renormalization. But these are **new representations chosen after seeing the original result**, and the null is reused on the same 24 pots. They are not separate replications, nor do their E values add up to the complete-community E.

Exploratory individual-species raw-gram associations after 15-species BH correction include *Solidago canadensis* (−3.90 g per pot, q≈0.00075), *Ranunculus pennsylvanicus* (−0.788 g, q≈0.027) and *Hypericum punctatum* (−0.185 g, q≈0.038), in addition to *Rumex* and *Dactylis*. These are candidates, not new causal claims.

## 5. Initial-size confounding is not solved by regression
The source records *Rumex* height just before planting: ambient mean **0.783 cm** versus warm-order mean **2.208 cm**. The two arms' observed height ranges overlap only from **0.8 to 1.4 cm**: **4 ambient and 2 warm-order pots** fall within this common-support range. The height variable may be affected by the germination/timing schedule, so adjusting for it conditions on a post-intervention state with very poor overlap; adjusted coefficients cannot identify a competition effect.

## Decision and next discriminating experiment
- **Retain:** arrival order changes species identity/composition far more clearly than unlabeled rank architecture in this experiment; top-three winner switching is large, and subordinate taxa also carry a history-associated compositional signature.
- **Downgrade:** specifically antagonistic *Rumex*–grass competition, true biomass compensation mechanisms, and equality of total biomass are **not identified**.
- **Next intervention:** independently randomize arrival schedule and initial plant size, cross with partner presence/absence, include single-species and simultaneous-arrival controls, measure species-specific **absolute** biomass and soil-resource trajectories over time. A competitor-specific inhibition term beyond the single-species age/size effects would separate interaction from an age-only mechanism.
- **Boundary:** Do not append this exploratory dataset to the pre-submission EL letter; keep Draft PR #3 unmerged until source-method and author review.

