# KNK interaction formation: source-verified mechanistic pilot (2026-11-23)

**Status: independent public two-species mechanistic anchor, already-established original discovery, not a new KNK mechanistic causal identification.** This stays **unmerged** and wholly separate from frozen v28–v38 (Ecology Letters), Dawson-Glass Draft PR #3 and macro/meta Draft PR #4.

## Biological question
Not just *where* does assembly history appear in final abundance composition, but *how* does arrival history alter the **realized per-capita effect of another species**, its spatial context, and therefore possible coexistence?

In notation, for species `i` and `j` and assembly history `H`, an experimental interaction estimand is

```
I_(j→i)(H) =
  [log(N_i,end / N_i,start) with j present, H]
− [log(N_i,end / N_i,start) with j absent, age- and density-matched H]

ΔI_(j→i) = I_(j→i)(j-first) − I_(j→i)(i-first).
```

Such an estimand **requires the matched partner-absence controls**, so the current two-species abundance-ratio comparison alone is not an estimate of `I`, `alpha_ij`, competition network adjacency, or long-term stability.

## Independent source
Fragata et al. (2022) **Specific sequence of arrival promotes coexistence via spatial niche pre-emption by the weak competitor**, *Ecology Letters* 25, 1629–1639, DOI https://doi.org/10.1111/ele.14021. Original raw files and analysis code at https://github.com/irfragata/order_arrival_niche_preemption, immutable commit `f1568dfb4a0da7bc9811985e416daec86447d0f3`.

Original study **already established** order-dependent competitive outcomes and, using additional density contrasts plus a Beverton–Holt model, predicted coexistence if weaker *Tetranychus urticae* (Tu) arrives before stronger *T. evansi* (Te); leaf occupancy shifts support spatial niche pre-emption. Do not present that published finding as novel.

### Raw source and unit audit
90 **whole box** units, nine treatments (n=10 each, five per block × two chronological blocks), 3,214 nested leaf/leaflet × leaf-surface observations. Git object checks for corrected raw CSV `dd734646364eb75464ab16b0169afb4b06b6e6d4` and leaf age map `354dd41ece632b8e5f8bd0266233f5a88f1517e5` passed. Original source R code explicitly maps `10Te10Tu48h` to **Tu-first** and `10Te48h10Tu` to **Te-first**.

Only the initial-10:10 balanced competition cases were included in the principal reanalysis, with ten whole boxes per schedule.

## Recomputed exact box-level effects

| Equal initial numbers | Mean final Te fraction | Mean Te adult females | Mean Tu adult females |
|---|---:|---:|---:|
| Te arrived 48 h before Tu | **0.8106** | 354.0 | 74.5 |
| Arrived simultaneously | 0.5957 | 236.8 | 149.2 |
| Tu arrived 48 h before Te | 0.4104 | 153.1 | 204.7 |

Difference (Te-first minus Tu-first), final Te fraction **+0.40021** = **+40.0 percentage points**. Because the balanced pair contains 5 boxes from each of two blocks, **all 252²=63,504** valid within-block 5-vs-5 assignments were enumerated; exact two-sided p=**0.00006299**. Independent box bootstrap stratified by block×order (20,000 resamples) gives a descriptive percentile **95% interval +0.2991 to +0.4921**.

This is consistent with published order-dependent competitive outcomes **but does not identify which interspecific coefficient changed**, distinguish growth-time head starts from suppression, or establish long-term coexistence from final abundances alone.

### Negative / interpretation-sensitive check

Per whole box, the fraction of Te adults on leaf 4 averaged:
- Te-first: **0.3796**
- simultaneous: **0.4614**
- Tu-first: **0.3398**.

In this **simple 10:10, box-level, unadjusted leaf-4 fraction** diagnostic, Te-first minus Tu-first = **+0.03984**, exact two-sided p=**0.6140**. Tu-first minus simultaneous = −0.1216, p=0.1323. None survives three-pair Holm correction. Thus the current reduced pilot does **not** independently identify the reported spatial mechanism. The source authors used more complete leaf, location and initial-frequency analyses and coexistence modelling. Do not claim that a large overall abundance effect proves a resource-occupancy shift.

## How to test newly formed interactions instead of just final abundance

**Mechanistic prediction A — competition changes with history.** Manipulate species identity, arrival order, initial density and species-partner absence. Compute `I_(j→i)` as partner-dependent change in absolute per-capita growth at matched time exposure. A statistically and biologically meaningful `ΔI` would show order dependence of realized interaction strength.

**Mechanistic prediction B — spatial interaction gates change.** Repeat the pair factorial with compartment access/leaf strata barriers or resource patch pre-emption, and measure absolute location-specific population densities. If a resource-access intervention removes the order × partner effect, it supports spatial mediation beyond the end-state correlation.

**Mechanistic prediction C — guild/trait-dependent generalization.** In communities with 3+ interacting species, independently measure some partner-removal interactions, then ask whether history-dependent changes in interaction strengths predict name-based vs identity-free rank retention. Avoid simply inferring a full interaction network from covariances with n=24 pots (Dawson-Glass) or 20 POEM plots per cohort.

Mechanism-ready plant example: Yan et al. (2024), DOI https://doi.org/10.1186/s13717-024-00483-y (arrival order × root isolation × partner removal). Another design: van Steijn et al. (2025), https://doi.org/10.1002/ecy.4528 (pairwise species arrival and competitive response). Their published mechanisms are background; they require separate raw-data audits before new KNK comparisons.

### Meta-analysis taxonomy

Separate `direct manipulated neighbor effect` (interaction-strength lane), `spatial resource-use` (mediation lane), `final community labeled/rank state` (storage lane), and `long-run coexistence` (demographic lane). A paper/experiment may contribute correlated measurements to more than one lane, but **not independent meta-study counts**. No cross-system pooled interaction estimate exists yet.

## Reproduction and code status

- `benchmarks/audit_fragata_interaction_order.py` — source-locked whole-box reconstruction, exact blocked permutation, 20,000 within-cell whole-box bootstrap.
- `tests/test_fragata_interaction_order.py` — synthetic design and null guards.
- `validation/interaction_formation/fragata_2022_exact_box_summary_v1.json` — aggregate-only scored values, no raw data redistribution.
- GitHub Actions source/score/test run: https://github.com/zuizui0223/knk/actions/runs/45094006190 — **success**.
- Frozen KNK manuscript, v28–v38, and all existing experiment-specific figures remain **unchanged**.

**Verdict:** Historical order can alter the outcome of realized competition; explaining the *formation of interaction strength* in the more complex KNK grassland/fungal systems requires stronger experiments than relative-abundance reanalysis.
