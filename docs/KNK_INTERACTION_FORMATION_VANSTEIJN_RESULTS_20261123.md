# KNK: from assembly-history memory to biotic interaction formation
*2026-11-23 • mechanism track, post-submission development only*

## Evidence grade

**Strong new-to-KNK quantitative reanalysis**, NOT a new original experimental demonstration. The original papers established both order-dependent competitive outcomes (Fragata et al. 2022) and trait-dependent partner effects (van Steijn et al. 2025). KNK adds an explicitly audited **directed interaction-response matrix**, independent species-vs-guild label controls, uncertainty and selection checks, and the bridge to the rank-vs-label storage problem. **It has not identified a dynamic interspecific interaction coefficient.**

### Separate what has and has not been measured
- **Observed experimental schedule response**: species-pair-specific final biomass when the focal is planted **second** compared with the same species pairing planted **simultaneously**.
- **Potential mechanistic quantity**: history-dependent direct neighbor effect on age-/initial-size-matched focal growth, `I_{j→i}(H)=g_i(j present,H)-g_i(j absent,H)`, which **requires additional controls absent from the current raw files**.
- **Community-storage outcome** (previous KNK): named species composition versus the same species abundances sorted with identities discarded. The direct-pair dataset does not measure multispecies architecture on an experimental unit and **cannot be pooled** with that outcome.

**Source-method correction (2026-11-23):** the original Figure 1 and Methods confirm 10-week focal growth in solo, simultaneous and sequential treatments. Earlier wording about a focal-age confound was incorrect; neighbor pre-establishment is the actual manipulated difference. This strengthens the identification of a *realized age-matched net neighbor effect*, while pre-emption versus environmental conditioning remain unresolved. See [erratum](VANSTEIJN_FOCAL_AGE_METHOD_CORRECTION_AND_NET_NEIGHBOR_EFFECT_20261123.md).

## Source provenance and reconstruction

van Steijn et al. (2025), *Ecology* 106 e4528, DOI https://doi.org/10.1002/ecy.4528; source Zenodo **v2** https://doi.org/10.5281/zenodo.10512662 (5 files, original MD5s checked). 15 plant species: five grass, five forb, five legume. 74 solo-control pots (~5 per species; one species n=4), 206 focal–neighbor directed simultaneous responses and 206 later-planted responses. **Only 202 directed pair keys match one-to-one**; four keys occur only in each treatment file. Each directed-pair combination normally has one pot, so there are **15 focal clusters, not 202 independent species replicates**.

From the exact original R analysis (lines 60–104), compute:
```
RII_sim(i,j)  = (Biomass_sim(i,j) − mean(Biomass_alone(i))) /
                (Biomass_sim(i,j) + mean(Biomass_alone(i)))
RII_late(i,j) = (Biomass_late(i,j) − Biomass_sim(i,j)) /
                (Biomass_late(i,j) + Biomass_sim(i,j))
```
Each also computed separately for roots and shoots. Note `RII_late` is **not** an experimentally isolated causal competition coefficient: the neighbor has a 3-week establishment advantage while the focal's 10-week age is matched; this is a total pre-establishment effect, not a mechanistic decomposition. The paper already analyzed RII responses and traits; KNK's reanalysis tests the structure of a directed 202-edge schedule-response matrix, subject to one pot per directed combination.

## Confirmed quantitative results

### Global response and anatomy
- Late planting produces lower focal total biomass in **200 of 202** matched directed pairs; 2 show modest positive changes (not sufficient evidence of facilitation).
- Across the **15 focal species means**, average `RII_late,total=−0.7309`, species-level bootstrap 95% `[−0.8144,−0.6345]`; all 15 focal means negative. Two-sided focal-species sign-flip diagnostic p=0.000061. This is primarily confirmation of the source's known arrival-time disadvantage, not a new mechanism.
- Shoot-vs-root RII_late difference across 15 focal species `−0.00476`, 95% bootstrap `[−0.02786,+0.02026]`; sign-flip p=0.712. No global organ-specific contrast supported.

### Who arrives first matters beyond which functional group
- After subtracting each focal species' average late response, the **same early-neighbor name consistently explains 47.74%** of residual pair-to-pair variation (`R²_neighbor|focal`); conditional within-focal label-null median **7.20%**, p=0.0001.
- After adjusting for both focal species and early-neighbor **functional guild** (grass/forb/legume), a species-specific partner label still tracks **41.15%** of the residual variation, within-focal×guild shuffled-label null median **5.79%**, p=0.0001.
- By contrast the observed guild grouping captures **11.51%** of focal-centered residual variation, and its rank among random balanced **5/5/5** species-to-guild partitions is not unusual (upper fraction **16.9%**). This is a *mapping specificity* observation, not an experimental randomized null of guild importance.
- Removing the two most positive partner identities (*Antennaria dioica* and *Briza media*) leaves partner-identity residual R² **20.43%**, label-permutation p=0.0005 (1,999 permutations). The identity result is not explained solely by those two.
- In leave-one-*unordered*-pair-out predictions (both directions held out at once, fixed ridge penalty 4), model with focal identity only predicts `RII_late` at RMSE **0.24335**; adding early neighbor species identity improves to **0.19145** (≈21.3% lower RMSE). This is within-study held-out *pair* prediction, not external species generalization.

### Two extreme partner types, corrected for discovering the extremes
Effects below are **deviations from the focal plant's own average late-response penalty**, not measured benefits or intrinsic `alpha_ij`:

| Earlier neighbor | Focal-adjusted RII shift | Direction consistency | 15-partner maxT mapping-FWER p |
|---|---:|---:|---:|
| *Antennaria dioica* | **+0.455** | lower cost in 13/13 focal species | **0.0001** |
| *Briza media* | +0.240 | lower cost in 11/13 species | 0.0044 |
| *Trifolium repens* | **−0.187** | stronger cost in 14/14 species | 0.0456 |

maxT controls **post hoc choice among 15 named neighbors** under a within-focal neighbor-label exchangeability null. It is NOT a randomized-treatment causal p-value.

### Interaction-direction asymmetry is common, but sign reversal rare
99 unordered species pairs have both directed orientations. Mean absolute difference between the two `RII_late` responses is **0.3875**, correlation between orientations **−0.230**, and only **2/99** have opposite `RII_late` signs (because nearly all late responses are negative). Magnitude asymmetry must not be mislabeled as a switch from antagonism to facilitation.

## Biological synthesis across KNK's research lines

The combined body of evidence suggests a **testable hierarchy**, not a proven universal mechanism:

1. **Species-specific neighbor filters**: a given earlier colonist can impose a repeatable schedule-dependent constraint (van Steijn, Fragata).
2. **Identity-rich storage**: those filters can alter *which species wins* without much changing relative rank shape (Dawson-Glass; some POEM periods).
3. **Architecture-rich storage**: conditions altering strength of dominance or access to preferred spatial resources can instead change unlabeled abundance inequality (Populus fungi; other POEM periods).
4. **Observable function**: total biomass, rust, and root/space outcomes may or may not reflect these past changes. Neither correlation nor the current pilot demonstrates causal mediation from interaction formation to storage.

This is compatible with source-reviewed priority-effect theory. The 2024 study Powell-Romero, Wells & Clark (*Ecology Letters*, DOI https://doi.org/10.1111/ele.14509) already established by simulation that inferring asymmetric species interactions from cross-sectional co-occurrence without arrival histories can fail. KNK must explicitly acknowledge this precedent, rather than claim the insight as novel.

## Next discriminating experiment

Choose four focal species from contrasting functional guilds and two candidate neighbors from the opposite extremes (**Antennaria** versus **Trifolium repens**, after checking no self-pair). Keep **focal age and initial density fixed** at evaluation. Randomize three neighbor-installation schedules (before, same time, after focal introduction) and collect matched **partner-absent focal controls** at corresponding focal measurement times:
- Four focal taxa × two neighbors × three schedules × six replicate pots = **144 mixture pots**
- Four focal taxa × three matched schedules × six replicate pots = **72 focal-only control pots**
- **216 pot units**, not counting optional spatial-barrier and environmental-conditioning contrasts.
- Predefine response `ΔI_{j→i}` based on **focal per-capita absolute biomass growth**, ideally repeatedly observed, rather than relative final biomass. Measure partner biomass, canopy/root space occupation, soil nitrogen and moisture to distinguish age advantage, resource pre-emption, environmental conditioning, and competition.
- Test out-of-sample across entirely withheld **neighbor species** (not merely withheld pots involving the same two focal neighbor species). With only two extreme neighbors, this latter generality gate **cannot** be satisfied; independent additional partner species are needed before making a universal functional-trait law.
- A resource access or partner-removal intervention after establishment is required to test mediation; otherwise a chronological treatment difference remains ambiguous.

## Reproducibility and stopping rule
- Immutable public original MD5s captured in `benchmarks/probe_vansteijn_competition.py`
- Primary response matrix: `benchmarks/score_vansteijn_partner_effects.py`
- Species-vs-guild specificity, maxT and influential-partner deletion: `benchmarks/audit_vansteijn_partner_specificity.py`
- Sources and all aggregate JSON in `validation/interaction_formation/`
- Source-audited GitHub Actions runs [primary](https://github.com/zuizui0223/knk/actions/runs/45095129265) and [selection-corrected follow-up](https://github.com/zuizui0223/knk/actions/runs/45095655693), both successful.

**No modification of existing KNK main manuscript, frozen v28–v38, or other Draft PRs. Do not claim this observational decomposition alone identifies a formed causal interaction network.**
