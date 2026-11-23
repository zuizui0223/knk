# KNK interaction formation — independent mechanistic anchor (2026-11-23)

**Track:** separate mechanism-development branch. **No changes** to frozen KNK v28–v38 or ongoing Draft PR #3 (Dawson-Glass) / #4 (macro-meta). This is not a novel discovery of mite coexistence; it is a grounded experimental contrast and a mechanism-audit anchor for a future cross-system programme.

## Core conceptual distinction

Historical state storage (KNK): `assembly schedule H -> present community structure C`.

Interaction formation (next work): `H -> realized interaction matrix A(H,t,space) -> abundance, coexistence, function`. The matrix can depend on species identities, arrival times, density, spatial resource use, and environmental/host feedback.

A difference in final species abundances **does not identify** the interaction matrix or prove competitive suppression. To estimate `alpha_ij(H)`, one needs species alone controls, density/frequency contrasts, repeated absolute demographic measurements and/or manipulated spatial overlap. Avoid inferring network edges from correlation in a 15-species, 24-pot cross-section (Dawson-Glass) or repeated POEM field plots without partner-removal controls.

## External positive mechanism result already known in the source study

Fragata et al. (2022), *Ecology Letters* 25:1629–1639, DOI https://doi.org/10.1111/ele.14021, experimentally varied order and starting frequency of two spider-mite species, *Tetranychus evansi* (generally superior competitor, Te) and *T. urticae* (Tu).

Original paper's key conclusion (not KNK novelty): Tu-first arrival could enable predicted coexistence by spatial niche pre-emption of preferred tomato-leaf locations, reducing asymmetry between competitors; the same order effect was **not** an example of positive frequency-dependent priority effects. Original source also measured plant-defense gene expression, suggesting this outcome was not explained by differential induction of plant defenses.

GitHub source: https://github.com/irfragata/order_arrival_niche_preemption ; pinned commit `f1568dfb4a0da7bc9811985e416daec86447d0f3`, public corrected data blob `dd734646364eb75464ab16b0169afb4b06b6e6d4`, leaf-age map blob `354dd41ece632b8e5f8bd0266233f5a88f1517e5`. Verify Git object identity before analysis; raw CSV is only temporarily downloaded by CI, never committed to KNK.

## Independent KNK quantitative audit — proposed before CI scoring

**Sample unit: entire experimental box** (not leaves/leaflets/upper-lower leaf sides). Nine manipulation categories, two chronological blocks, 5 independent boxes per category per block (90 boxes total). Corrected source has 3,214 within-box leaf/leaflet count records.

Primary targeted contrast restricts to **matched initial frequency 10 Te + 10 Tu**, three order conditions, 10 boxes each:
- `10Te10Tu` = simultaneous;
- `10Te10Tu48h` = Tu-first (mapping verified directly from source analysis R code);
- `10Te48h10Tu` = Te-first.

Primary descriptive outcome: proportion of adult females Te in each box after two generations, `pTe=Te/(Te+Tu)`. Pairwise contrasts use **exact block-stratified permutations**: for two order treatments, 5 vs 5 boxes in each block, `choose(10,5)^2=63,504` unique assignments; two-sided p by absolute mean difference. Bootstrap 20,000 whole boxes resampled within each treatment×block for a 95% conditional box-level interval.

Secondary post hoc: Te and Tu absolute abundance, the per-box fraction of Te occupying leaf 4, and leaf spatial-overlap index. These describe a spatial **association**, not an intervention on niche access or a direct estimate of competition coefficients. The published paper used further initial frequency conditions and Beverton–Holt modelling to infer niche and fitness differences; do not reproduce that mechanistic claim from a single aggregate contrast.

### Competing hypotheses for KNK mechanism synthesis

| Prediction | Competitive priority | Spatial niche pre-emption/equalisation | Growth duration alone |
|---|---|---|---|
| Abundance changes with order | yes | yes | yes |
| Spatial occupancy / resource distribution shifts with order | possible | predicted | not necessarily |
| Extra effect beyond single-species age-matched growth | predicted | predicted | no |
| Competitor-removal or barrier rescues suppression | predicted | predicted if relevant barrier | no |
| Source dynamics yield stable coexistence | not guaranteed | possible with reduced fitness asymmetry | not guaranteed |

For plant follow-up, source designs explicitly manipulating root barriers and partner removal (e.g., Yan et al. 2024, https://doi.org/10.1186/s13717-024-00483-y) offer more direct leverage than pairwise correlations in Dawson-Glass. The source published this design and its positive biomass/biodiversity and root-interaction findings; these are not new KNK findings.

## What would count as novel cross-system evidence

Predefine interaction contrasts such as `I_ij(H) = log[Y_i(mixed, H)/Y_i(alone, matched age and density)]`, where positive suggests net facilitation/negative suggests suppression only under proper matched controls; or causal growth response with/without partner, controlling exposure duration. Compare **within pair** the shift in interaction strength when arrival order reverses, `ΔI_ij`. Then link experimentally estimated `ΔI` to history-storage location (identity-rich vs architecture-rich) and coexistence predictions on **independent source studies**, with publication-cluster dependence accounted for.

Do not treat a spider-mite two-species pair as an independent multispecies rank-architecture replication or mix these effect sizes in the same meta-analysis without a clearly different lane.

## Boundaries

- Prioritize a mechanism-focused follow-up **separate** from the Ecology Letters submitted-candidate letter; the old KNK mechanism ceiling was explicit.
- Leaf-location aggregates alone do not identify root/resource competition or niche modification; no assertion of an interaction network based only on correlations.
- CI's role is deterministic source audit, source study replication of an established result, and preserving provenance/negative contrasts.
