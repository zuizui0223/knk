# Method correction and new experimental interaction estimand (2026-11-23)

**This is an explicit, post-primary-analysis SOURCE-METHOD CORRECTION.** Earlier `docs/VANSTEIJN_PAIR_INTERACTION_PRE_SCORE_20261123.md`, `docs/KNK_INTERACTION_FORMATION_VANSTEIJN_RESULTS_20261123.md` and aggregate summary included a mistaken suggestion that the focal's growth duration differed between the sequential and simultaneous conditions. That claim is **superseded** and cannot be reused.

## Source correction: focal growing time is age matched

van Steijn et al. (2025), *Ecology* 106:e4528, DOI https://doi.org/10.1002/ecy.4528; published methods and Figure 1 (public full text https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecy.4528).

- **Simultaneous mixed pair:** focal and neighbor both introduced in week 3, harvested week 13 -> focal receives **10 weeks** growth.
- **Sequential mixed pair:** early neighbor introduced in week 1, later focal in week 3; both harvested in week 13 -> focal receives **10 weeks** growth, early neighbor **13 weeks**.
- **Focal alone control:** focal planted in week 3, harvested in week 13 -> focal receives **10 weeks** growth.

Thus the authors deliberately match the focal's elapsed growth period across all three treatment types. A history-dependent change in *net neighbor effect* on focal biomass **is identifiable at the measured endpoint**, conditional on validity of single-pot species combinations and the control sampling. The experiment does **not** independently manipulate the early neighbor's size, resource depletion and habitat conditioning; those are bundled mechanisms of the three-week head start. It also cannot identify a time-varying per-capita dynamical interaction parameter \(\alpha_{ij}(t,H)\) from one final harvest.

The original authors' exact R code already defines RII_seq=(B_late−B_sim)/(B_late+B_sim) and RII_sim=(B_sim−B_alone)/(B_sim+B_alone). **Do not misinterpret RII_seq as a direct interspecific coefficient**, but the matched solo control allows an additional empirically meaningful response:

## Additive new-to-KNK net neighbor effect audit — post hoc after source correction

For each of the **202 aligned directed species pairs**, with solo mean \(A_i\), simultaneous focal biomass \(S_{ij}\), and 3-week-preestablished-neighbor focal biomass \(L_{ij}\):

\[
C_{ij}^{sim}=(S_{ij}-A_i)/(S_{ij}+A_i),
\quad
C_{ij}^{pre}=(L_{ij}-A_i)/(L_{ij}+A_i)
\]
\[
\Delta C_{ij}=C_{ij}^{pre}-C_{ij}^{sim}.
\]

This is a **history-specific realized net neighbor effect on final focal biomass** (relative to an age-matched solo focal). It is not a mechanism-specific direct pairwise competition coefficient. Positive C indicates a focal grew more with that neighbor than the solo-species mean, but *one pot per directed pair* cannot prove mutualistic facilitation versus stochastic variance.

The post hoc tests must:
1. Reuse the identical source MD5 checks, 202 exact pair join and 15 focal species clusters.
2. Report counts of focal combinations with positive vs negative C in both scenarios, sign transitions, and per-focal means of ΔC, with **15-species exact sign flip** plus 20,000 species-cluster bootstrap (do NOT treat 202 dyads as independent species-level replicates).
3. Repeat the prior within-focal **neighbor-identity pseudo-null** and the more conservative within-focal×neighbor-guild identity pseudo-null for ΔC, with the same 9,999 shuffles; a guild mapping placebo is descriptive.
4. Report **global interaction formation effect** \(\Delta C\) and partner-specific repeatability separately; avoid interpreting label-permutation p as randomized-neighbor causal p.
5. Assess how much the changed net neighbor effect is predictable out of pairs from focal-only vs focal+neighbor species labels using the same held-out unordered-pair folds and fixed ridge penalty.
6. No retrospective declaring this newly chosen endpoint “preregistered”, and no editing KNK frozen v28–v38.

This amendment strengthens the interpretation of source experiment's history-dependent **neighbor impact** but does not establish niche pre-emption vs modification, which need measurements of resource occupancy, early-neighbor abundance and intervening host/soil state.
