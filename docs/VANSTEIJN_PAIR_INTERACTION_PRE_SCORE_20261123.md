# Van Steijn 2025 pairwise-plant formation-response contract (fixed before scoring)

> **Source-method erratum (2026-11-23, after scoring):** The original article's Methods and Figure 1 explicitly establish that solo, simultaneous and late focal plants **each grew 10 weeks**; only the pre-established neighbor grew 13 weeks. Any statement below implying focal growth duration differs is **superseded** by [the methods correction](VANSTEIJN_FOCAL_AGE_METHOD_CORRECTION_AND_NET_NEIGHBOR_EFFECT_20261123.md). The original statistical protocol and date ordering are preserved for auditability.

## Source and strict boundary

van Steijn et al. (2025), *Ecology* 106 e4528, DOI https://doi.org/10.1002/ecy.4528. Original public Zenodo **v2** https://doi.org/10.5281/zenodo.10512662 (2024-01-15). Exactly five files are source-verified in Actions https://github.com/zuizui0223/knk/actions/runs/45094845862: `alone.csv` (74 rows), `simultaneous_arrival.csv` (206 rows), `second_arrival.csv` (206 rows), `traits.csv` (15 rows), `R_code.Rmd` (585 lines). Respective MD5s are pinned in `benchmarks/probe_vansteijn_competition.py`. Original data source has 5 Grass, 5 Forb, and 5 Legume species. **Never commit raw data**.

## Scientific question distinct from KNK priority-storage results
Within the **directed plant-pair response matrix**, how often and how consistently does arrival timing alter measured focal biomass relative to simultaneous planting? Is the timing response chiefly a property of the focal species, of the early neighbor, or of particular ordered pairings? Does the arrival penalty differ between shoot and root tissues?

The source authors already published the strong role of competitive traits and functional-guild identities. The additional KNK contribution is a source-audited, matrix-oriented *comparative bridge* to the assembly-history storage theory, not a claim of discovering priority effects or direct competition.

## Units and estimands

**Unit:** one directed focal species `i` with specified neighboring species `j` under one planting condition; each directed pair normally has **one pot**, so species-pair rows do not provide independent replicates for a causal inference. Controls `alone.csv` have ~5 pots per species; use focal-specific mean of numeric, nonnegative `biomass.total.CA`. Two concurrent experiment outcomes (named partner) use same field level names and species roster. Do not include `i=j` self pairs. Inspect overlap of directed pairs in simultaneous vs second-arrival files and STOP if not sufficient / duplicated. If any total biomass = 0, use the authors' RII where denominator > 0 and flag undefined, never arbitrarily impute.

**Published response definitions (verify R_code.Rmd lines 60–104):**

\[
RII_{sim,ij} = (B_{sim,ij} - \overline{B}_{alone,i})/(B_{sim,ij} + \overline{B}_{alone,i})
\]
\[
RII_{late,ij} = (B_{second,ij} - B_{sim,ij})/(B_{second,ij} + B_{sim,ij})
\]

Calculate total, aboveground and belowground variants with same matched species-pair keys. **Important:** `RII_late` is the *additional realized effect of planting second relative to the same focal/neighbor pair planted simultaneously*. It **is not a direct interspecific interaction coefficient**, since focal age/timing and neighbor presence both change. `RII_sim` is a competitive outcome relative to an alone control, but if age, initial density and treatment scheduling are not matched for all interventions, do not claim an unbiased direct causal neighbor effect across schedules.

## Primary outcome and falsifiable checks

1. Exact one-to-one directed pair join, 15-species trait/guild join, check original R-code's means and individual denominators. Outcome `RII_late_total` for observed complete pairs. Report distribution, proportion negative/positive/zero and count, plus per-focal species mean and per-neighbor species mean.
2. **Focal-cluster aggregation**: compute mean `RII_late_total` for each of 15 focal species, each based on ~14 neighbors. To test if the mean schedule penalty is negative, use an exact **2^15 sign-flip randomization** on these 15 focal-species means (one-sided only under symmetric species-cluster null); report paired-species 95% bootstrap CI, 20,000 species resamples, source diversity n=15. Do not treat ~206 dyad rows as 206 independent p-value units. The 15 species are one nonrandom experimental pool.
3. **Shoot vs root differential**: compute each directed pair's `RII_late_above - RII_late_below`, average within focal, same sign-flip logic and species-level bootstrap. Interpreted as a difference in organ-specific *timing responses*, not absolute allocation or change of interaction direction.
4. **Early-neighbor identity test**: remove focal-specific mean late penalty and ask whether same neighbor identity recurs across focal taxa, using explained fraction of residual sum of squares from neighbor means. Null: shuffle **neighbor labels within each focal row** (9,999, seed 20261123), preserving observed values and number of pairs. This is an *exchangeability/consistency diagnostic*, **not a randomized manipulation p-value**. It tests transferable identity-specific neighbor influence.
5. **Guild-level analysis**, exploratory: 5 grass/forb/legume classes for both focal and neighbor; report a 3×3 focal×neighbor mean outcome table and compare neighbor guilds after focal centering with within-focal shuffle. Do not interpret as a mechanistic trait interaction coefficient.
6. **Reciprocal directed pairs**: count unordered pairs with both directions; report asymmetry \(RII_{ij}-RII_{ji}\), proportion differing in sign, and correlation between directed responses. These are descriptive single-pot dyad comparisons, no dyad-level causal p-value.
7. **Prediction check**: compare leave-one-unordered-pair-out baseline (focal-only) vs additive (focal + neighbor) shrinkage prediction, with all preprocessing/fitting inside held-out pair folds. Stable positive gain would be useful as a robustness diagnostic of partner identity contribution, but is not independent external validation or novel evidence of interaction formation.

## Strongest possible, and impossible, conclusions

Even a significant negative mean `RII_late` establishes a widely known growth cost of later planting **in the presence of neighbors**; it does not isolate competition from time available to grow. Partner identity, reciprocal asymmetry and organ allocation patterns offer **an operational map of context-dependent competitive outcomes** but not direct `alpha_{ij}(H)\) of a dynamical population model. To directly estimate a formed interaction, run an age-/density-matched factorial **partner-present vs partner-absent × partner-arrival-order** (and resource spatial-access manipulation), ideally collecting absolute growth trajectories. Compare to KNK rank-vs-label storage only as a mechanistic *hypothesis* unless the same experimental panel measured both.
 
## Repository boundaries

This is a separate post-submission exploratory branch: `analysis/interaction-formation-20261123`. Frozen KNK v28–v38, Ecology Letters submission files, source fingerprints, and PR #3 Dawson-Glass / PR #4 POEM remain unchanged.
