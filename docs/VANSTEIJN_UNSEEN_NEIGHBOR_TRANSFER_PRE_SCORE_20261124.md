# Transfer prediction of history-dependent neighbor effects to unseen plant species

**Analysis-plan freeze: 2026-11-24; this is a new exploratory question motivated by prior results, not a preregistered independent experiment. This document is committed BEFORE the strict species-out prediction code is scored.** Source and original 2025 study already link competitive responses to traits; novelty would be stricter prediction validation, not trait discovery. No edits to KNK v28–v38, the Ecology Letters letter or other Draft PRs.

## Fixed source and outcome

Van Steijn et al. (2025), *Ecology* DOI https://doi.org/10.1002/ecy.4528, primary data [Zenodo v2](https://doi.org/10.5281/zenodo.10512662), five MD5-pinned original source files already audited in `benchmarks/probe_vansteijn_competition.py`. Directed 202 matched focal×neighbor pairs, 15 plant species: five grass, five forb, five legume. Each directed pair normally has a single pot. Solo controls are approximately five per focal species.

**Outcome** is age-matched, net neighbor pre-establishment effect on focal total biomass, matching `benchmarks/score_vansteijn_age_matched_neighbor_effect.py`:
`C_sim=(B_sim−B_alone)/(B_sim+B_alone)`, `C_pre=(B_pre−B_alone)/(B_pre+B_alone)`, `deltaC=C_pre−C_sim`. Focal growth duration was 10 weeks in both mixed conditions and solo controls. Early neighbor received an extra 3-week establishment advantage. This is a total priority-history response, NOT a dynamical per-capita competition coefficient \(\alpha_{ij}\).

## Truly out-of-species prediction (not leave-one-pot-out)

For each of the **15 species** as prospective early neighbor `j`, exclude **all directed-pair rows in which `j` appears in EITHER role** from training. Predict only held-out rows where `neighbor=j`; each observed directed pair serves as a test outcome **once**, no reverse-pair contamination. The held-out species' fixed trait values and guild label may be used as predictors. No response (neither focal nor neighbor) involving that species may be used in fitting. Report per-neighbor test size, per-neighbor MSE and overall paired comparisons; species-fold CV is not independent external data and 15 folds share training data.

**Do NOT compare this evaluation to prior leave-one-unordered-pair-out RMSE**, which allows the same species to occur in training and measures a different transfer problem.

## Fixed traits, transformations, and models

Primary fixed **morphology/strategy trait set** from public `traits.csv`: `SLA`, `LA`, `RMF`, `Plant height`, `Lateral spread`, `LDMC` (six source columns, no replacement or selection after looking at outcomes). Use `log1p(x)` on nonnegative `SLA`, `LA`, `Plant height`, `Lateral spread`; `RMF` and `LDMC` stay on original scale. All species traits must be finite/nonnegative. Fit per-column standardization **using the 14 training species, never the test species**. Any zero-variance column gets zero standardized values; record it rather than use a test-derived SD.

**Exclude the source `GR` (growth rate) from the confirmatory predictor set**, because GR is derived from the same solo-grown biomass that enters the \(\Delta C\) denominator; its inclusion could exploit mathematical coupling or post-experiment performance. Include it only as a clearly marked post hoc sensitivity, log1p transformed.

Fit the following fixed linear ridge models to training directed pairs, with an unpenalized intercept and a fixed penalty **10** on all non-intercept coefficients. No hyperparameter tuning or model selection using test folds. For every holdout species use the *same test rows* for all models.

1. `intercept`: training-response mean
2. `guild`: focal and neighbor three-level grass/forb/legume labels, one-hot reference-coded
3. `focal_morph`: six standardized focal traits
4. `both_morph`: six focal + six earlier-neighbor standardized traits
5. `both_morph_abs_gap`: those 12 plus six absolute standardized trait differences, a proxy for non-additive trait matching
6. `guild_both_morph_abs_gap`: guild one-hot plus the 18 trait values

Targets: primary `deltaC_total`. Predeclared sensitivity: the source's `RII_late_total` and `RII_sim_total` using the same folds/models; and GR included in a duplicate version of the trait models.

## Decision rules, controls and uncertainty

Evaluate squared errors with one observation per directed held-out neighbor pair. Report overall RMSE and species-equal mean squared error (each of 15 held-out early-neighbor names one group). Two essential **paired** comparisons: `both_morph_abs_gap vs guild` (does quantitative morphology beat coarse guilds on unseen neighbors?) and `both_morph_abs_gap vs focal_morph` (do traits of an unseen earlier neighbor improve over the focal alone?).

Construct a descriptive **15-neighbor-cluster bootstrap** (20,000 draws, seed 20261124) of the *difference in species-average MSE* for each main comparison; report 95% percentile intervals and fraction of 15 held-out neighbors individually improved. Neighbor species are nonrandom and phylogenetically nonindependent; do not call the interval a phylogenetically calibrated cross-flora CI or a causal p-value.

As an *explicit mapping specificity check*, permute the **earlier-neighbor six-trait vectors among the five species in each functional guild**, while keeping focal traits and both species' guilds unchanged; 199 independently shuffled global species mappings, seed 20261125. Repeat strict species-out CV for `both_morph_abs_gap`. Report the shuffled null distribution and fraction of 199 mapping trials with MSE ≤ observed. This is **not** a formal experimental randomization p-value, and the model is chosen based on prior result.

## Interpretation gates

- If unseen-neighbor morphology does not improve compared to guild or focal-only baseline, **do not claim a transferable trait-based interaction-formation law**. Strong in-sample neighbor identity can be idiosyncratic or require unmeasured traits/early-neighbor size.
- If morphology improves, it supports **within-this-15-species-panel prediction to an unobserved neighbor-role species** under these source conditions, not a macroecological law or a causal mechanism.
- Fixed trait measurements are from original solo-grown control conditions, rather than remotely observable prior-to-experiment traits; state this domain/availability limit.
- One source publication; 202 directed pair records are **not 202 independent study effects**. No random-effects meta-analysis.
- The original paper already demonstrates the role of traits in in-sample priority effects; KNK only tests the more stringent species-out generalization.
