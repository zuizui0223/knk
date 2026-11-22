# Independent KNK validation (post-submission) — Dawson-Glass et al. 2025

**Updated status (2026-11-22): the original CSV was supplied by the user, audited and scored. Primary results are in `DAWSONGLASS_EXTERNAL_RESULT_V1.md`; the exploratory winner/biomass follow-up is in `DAWSONGLASS_WINNER_IDENTITY_FOLLOWUP_V1.md`. This file preserves the original pre-access analysis specification and its historical access blocker.**

This is an *unmerged*, non-v39 analysis candidate. It does not reopen or modify
KNK's v28–v38 frozen results, Ecology Letters manuscript, or submission package.

## Sources

- Experiment: Dawson-Glass et al. (2025), *Ecology*, e4504, DOI 10.1002/ecy.4504.
- Primary data: https://doi.org/10.5061/dryad.sbcc2frh4 (2024-11-13 version).
- Original analysis code: https://doi.org/10.5281/zenodo.13988289.
- Dryad describes the individual-plant file:
  Dawson-Glass_ECY24-0464_assembly_data.csv
  with Pot, Species, init_height, Sp_code, Treatment, Biomass.
- Design published: 24 mesocosms (12 ambient assembly, 12 warmed assembly),
  15 species per mesocosm. Two order schedules were derived from temperature
  manipulations in a *separate germination experiment*. The mesocosms themselves
  were grown in shared greenhouse conditions, not warmed vs ambient temperatures.

## Historical access STOP (resolved later on 2026-11-22)

- The Dryad landing page and file manifest were readable on 2026-11-22.
- Direct requests for individual CSV file_stream/3629006 returned HTTP 403.
- At that initial stage the source data were **not yet locally accessible**. Subsequently, the user supplied the raw CSV (SHA-256 `1514d4501a0bdf7a61ace6c250566c5a0b8e4a24a0f4737646d182ebad7b8230`) and all audit gates passed. **The initial STOP is resolved, not a current blocker.**
- Raw CSV must remain local and must not be redistributed into this repository.
- Source checks performed *by the script once data are present*:
  24 unique pots, 15 unique species, 12 pots per order, consistent treatment
  per pot, finite nonnegative biomass, a complete pot × species panel, no
  duplicate pot/species entries, positive total biomass, and SHA-256 of raw input.
  A STOP on missing or duplicated records requires human reconciliation with
  the actual source; do not silently sum/zero-fill.

## Reproducible execution from the audited local CSV

Run from the KNK checkout with Python >= 3.10 and NumPy:

    python benchmarks/run_dawsonglass_storage_post_submission.py \
      --assembly-csv /path/to/Dawson-Glass_ECY24-0464_assembly_data.csv \
      --output /path/to/private_results/dawsonglass_report.json \
      --permutations 9999 --seed 20261122

Do not commit raw data or preliminary results to main.

## Predefined comparisons

1. **Primary (within mesocosm):** labeled 15-species biomass proportions
   versus the *same* proportions sorted descending with all species labels
   discarded. For each target, compute Bray–Curtis, distance-based partial R2
   for assembly-order assignment, and design-respecting null.
   Null-calibrated excess: E = R2(observed) − median R2(permuted).
2. Shared stream of 9999 whole-pot permutations, seed 20261122. The paper does
   not report blocks; verify actual assignment before interpreting unblocked
   permutations as causal randomization. Same permutation vector used for both
   representations. Report all upper-tail p values and null 95% intervals.
3. **Predictive companion:** fixed-penalty logistic model, each of 24 pots held
   out once, preprocessing within fold. Report held-out log-loss improvement
   over training-fold prevalence baseline for both representations. No naive
   independent test/p-value is attached to these overlapping folds.
4. **Sensitivity only:** original raw biomass matrix, pot-total biomass,
   identity-free rank-1 maximum and Shannon entropy. Do not search for
   additional positive endpoints after scoring.

## Claims and limitations

The source's raw-biomass composition result (PERMANOVA pseudo-F 16.61,
p < .001) is a prior published reference, **not a new KNK result**.
Our proportions are a different estimand and need not reproduce its statistic.

**Do not equate E_labeled / E_sorted with percentage information retained.**
Sorting is deterministic: genuine mutual information cannot increase after
sorting, although null-calibrated R2 can rise due to changed metric geometry.
This distinguishes predictive historical recoverability from ecological
mechanism.

Even if results are significant, experimental arrival order also changes time
available for growth. The source paper acknowledges that these data cannot
distinguish competitive effects from timing-only growth advantages. No
niche-preemption / niche-modification mechanism is identified by this reanalysis.

## Current gate

- [x] New independent dataset identified and source publication audited.
- [x] Source schema verified from Dryad README.
- [x] Independent analysis code committed to a separate branch.
- [x] Synthetic tests to verify identity-only and architecture-rich signatures.
- [x] User-supplied original CSV audited and SHA-256 recorded.
- [x] Pot-level audit PASS.
- [x] Original-source raw-biomass PERMANOVA pseudo-F matched to 16.61.
- [x] Predefined real-data endpoints scored; see standalone result report.
- [ ] Independent result/interpretation review.
- [ ] Submission-ready main manuscript completed before using this as an additional-study result.
