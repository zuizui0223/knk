# Gibbs 2026 grouped-target sensitivity — verified 2026-11-27

Post-score diagnostic frozen 67a7d2fc5862350eace03602b58bdfe2de3f88f1 before new outcomes; prior per-label scores and source audit replayed exactly. Run 45947718614 SUCCESS at bf3d195cc71cf016484c28cb62b796fcd2033180, job 140336140010, artifact 15175188406 SHA256 b1f34a9cc139ec1452f416c0e3829f1748b3afe199f03ce84f881bb8772edc82. Same 430 complete matched label contrasts, 256 fixed partitions of 15 labels into 10 training and 5 held; no skipped folds. Each held group has 1–5 complete contrasts. 30 focal×competitor-pair groups and 7680 dependent group-fold outcomes, NOT independent observations or new study. For AC/UR-only sensitivity, 20 groups, 5120 dependent outcomes and 2–5 held contrasts/group. 

| Measurement and target | No change MAE | Other-pair mean | Known pair | Same pair other focal | Held-pair additive identity |
|---|---:|---:|---:|---:|---:|
| Author raw, grouped mean | .314856 | .276519 | .294464 | .283129 | .263404 |
| Author raw, density-standardized mean | .318156 | .261660 | .298887 | .287364 | .252329 |
| Hypothetical SA ×10 diameter, raw | .428168 | .364354 | .373807 | .342016 | .338293 |
| Hypothetical SA ×10 diameter, standardized | .426459 | .348279 | .369654 | .337338 | .327703 |
| AC/UR only, raw | .401849 | .326201 | .368528 | .309177 | .322948 |
| AC/UR only, standardized | .401133 | .305317 | .363460 | .310703 | .306857 |

All fixed nonzero group-mean forecasts beat zero in all six configurations, unlike single-label raw prediction. Held-pair additive identity beats no-change by 16.34% and anonymous other-pair mean by 4.74% under primary measurement/grouped-target diagnostic; it is not consistently better than other models (AC/UR-only raw cross-focal better; standardized anonymous slightly better). The held-pair model is a fixed ridge on the five constituent competitor identities with intercept, trained on other pairs and other replicate labels; does NOT hold out constituent species. Density-standardized target uses slopes estimated only from training labels/other pairs; not causal control or independently measured baseline. Averaging multiple held labels reduces unpredictable variation and changes the estimand; this is not an independent confirmation of an ecological scale mechanism. No inference from repeated folds; no claim of general cross-system prediction.

## Measurement caveat
Salvia calibration table explicitly labels diameter mm, range 4.18–23.68 (median 9.95); nonbare head data have no explicit diameter unit, range .2–2.7 (median .9), 345 heads in matched plot set. The original author R code applies calibration to recorded values without conversion. No verified unit crosswalk; discrepancy could be unit mismatch or differences in sampled head size. Preserve original as written, and report fixed hypothetical ×10 diameter and exclusion of Salvia, neither selected by fit. The group-averaging gain remains in AC/UR-only sensitivity, but identity advantage is small/inconsistent. No raw values or source code republished. 

All 30 model×lane×scale records archived in database/ahif_mechanism_v0/gibbs_group_mean_sensitivity_20261127.csv; full JSON validation/ahif_v1/gibbs2026_group_mean_summary_20261127.json. Prior SQLite v0.8 not updated in this PR. No signed coefficients, causal mechanism, manuscript modification or formal meta-analysis.