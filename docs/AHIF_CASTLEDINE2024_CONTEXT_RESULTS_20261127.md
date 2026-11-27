# Castledine 2024: finite density-response portability — verified 2026-11-27

Source DOI 10.1099/mic.0.001489, author commit 14c9ac7e3b8ab05447fbe19519b1e3c8f94afea3. Corrigendum 10.1099/mic.0.001528 corrects genome summary. Frozen prescore 1bb11d322a71b1dc9993fcb58b6a56c31085c98c. Verified [run 45939572589](https://github.com/zuizui0223/knk/actions/runs/45939572589) job 140310529740, commit eda4153e31c1bc9b1b11b50d8a8a23afa90aef68, artifact 15172062623 SHA256 a56b36b36a1a5a7f2a8965400d664ddb53ef1d5f6709c3ef2d716dcb19b3bf3a. Exact processed interaction CSV and author analysis script Git blob hashes verified. All 75 focal×nonempty-partner subsets unique and complete, no zero densities; author higher-order summary identity observed−expected reproduced. Data are author's processed CFU-relative-to-monoculture summaries, not raw biological replicates. A value of 1 for empty background is definitional, not additional observation.

## Context-dependent finite marginal response
For focal i, added species j, pre-existing background S (excluding i,j), define M(i,j,S)=w_i(S+j)−w_i(S), where w is author relative density. 20 pair-alone contrasts; 140 nonempty-background contrasts sharing the same endpoints and monoculture references. Relative density is not an instantaneous growth rate, per-capita coefficient or identified direct mechanism. Starting focal inoculum nominally equal across combinations; total inoculum changes with richness, and indirect effects are included. 

| Background size | Same sign | Opposite | Agreement | Same sign, threshold .05 | Opposite | Unresolved | Agreement among defined |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 47 | 13 | 78.33% | 41 | 10 | 9 | 80.39% |
| 2 | 42 | 18 | 70.00% | 33 | 15 | 12 | 68.75% |
| 3 | 12 | 8 | 60.00% | 8 | 6 | 6 | 57.14% |
| All | 101 | 39 | 72.14% | 82 | 31 | 27 | 72.57% |

12/20 directed pairs reverse at least once (9/20 after fixed .05 absolute relative-density effect threshold). Threshold not inferential equivalence. Four positive and 16 negative pair-alone point estimates; no statistical significance claimed. The decrease in agreement with background richness is descriptive within one five-species system. No independent edges or iid uncertainty.

## Partner-identity transfer
Fixed forecasts for target w_i(S+j) using known observed w_i(S). Pair-balanced MAE, lower is better. Primary 120 rows at background sizes 1–2; size3 has only one possible added species and cannot assess identity. 

| Forecast | MAE |
|---|---:|
| No change from background | .2725546064 |
| Add pair-alone delta, clip >=0 | .2301438718 |
| Multiply by actual partner pair-alone response | .2242231984 |
| Multiply by mean response over remaining candidate partners | .2155356480 |

The actual-partner multiplicative forecast improves on no-change but is worse than a focal/richness/background-availability matched anonymous partner baseline by .0086875504 MAE; no assertion of universal superiority for anonymous effects. Actual identity beats anonymous in only one of five focal averages. In the secondary full endpoint forecast, using only pair data, all 50 identity-comparable 3/4-species focal endpoints: no response .5314309033, actual pair product .2551914064, anonymous geometric-mean product .2413103307. Five-species full endpoint identity comparison tautological and kept separate. No fitting, no independent validation, no p-values. The source-specific partner detection result from earlier studies concerns exposure/arrival, not growth-effect portability; their gains are not the same metric. This study is evidence that contact opportunity and net response are distinct prediction stages. It does not establish that a specific species interaction is reversed at the molecular level or that higher-order effects alone cause the discrepancy. Original paper already examines higher-order deviations; contribution here is explicit background marginal portability and matched anonymous competitor comparator.

## Registration boundary
Derived artifact contains finite_density_endpoints.csv (75), context_marginal_responses.csv (160, 20 pair-alone + 140 background-present), full_endpoint_forecasts.csv (55), summary.json. Exact source/run receipt validation/ahif_v1/castledine2024_context_summary_20261127.json. Evidence class: author_processed_CFU_relative_density. Prior v0.6 signed interaction coefficient table unchanged and empty. No formal meta-analysis eligible effects or raw replicate uncertainty. Venturelli EV1/EV3 source gate 45939324000: both official workbooks verified, but only relative composition summaries; no matched OD and monoculture source matrix, absolute response STOP. Do not treat this as biological failure or claim that raw data are unavailable everywhere. No manuscript or main branch change.