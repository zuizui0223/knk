# Kansas same-target role forecasts: direct aggregation does not beat simple historical averaging

2026-11-25 JST. Source-locked scorer 8b483aad38712da96060b380046df5e6f6b23875, run https://github.com/zuizui0223/knk/actions/runs/45544096731, job 138969396575; three synthetic tests, original curated file hash checks, exact previously frozen taxonomic mapping and 300 prior species model-by-cell SSE comparisons (max difference 0), real scoring and artifact upload passed. Full derived JSON (all four representations ×60 prediction cells) in [artifact 14990766503](https://github.com/zuizui0223/knk/actions/runs/45544096731/artifacts/14990766503), sha256:7572ab588ceee5d192517281604002075eb976c6ca827db6e0df6643cfef9d6d, 134849 bytes; Git-tracked summary at validation/ahif_v1/kansas2020_group_forecast_summary_v1.json. Full artifact expires under Actions retention; exact frozen code can reproduce it. The Git summary deliberately excludes cell-level predictions, not recovered from the artifact. Same 120 direct curated observations, 100 taxon labels, 4 matched source blocks, 15 source age codes. Hold out 2 entire blocks and all current/future time data; 6 overlapping oriented splits ×10 expanding-window forecast times (codes 6–15). Not 60 independent samples. Original EDI source unavailable; 27 published sowing list matches (23 literal punctuation-normalized, four synonym/orthography), 62 other identified labels, 11 unresolved. Source paper Foster 2009 DOI 10.1890/08-0849.1; curated Ladouceur 2020 DOI 10.1111/ele.13566. Original lineage probable, not independently audited. Same target vector within each representation; no comparison of raw species SSE against grouped SSE. No independent study, p values or causal inference.

## Primary broad3 target, total held SSE (same mass units and held cells)
| Model | SSE | Gain vs zero | Whole-community total SSE | Opposition forecast balanced accuracy |
|---|---:|---:|---:|---:|
| zero | 452466.655 | 0.000% | 280348.768 | 0.5000 |
| previous | 478596.749 | -5.775% | 498322.472 | 0.5686 |
| cumulative | 386858.615 | 14.500% | 354680.069 | 0.5327 |
| bottom_history | 422639.322 | 6.592% | 383484.075 | 0.4902 |
| bottom_state | 469894.566 | -3.852% | 492640.177 | 0.5523 |
| direct_history | 460568.851 | -1.791% | 448767.744 | 0.6111 |
| direct_state | 463496.232 | -2.438% | 508260.417 | 0.5948 |
| direct_history_scaled | 615276.181 | -35.983% | 670665.334 | 0.5163 |
| direct_state_scaled | 713443.695 | -57.679% | 762015.203 | 0.4510 |

Primary direct history ridge vs bottom-up species history ridge gain -8.974%: FAIL. Direct state ridge vs bottom-up state ridge +1.362% but direct state worse than zero by 2.438% and cumulative by 19.810%. Reduced per-row-scaled penalties worsen both comparisons. Historical cumulative role response is best on every representation (not a prespecified universal best method). For previous and cumulative predictions, summing individual taxon predictions exactly equals direct group means by linearity; cannot interpret cumulative success as evidence of group-specific modeling advantage. Broad cumulative SSE gain versus zero 14.500%; sown SSE reduced, other-named SSE increased. Whole-community total SSE best zero across all nine tested models, so no improved total mass forecast or biomass equivalence.

## Resolution and taxonomic sensitivity
| Target | Group counts | Direct history gain vs bottom-up history | Direct state gain vs bottom-up state | Scaled history gain | Scaled state gain | Best tested magnitude model | Cumulative gain vs zero |
|---|---|---:|---:|---:|---:|---|---:|
| broad3 | 27/62/11 | -8.974% | 1.362% | -45.579% | -51.831% | cumulative | 14.500% |
| exact23 | 23/62/15 | -0.800% | 13.561% | -16.376% | -18.091% | cumulative | 8.494% |
| sown_guilds6 | 1/6/5/15/62/11 | -5.307% | 5.529% | -25.037% | -22.963% | cumulative | 6.438% |
| guilds_to_broad3 | 27/62/11 | -5.449% | 1.622% | -28.055% | -36.228% | cumulative | 14.500% |

The six-category partition includes published functional guilds of sown taxa only (AF=1, C4-PG=6, PLF=5, PNLF=15), plus the other-named and unresolved roles. It is not a whole-community functional trait partition. Collapsing the fine resolution onto identical broad3 target reduces direct history SSE relative to coarse direct history but still worse than bottom-up history and cumulative baseline. Exact23 sensitivity moves four alias matches to unresolved. Mild positive state-specific gains are model/penalty dependent, do not beat cumulative. No arbitrary-group control or test of biological role-specific superiority. 

## Direction forecast (distinct from conditional held sign reproduction)
Observed broad opposition (sown positive and other negative) in 51/60 dependent target cells. Always predicting opposition gives 85% ordinary accuracy and 0.5 balanced accuracy. Previous-age sign predicts sown sign correctly in 60/60, other sign 44/60, opposition 44/60; balanced accuracy 0.5686. Direct history balanced accuracy 0.6111 (TP34 FN17 TN5 FP4), direct state 0.5948; exact23 direct state 0.6405 secondary. High prevalence and repeated cells mean no significance/generalization claim. This is forward sign forecast from past other blocks, not the earlier same-time fixed-role descriptive score. Previous opposition occurrence itself was not uniformly forecast, and response magnitude remains difficult.

## Conclusion / stop
Fixed treatment roles show persistent observed opposition at source codes 7–15, but direct group ridge forecasting did not outperform the historical mean on any tested target. This is a falsification of automatic forecast improvement from ecological grouping in this panel, not proof that ecological roles never improve prediction. Freeze this source's predictive exploration; no further outcome-tuned model search. Need original independently verified arrival/recruitment manipulation and complete trait classification before a general mechanistic statement. Keep primary AHIF-DB v1.0 unchanged (36 study records, 38 panels, 161 effects, 14 original outcome audits); no new study or causal bridge. Main and frozen PRs3–5 unchanged.