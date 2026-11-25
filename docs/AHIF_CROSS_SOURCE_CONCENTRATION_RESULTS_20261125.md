# Cross-source response concentration: lagged response history transfers better than reference dominance

2026-11-25 JST. Source-locked scoring run https://github.com/zuizui0223/knk/actions/runs/45560465895, job 139022660041; all synthetic tests, source hashes and grids, archived contrast checks, real scoring and derived artifact upload succeeded. Prescore commit e59dc577ebb21de50a03979caf61f6d48e42a534. Full features and prediction cells in artifact 14996984477 (sha256:edacd903473b22dd60810f5a5fec48c659e00de8549e15fd0ccd9476bee23a15); exact aggregate log summary archived validation/ahif_v1/cross_source_concentration_summary_v1.json. All four sources already represented in prior research diagnostics, no new independent study or ecological effect row. N=25,680 split-by-time-by-role-by-direction cells, 2,957 undefined because no supporting gross response; 22,723 eligible. These are correlated derived cells, NOT independent observations. 19 source contrast panels, four experiments (Kansas, Narwani ponds, POEM, van Steijn). Independent within-arm splits for Narwani/POEM, matched block splits Kansas/van Steijn. Whole sources held out, equal source and within-source stratum weights; directly checked equal-source aggregate MSE arithmetic. Same 4 source families as prescore. Kansas EDI access still blocked; published roster linkage provisional. POEM and van Steijn potentially related research groups, so conservative grouped holdout included. Original source provenance preserved by pinned parsers.

## Primary and sensitivity scores

Top supporting contributor fraction, conditional on defined positive gross support in either signed direction. M0 intercept; M1 nominal and previous reference richness, absence flag and direction; M2 adds lag reference top1 and Simpson; M3 adds lag training response support count/top1/absence to M1; M4 adds lag reference dominance to M3. Same training-only fixed ridge lambda .1, clipped predictions. Positive gain is lower held MSE. Both-arm profile sensitivity is not pretreatment. Restriction to >=2 supporting taxa uses observed target eligibility, descriptive only.

| Validation | Eligible cells | M0 | M1 | M2 | M3 | M4 | M2 vs M1 gain | M3 vs M1 gain | M4 vs M3 gain |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| control_source_all | 22723 | 0.064128 | 0.059145 | 0.059270 | 0.047820 | 0.048629 | -0.21% | 19.15% | -1.69% |
| both_source_all | 22723 | 0.064128 | 0.053440 | 0.051936 | 0.046577 | 0.045923 | 2.81% | 12.84% | 1.40% |
| control_cluster_all | 22723 | 0.062765 | 0.072350 | 0.073470 | 0.054376 | 0.056329 | -1.55% | 24.84% | -3.59% |
| control_source_multi | 18273 | 0.050949 | 0.052268 | 0.051446 | 0.044848 | 0.045142 | 1.57% | 14.20% | -0.66% |

Primary reference-only dominance gives -0.21% gain over richness/direction and -1.69% beyond lag response structure. Lag response structure gives +19.15% over M1. Conservative joint research-group holdout retains +24.84% for M3 over M1, dominance -1.55%/-3.59%. Restricting to >=2 supports retains +14.20% for M3, dominance +1.57% over M1 and -0.66% beyond M3. Mean-arm reference sensitivity gives +12.84% for M3, dominance +2.81% and +1.40%, small and not consistent with the primary. M1 vs intercept not robust: +7.77% primary, -15.27% grouped holdout, -2.59% multi-support. Do not attribute M3 improvement uniquely to response concentration; lag support richness, absence and concentration enter jointly. No inference of significant difference or generality from four source families.

| Held experiment | M1 vs M0 | M2 vs M1 | M3 vs M1 | M4 vs M3 |
|---|---:|---:|---:|---:|
| kansas | -9.94% | -2.80% | 33.53% | -4.02% |
| pond | 27.64% | 4.41% | 14.89% | 1.52% |
| poem | 13.93% | -2.61% | -3.31% | -1.22% |
| field | 2.06% | 0.74% | 20.92% | -2.34% |

M3 improves on M1 in 3/4 primary held experiments, POEM -3.31%. In grouped holdout all four relative gains positive, but reference-only baseline weaker than intercept in both grassland experiments; no claim of universally best model. Primary source-equal M3 MSE .0478196103 vs M0 .0641276264. Dominance signal not robustly transferable.

## Descriptive concentration and equal nominal roster fraction exclusion

All contrasts, observed times following first observation, both response directions scored; selected biologically interpretable directions displayed here, full table in JSON. Equal fraction rule is ceil(.1*nominal roster), up to that many positively supporting taxa selected from previous training data. Thus actual removed proportion can exceed 10% in small groups or be less if little prior support. Conditional net sign retention among positive held net responses, not a forecast accuracy or physical intervention. Pooled unweighted split-cell descriptive values; no cross-source pooling of units or claim of formal comparative effect.

| Source | Role direction | Defined gross support | Mean support count | Mean top1 gross share | Retain net sign after exclusion |
|---|---|---:|---:|---:|---:|
| kansas | other_named − | 84/84 | 10.702 | 55.64% | 57/60 (95.00%) |
| kansas | sown_list + | 75/84 | 9.147 | 47.46% | 73/74 (98.65%) |
| pond | cyanobacteria + | 1505/2016 | 1.383 | 96.45% | 577/1235 (46.72%) |
| pond | other_named_algae − | 1946/2016 | 4.956 | 72.24% | 908/1264 (71.84%) |
| poem | early_sown + | 2275/2400 | 3.057 | 75.61% | 1672/1847 (90.53%) |
| poem | nonsown + | 2398/2400 | 6.192 | 63.44% | 1859/1981 (93.84%) |
| poem | other_sown − | 2400/2400 | 5.469 | 64.64% | 1715/1853 (92.55%) |
| field | early_sown + | 480/480 | 5.075 | 47.92% | 445/462 (96.32%) |
| field | nonsown + | 464/480 | 4.287 | 64.39% | 105/204 (51.47%) |
| field | other_sown − | 480/480 | 4.744 | 52.67% | 440/472 (93.22%) |

Pond cyan increases often have only one supporting morphotype (1057/1505 defined positive-direction cells across all pond contrasts). Kansas sown increase 3/75 single-support cells. POEM early-sown increase has 321/2275 single-support cells; Swedish early-sown increase 0/480. Different taxonomic resolution/nominal roster and units; do not interpret these as universal categories or infer actual biological replacement. Normalized excess and exact-count bins in archive. Relative exclusion attenuates dependence on the arbitrary top-three comparison but does not eliminate taxonomic grain or rounding effects. No harmonized functional trait table, no test of functional diversity. Prior control state is generally measured after intervention started. No causal effect-modification claim or meta-analysis. 

## Archive scope / deviation
Full feature and prediction cells exist in verified uploaded artifact; fixed code and hash-pinned sources permit reproduction after expiration. Prescore mentioned retaining model coefficients; scorer calculates them but does not store them. Do not claim a coefficient audit or direction-of-slope inference. All model losses, all sensitivity summaries and conditional-support coverage retained. No tuning after held outcomes. Existing v1.0 effect ledger, main, frozen PRs3–5 unchanged. Next independent source priority is randomized formation/recruitment or true intervention perturbation, with common functional-trait coverage, rather than further fitting these held sources.
