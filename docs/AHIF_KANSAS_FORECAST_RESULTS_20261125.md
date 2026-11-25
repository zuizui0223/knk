# Kansas hay meadow — directional recurrence but poor magnitude forecasts

2026-11-25 JST. Curated original two hash-pinned Ladouceur synthesis files, article https://doi.org/10.1111/ele.13566 and Figshare https://doi.org/10.6084/m9.figshare.12319682.v1. Verified run https://github.com/zuizui0223/knk/actions/runs/45541282653, job 138960137421, commit a7cdd253cf0441d7d1e2bbd3eff4f2e19108c97d; three unit tests, source checks, grid checks, score and artifact upload passed. Prescore docs/AHIF_KANSAS_FORECAST_PRESCORE_20261125.md. Original Table 1 of Ladouceur identifies Kansas Hay Meadow (Foster 2016a), abandoned hay meadow, 41 added taxa, one site, years 1–15. Archive source age codes agree with that table for this selected experiment; absolute calendar dates and original experimental design/plot linkage have not independently been verified. Do not identify other Foster hayfield publications as the same experiment merely from the title. 120 direct-mode plot-state records, four matched site/block labels, one seeded and one control state per block and source age; 100 archived taxa. Six oriented 2+2 partitions, 10 forward forecast times (source codes 6–15), 60 dependent scoring cells. No current/held block data used in fits, features or seen-taxon mask. Fixed retrospective species roster, unseen-in-training taxa forecast zero. Model comparisons descriptive, no p-values. Selected this experiment after prior results, not independent confirmation.

| Model | Total held SSE | Gain versus zero | Mean held increase | Mean held decrease | Mean held min | Passing ages among forecast 6–15 |
|---|---:|---:|---:|---:|---:|---|
| zero | 182196.692 | 0.000% | 0.000 | 0.000 | 0.000 | none |
| previous | 312648.088 | -71.599% | 26.731 | 14.606 | -0.794 | 9,10,12,13,14 |
| cumulative | 239189.037 | -31.281% | 25.351 | 10.773 | -3.271 | 9,10,12,13,14 |
| history_ridge | 217467.471 | -19.359% | 20.202 | 5.624 | -8.420 | 9,10,12,13,14 |
| state_ridge | 237667.702 | -30.446% | 22.497 | 7.919 | -6.125 | 9,10,12,13,14 |
| concurrent_sign_reference | not a forecast | NA | 27.092 | 11.968 | -3.926 | 9,10,12,13,14 |

All four nonzero forecasting methods meet strict descriptive directional gate at the same five ages 9,10,12,13,14, fail at 6–8,11,15. Same-age held-block sign reference shows same pattern, with 13 marginal (mean minimum +.629). Exhaustive split enumeration revises earlier 256-random-split result from 4 passing Kansas ages to 5; code13 was close to zero and sensitive to split weights. This is not an independent positive result or threshold redefinition; all six orientations have equal weight now. No inference of time trend or stationary species identities. Shared training/held blocks and repeated years preclude treating these 60 cells or species as independent. Directional recurrence does not establish that the same individual taxa remained increasing/decreasing across all years, or that prospective timing of gate passage was predicted without current outcomes: gate evaluation uses held response. 

Magnitude prediction **fails overall**: state ridge versus history ridge gain -9.2889%; state ridge vs last-response improvement 23.9824%, yet vs zero -30.4457%. The zero forecast outperforms every nonzero model in aggregate SSE, despite some ages having positive signed direction prediction. State covariates improve history-ridge SSE in only 2/10 ages (10 and 12), worsen remaining eight. Direction and magnitude are distinct objectives; cannot claim added control-state information explains timing or improves overall forecasting. Neither aggregate strong-direction gate nor final age passes. No further tuning on the held ages. 

| Age | State gain over history | Prior-year sign mean min | State sign mean min | Concurrent sign mean min |
|---:|---:|---:|---:|---:|
| 6 | -9.0771% | -22.5192 | -44.4258 | -35.6342 |
| 7 | -13.3006% | -25.0717 | -51.4642 | -23.7008 |
| 8 | -8.2793% | -35.2792 | -36.4050 | -20.1717 |
| 9 | -7.1516% | 22.6517 | 7.2367 | 15.3892 |
| 10 | 11.3972% | 25.6158 | 28.1067 | 25.5600 |
| 11 | -49.5366% | -6.7058 | -0.2858 | -11.9708 |
| 12 | 3.2623% | 19.8233 | 24.8842 | 17.9558 |
| 13 | -30.6105% | 9.4908 | 8.4683 | 0.6292 |
| 14 | -2.6287% | 12.3767 | 13.1017 | 12.5117 |
| 15 | -14.4610% | -8.3258 | -10.4667 | -19.8292 |

All results at validation/ahif_v1/kansas2020_time_block_forecast_v1.json. Lagged controls describe the untreated reference after treatment start, not pretreatment state of seeded plots. No seed membership or independently measured recruitment opportunity, no causal compensation, equivalence, climate trend, new raw source audit, or formal meta-analysis. Same curated source; existing AHIF v0.9 counts unchanged. Freeze Kansas predictive follow-up. Next priority obtain original factorial recruitment-source dataset or original species-by-plot biomass/seed identity and assignment metadata, not further optimizing this one selected series. Main and PRs3–5 unchanged.