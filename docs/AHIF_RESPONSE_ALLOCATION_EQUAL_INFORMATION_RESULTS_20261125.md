# AHIF matched-information recency — verified retrospective result

2026-11-25. Run https://github.com/zuizui0223/knk/actions/runs/45572257570 succeeded; job 139061521350; commit d2e9e42fcabf6909144d39a20f98e8b72d0a9ae3. Artifact 15002342538 sha256:7e6715b1cdc5f1d592ae6355cc7328047ba72a7002f71e36e997d21447af8648. Original source rechecks and prior contrast numeric replay passed across all 19 panels. Matched recency estimand: 13 eligible panels, 1,960 dependent fold-time cells, 13,820 matched candidate contrasts. POEM 6 panels have three dates, structurally ineligible for this t>=3 gate. These counts are not independent experiments.

| Source | Matched older TV | Matched recent TV | Gain pp | H vs M2 gain pp | H vs V gain pp | Current dominant top1 overlap | Lag response top1 overlap |
|---|---:|---:|---:|---:|---:|---:|---:|
| Narwani pond | .572342 | .456998 | 11.5344 | 10.6465 | -1.1469 | .412616 | .592014 |
| Kansas | .475704 | .447715 | 2.7989 | 8.6249 | 1.8849 | .277778 | .277778 |
| van Steijn field | .517838 | .376641 | 14.1197 | 13.3257 | 18.3542 | .000000 | .643750 |
| Equal-source mean | .521961 | .427118 | 9.4843 | 10.8657 | 6.3641 | .230131 | .504514 |
| Excluding Kansas | .545090 | .416820 | 12.8270 | 11.9861 | 8.6037 | .206308 | .617882 |

Positive gain means lower held TV; source and within-source weighting as prescore. H is mean historical absolute treatment-control contrast, M2 square root control abundance, V training variance proxy. L alone equal-source TV .411481, H .451235 (same complete cells); Kansas L is worse than H. Separate earlier Portal observation result H-vs-recent mixture near zero; Portal target was temporal change in captures, while AHIF target here is a contemporaneous treatment-control contrast at a future sampled date. These are different estimands, so no demonstrated observation-vs-experiment effect modification. No time-varying causal interaction, new independent dataset, formal meta-analysis, or inference from overlapping folds. Recent candidate may predict merely persistent treatment effects or environmental context. Absolute contrast does not identify sign. Next diagnostic compares treatment contrast with matched latest state/noise and evaluates signed response persistence, without choosing favorable sources. No edits to frozen AHIF v1.0, prior PR6 or manuscript.
