# Saine 2025 — site-held-out partner identity transfer, verified negative result

2026-11-24. Source 10.5281/zenodo.13819389, pinned RData MD5 4dbc27f81ffdcf1925887c7e0805b031, run https://github.com/zuizui0223/knk/actions/runs/45388257428, job 138424791903, commit 43a837e9f5d9e02d12f3c410c547e52da6388a01. All source joins and four synthetic tests passed. 446 complete logs (347 treated, 99 controls) from the same 460-log experiment; not independent new publication evidence.

| Response | Pooled gain vs zero | Nested identity gain vs pooled | Full identity gain vs pooled | Conditional p | Holm p |
|---|---:|---:|---:|---:|---:|
| named | -0.041132 | 0.000000 | -0.162243 | 0.984 | 1 |
| rank | -0.068633 | -0.007761 | -0.110137 | 0.712 | 1 |

Outer folds, no held-out site used in shrinkage selection:

| Representation | Held site | Alpha | Zero SSE | Pooled SSE | Nested identity SSE | Full identity SSE |
|---|---|---:|---:|---:|---:|---:|
| named | Kesijarvi | 0 | 3.31145 | 3.48327 | 3.48327 | 3.84324 |
| named | Lapinjarvi | 0 | 2.03417 | 2.35886 | 2.35886 | 2.83560 |
| named | Luukki | 0 | 2.70764 | 2.64554 | 2.64554 | 3.06350 |
| named | Saajarvi | 0 | 1.96383 | 2.05510 | 2.05510 | 2.42948 |
| named | Seitseminen | 0 | 2.74624 | 2.74553 | 2.74553 | 3.27243 |
| rank | Kesijarvi | 0 | 0.64596 | 0.61054 | 0.61054 | 0.68705 |
| rank | Lapinjarvi | 0.25 | 0.35898 | 0.41084 | 0.43940 | 0.60774 |
| rank | Luukki | 0 | 0.93174 | 0.91453 | 0.91453 | 1.08014 |
| rank | Saajarvi | 0 | 0.47782 | 0.43302 | 0.43302 | 0.35251 |
| rank | Seitseminen | 0 | 1.02931 | 1.31125 | 1.31125 | 1.35806 |

Neither the pooled response nor introduced-identity-specific response generalizes under the prescored year-1 training → year-2 held-site group-response squared loss. Named alpha=0 across all sites; rank alpha=.25 for Lapinjarvi and 0 otherwise; rank nested prediction is slightly worse than pooled. All nulls preserved. The pooled-vs-zero gain is descriptive, not a tested null; p-values test identity gain only. The negative predictive result combines space and time transfer, cannot separate site dependence, temporal evolution, and sampling variance. No inference of absence of interactions. No proof of buffering. The small binary named centroid effect at year two remains source-verified and depth-robust, but does not imply useful prediction of group-response vectors. No independent mediator or fitted alpha competition coefficient.

Future narrow question: does pre-introduction community state condition treatment responses across sites? Test baseline richness/dominance/read-depth features and natural/felled context with a held-site randomized treatment contrast loss. Nulls accepted; no more claims from separate significance levels. Same source experiment; no independent study count increment, 161 DB effects remain unchanged. No merge or main/manuscript changes. Source aggregate validation/ahif_v1/saine2025_site_transfer_v1.json.
