# Fixed-role response persistence: dominant identity vs distributed support
2026-11-25 JST. Source-pinned run https://github.com/zuizui0223/knk/actions/runs/45559547234, job 139019671701, scorer/prescore commit 7ba6523e4848da62f3e13eb03a04abe28b66458a; synthetic tests, input hash/grid/mapping checks, every archived pooled role contrast and sign flag reproduced, full real scoring and artifact upload passed. See validation/ahif_v1/fixed_role_contributor_summary_v1.json; full per-time/split results in artifact 14996848280, sha256:4d11a83b6d301533565a3ec53a8ab75591638e6c21c119c294a000b0ce38522f, 660008 bytes. Same already scored Narwani and Kansas experiments, no new independent source or effect. All dates and both signs retained. Source references: Narwani 2019 DOI 10.1098/rspb.2019.1857 / 10.5061/dryad.qv9s4mw99 (Zenodo 4941856); Ladouceur 2020 DOI 10.1111/ele.13566 / 10.6084/m9.figshare.12319682.v1; Foster 2009 DOI 10.1890/08-0849.1. Kansas original cited EDI still unavailable and membership link conditional. Narwani treatment manipulates foundation species, not arrival history.
## Directional gross support concentration
For each fixed group, orient the signed difference by expected direction (+ first, − second). Define q=max(oriented species/morphotype contrast,0). Maximum and top-three shares are relative to positive gross support, not net effect, abundance, recruitment or a causal attributable fraction. Average over defined held cells; zero gross leaves share undefined. Synechococcales is an order-level morphotype label, not an identified individual species.
| Contrast | Role | Held gross defined | Top1 share | Top3 share | Most frequent top contributor (count) | Mean signed net |
|---|---|---:|---:|---:|---|---:|
| pond_MD-C | positive first | 340 | 99.74% | 100.00% | Synechococcales (338) | 4.524 |
| pond_MD-C | negative second | 360 | 76.63% | 98.30% | Coelastrum_2 (231) | 11.862 |
| pond_M-C | positive first | 236 | 98.45% | 100.00% | Synechococcales (194) | 0.092 |
| pond_M-C | negative second | 360 | 72.88% | 97.96% | Coelastrum_2 (234) | 14.500 |
| pond_D-C | positive first | 179 | 97.95% | 100.00% | Synechococcales (105) | -0.242 |
| pond_D-C | negative second | 357 | 78.72% | 98.07% | Coelastrum_2 (273) | 12.965 |
| pond_C-O | positive first | 360 | 96.41% | 100.00% | Synechococcales (360) | 0.510 |
| pond_C-O | negative second | 309 | 82.52% | 98.58% | Unknown_9 (131) | -19.740 |
| kansas | positive first | 60 | 43.07% | 78.60% | Sorghastrum.nutans (33) | 53.125 |
| kansas | negative second | 60 | 54.27% | 87.80% | Schedonorus.arundinaceus (33) | 34.867 |
## Excluding training-selected contributors from independent held contrasts
Top up-to-three positive supporting taxa selected in each role using training-only data. Arithmetic exclusion from an observed held vector, not actual removal, reassembly, intervention, or forecast of an ecological extinction. Conditional sign retention denominator = held cells with both expected signs before exclusion. Also preserve unconditional residual-opposition counts in JSON; conditional counts are not forecast accuracy.
| Contrast | Full opposition / all cells | Early top3 retained | Previous top3 retained | Expanding-history top3 retained | Same-time other-half top3 retained |
|---|---:|---:|---:|---:|---:|
| pond_MD-C | 280/360 | 95/280 (33.93%) | 54/280 (19.29%) | 58/280 (20.71%) | 44/280 (15.71%) |
| pond_M-C | 183/360 | 61/183 (33.33%) | 81/183 (44.26%) | 51/183 (27.87%) | 80/183 (43.72%) |
| pond_D-C | 52/360 | 38/52 (73.08%) | 36/52 (69.23%) | 26/52 (50.00%) | 38/52 (73.08%) |
| pond_C-O | 0/360 | NA (0 full opposition) | NA (0 full opposition) | NA (0 full opposition) | NA (0 full opposition) |
| kansas | 51/60 | 47/51 (92.16%) | 48/51 (94.12%) | 44/51 (86.27%) | 47/51 (92.16%) |
Primary MD-C: previous top3 in both groups retains 54/280 (19.29%), same-time other-half 44/280 (15.71%); Kansas 48/51 (94.12%) and 47/51 (92.16%). Expanding top3 retains 58/280 vs 44/51. Early top3 retains 95/280 vs 47/51. Selection up to three, not necessarily exactly three (see random expected capture in JSON); group sizes differ (pond 8/43, Kansas 27/62) and taxa are different grain. No formal between-system statistical comparison. For pond cyan, previous top3 captures 90.95% gross support on average, mean net 4.5239 -> .16457 mm^3/L; same-time top3 .04993. Kansas sown mean net 53.125 -> 27.6359 g/m², other named decrease oriented 34.86675 -> 21.213. Do not pool absolute magnitudes.
## Contributor identity overlap and temporal vs spatial disagreement
Normalized positive support overlap=sum min of profiles, 0 disjoint, 1 identical. Comparison on common eligible cells only (positive gross support in lag train, concurrent train and held current). A low lag overlap can result from block/pond heterogeneity, not only time turnover.
| Contrast | Role | Common cells | Lag training vs held current | Same-time training vs held current | Same-time minus lag |
|---|---|---:|---:|---:|---:|
| pond_MD-C | positive first | 304 | 0.981 | 0.985 | 0.004 |
| pond_MD-C | negative second | 360 | 0.433 | 0.459 | 0.026 |
| pond_M-C | positive first | 90 | 0.458 | 0.476 | 0.018 |
| pond_M-C | negative second | 360 | 0.572 | 0.616 | 0.044 |
| pond_D-C | positive first | 28 | 0.328 | 0.214 | -0.114 |
| pond_D-C | negative second | 351 | 0.547 | 0.580 | 0.034 |
| pond_C-O | positive first | 360 | 0.943 | 0.934 | -0.009 |
| pond_C-O | negative second | 220 | 0.237 | 0.390 | 0.153 |
| kansas | positive first | 60 | 0.460 | 0.502 | 0.042 |
| kansas | negative second | 60 | 0.321 | 0.318 | -0.003 |
MD-C cyan support dominated by Synechococcales in 338/340 defined cells, 99.74% mean largest share, both overlaps ~.98. The other algal group shows lower overlap and multiple major contributors; no causal turnover demonstrated. Kansas support has low overlaps but same-time between-block overlaps are also low (.502 and .318 vs lag .460 and .321). Thus distinguish distributed/heterogeneous support and response robustness from proof of temporal replacement. Prior-known windows are outcome-informed descriptive only.
## Conclusion
Two different forms of group-response persistence: dominant-label concentration in pond cyan response versus support distributed across multiple taxa in Kansas; shared role opposition need not imply interchangeable species. Stronger exclusion of trained dominant contributors erases much of the pond response while Kansas signed group differences generally remain. This does NOT identify functional redundancy, replacement mechanisms, population persistence, actual manipulation of contributor removal, or causal compensation. Six/36 overlapping splits and repeated dates dependent, all descriptive no p-values. Differences in group sizes, measurement (biovolume vs dry mass), intervention, and taxonomic resolution preclude general or causal cross-system inference. All four pond contrasts retained, including C-O with no full opposition in scored cells. Freeze follow-up; independent actual randomized arrival/recruitment or perturbation dataset needed for mechanism. Existing AHIF DB v1.0 and 161 effects, main and frozen PRs3–5 unchanged.
