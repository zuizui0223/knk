# Seed-addition synthesis: held-block two-direction validation

- Source: Ladouceur et al. (2020), DOI 10.1111/ele.13566; Figshare 12319682 v1 (two verified source files).
- Frozen diagnostic code: 0b0aeed347a2dff52249496dda9b8269674b1ef1
- Successful run: https://github.com/zuizui0223/knk/actions/runs/45539576248; job 138954498034.
- Same curated source as the previous analysis, not new raw studies.
- Hold out block groups; learn increasing and decreasing sets on training blocks, then evaluate their signed mass changes on held-out blocks. Strong descriptive gate requires mean held increase > 0, mean held decrease > 0, and mean split-wise minimum > 0. Repeated splits are not independent; neither a p-value nor a confidence interval. Gate devised after the initial opposition analysis; not a preregistered confirmation.

## Counts

| Mode | Scored experiment-age cells | Eligible | Both mean directions positive | Strong gate | Final eligible | Final strong |
|---|---:|---:|---:|---:|---:|---:|
| direct | 24 | 23 | 8 | 5 | 3 | 0 |
| indirect | 34 | 30 | 20 | 20 | 5 | 4 |

## Direct measurement: passing cells

| Experiment label | Source age code | Held increase | Held decrease magnitude | Mean minimum | Fraction both positive |
|---|---:|---:|---:|---:|---:|
| California_Invade | 4 | 18.195 | 29.247 | 13.815 | 0.7734 |
| Kansas_KUFS_LTER_Hay_Meadow_Exp_2 | 9 | 15.676 | 48.359 | 15.676 | 1.0000 |
| Kansas_KUFS_LTER_Hay_Meadow_Exp_2 | 10 | 26.588 | 39.229 | 25.534 | 1.0000 |
| Kansas_KUFS_LTER_Hay_Meadow_Exp_2 | 12 | 45.009 | 15.712 | 15.712 | 0.8145 |
| Kansas_KUFS_LTER_Hay_Meadow_Exp_2 | 14 | 28.359 | 23.866 | 11.693 | 0.8145 |

Direct cells passing the diagnostic are confined to two source experiments; none of the three eligible final direct cells passes. Indirect allocation of measured total biomass by cover is not independent support for opposing absolute biomass responses. Source age codes are not independently verified durations. The synthesis lacks plot-specific seed membership; added vs resident species cannot be identified. No causal recruitment compensation, proof of total biomass equivalence, independent meta-analysis, or new original-source audits are claimed. DB v0.9 primary effect and cancellation ledgers remain unchanged; preserve main and PRs #3–5.
