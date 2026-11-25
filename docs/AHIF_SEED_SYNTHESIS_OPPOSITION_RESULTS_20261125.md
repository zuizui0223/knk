# Ladouceur 2020: descriptive species opposition in curated seed-addition synthesis

2026-11-25 JST. Source-locked scoring run https://github.com/zuizui0223/knk/actions/runs/45539367672, job 138953805220, commit c22292ea96c561106abf122167a8f226a0b78b84. Both official v1 files match pinned MD5/SHA256; 3393 state rows, 58426 species rows, 12 experiment labels. Two synthetic tests and scorer/artifact upload passed. Original paper https://doi.org/10.1111/ele.13566 and data https://doi.org/10.6084/m9.figshare.12319682.v1. Secondary synthesis, not 12 independently audited raw source studies. Exact original author code/source methods audited previously. 
Two experiment labels STOPPED on physical-key/age ambiguity: Germany_Montane (116 duplicate rows and 40 treatment/dose-drifting keys) and Texas_Temple_Prarie (25 duplicate rows and 25 drifting keys). No arbitrary repair, no scoring from ambiguous units. Source Table 1 and data yr.trt range disagree in selected datasets (e.g. CCR_093 contains 0 and 21 beyond listed years, Jena2 3 vs table 2). Treat age as author table code pending original-source reconciliation, not verified absolute calendar duration. One Halle comparison lacks both treatment arms in matched strata. Final labels selected by design-only max matched source code, not outcome.
| Mode | Observed states | Matched states scored | Scored study-time comparisons | Crossfit available | Gap positive | Gap and signed gross positive |
|---|---:|---:|---:|---:|---:|---:|
| direct | 312 | 262 | 24 | 23 | 16 | 12 |
| indirect | 3081 | 2229 | 34 | 30 | 24 | 22 |
| Experiment | Mode | Final matched source age code | Matched blocks | Pooled C/G | Crossfit signed gap | Signed gross |
|---|---|---:|---:|---:|---:|---:|
| ASGA_Michigan | direct | 3 | 8 | 0.8243 | 1.5601 | -15.9590 |
| CCR_04 | indirect | 1 | 1 | 0.8319 | not scored | not scored |
| CCR_04 | direct | 2 | 1 | 0.9098 | not scored | not scored |
| CCR_093 | indirect | 21 | 30 | 0.9702 | 90.4402 | 82.7748 |
| California_Invade | direct | 5 | 8 | 0.5721 | -8.7387 | 20.1087 |
| California_Prop_Limi | indirect | 4 | 9 | 0.7607 | 11.3473 | 12.0175 |
| Halle | indirect | 3 | 10 | 0.8160 | 2.3648 | 3.5568 |
| Jena | indirect | 3 | 2 | 0.5208 | not scored | not scored |
| Jena2 | indirect | 3 | 6 | 0.8579 | 40.3383 | 38.2989 |
| Kansas_KUFS_LTER_Hay_Meadow_Exp_2 | direct | 15 | 4 | 0.7193 | 50.6086 | 2.8853 |
| Kansas_KUFS_LTER_Old_Field_Exp_1 | indirect | 14 | 6 | 0.3759 | 264.2295 | 873.6929 |

Original units g/m² as archived (some allocated from cover); species-sum contrasts not forced to agree with reported totals. Pooled gross-net gap can incorporate across-block cancellation; mean within-block metrics separate. Direct measurement results mixed. Positive signed gap alone can arise even with negative signed gross and/or unstable net sign; Michigan final positive gap but negative signed gross. Kansas Hay Meadow final positive gap but signed net negative (total direction unstable). Thus not evidence of robust simultaneous opposing species responses by itself. Many indirect records are constructed by allocating community mass to cover, so apparent compensation can be mechanical; source totals may be reused across time (CCR_093 identical reported total differences over many codes). Do not infer recruitment identity, ecological mediation, equivalence, causal opposition or formal meta result from these outputs. Need direct crossfit validation of both positive and negative learned species groups, without unstable total sign. All outcomes archived in validation/ahif_v1/ladouceur2020_species_opposition_v1.json; no independent source count or primary ledger change.
