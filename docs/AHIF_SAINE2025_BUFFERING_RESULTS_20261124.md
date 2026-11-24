# Saine: paired structural-buffering falsification and sequencing-depth sensitivity

2026-11-24. Source-verified run https://github.com/zuizui0223/knk/actions/runs/45387588766, job 138422482690, commit ca7b6858fc9c8f44b913b54af98e88752abe5ff1; original file checks and four new synthetic tests passed. Same 460-log publication family, not independent new study. Prescore contract `docs/AHIF_SAINE2025_BUFFERING_PRESCORE_20261124.md`; aggregate `validation/ahif_v1/saine2025_buffering_depth_v1.json`.

## Primary paired displacement

| Year | Endpoint | Treated minus control TV | Unadjusted conditional bootstrap 95% interval | raw p | Holm p (6 tests) |
|---|---|---:|---|---:|---:|
| 2020 | named_TV | 0.031736 | 0.00279, 0.06168 | 0.0183 | 0.1098 |
| 2020 | rank_TV | -0.007252 | -0.03855, 0.02458 | 0.6748 | 1.0000 |
| 2020 | identity_excess_TV | 0.038988 | 0.00473, 0.07335 | 0.0326 | 0.1630 |
| 2021 | named_TV | 0.022579 | -0.00082, 0.04565 | 0.0363 | 0.1630 |
| 2021 | rank_TV | 0.014714 | -0.01615, 0.04663 | 0.4103 | 1.0000 |
| 2021 | identity_excess_TV | 0.007866 | -0.03024, 0.04472 | 0.6738 | 1.0000 |

No endpoint survives six-test correction. In year one, identity-excess TV difference +.038988 (raw p=.0326, Holm .163); year two +.007866 (Holm 1). These do NOT establish structural buffering or equivalence. TV gap is mechanically nonnegative within each log; contrast between treatments is tested. Different significance of earlier named vs rank R2 cannot establish a difference either: earlier E named .002368, rank .002837. TV magnitude versus centroid response are distinct estimands.

## Fixed non-target depth

| Depth | Complete logs | Eligibility min Holm p | 2021 named E | named Holm p | 2021 rank E | rank Holm p | 2020 identity-excess TV contrast | gap Holm p |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1000 | 421 | 0.314 | 0.002372 | 0.0120 | 0.004290 | 0.3195 | 0.047797 | 0.0624 |
| 5000 | 380 | 0.343 | 0.002328 | 0.0240 | 0.002468 | 0.4815 | 0.046693 | 0.0924 |

At 1000 and 5000 reads, the binary two-year named centroid response persists (E=.002372/.002328, within-family Holm .012/.024). Rank signal remains undetected; E is not demonstrably smaller. Matched unrarefied subsets show closely similar effects; no significant paired TV result even after rarefaction. Eight rarefactions averaged before inference, not eight replicates. 25 and 66 of the 446 QC logs excluded, respectively. Treatment-dependent eligibility tests not significant, but do not demonstrate exchangeability after post-treatment filtering. Sequencing depth sensitivity is not absolute abundance standardization or independent replication.

## Interpretation and next gate

This is a boundary condition, not proof of structural buffering. Small taxonomic centroid shift survives depth matching. Need test whether introduced-species-specific first-year response vectors from other forests predict second-year responses better than pooled introduction responses and zero effect. Separate ecological community read responses from independently measured interaction coefficients. Same source experiment and five sites; no causal mediator, no random-effects meta-analysis. Archive all outcomes without adding repeated sensitivity endpoints to the 161-row ledger. No frozen main/manuscript changes.
