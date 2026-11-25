# BioTIME: source availability and stationary-noise diagnostic

2026-11-25. Separate exploratory macro branch; no change to frozen KNK manuscript or AHIF effect ledger.

## Observed source gate
- Full v2.0 release DOI 10.5281/zenodo.15222193: HTTP 403 before raw RDS acquisition, run 45565939719. No full-release result.
- Official BioTIMEr teaching subset commit 7c33f32f0235d2bb396ea8b71f7f39a91f0253f9: both source Git blob hashes verified; SHA256 data 3a84d0fd3b7ec1eac360a14675e3c9e0e294e353fd71fc80fab1a3912e537af4 and metadata b1629cf0c148ac3523fc59d7078507359e0d1b03e7d2fdb7c6a04eb54970e7fd.
- 81,084 rows, 12 metadata records. Original source-level bird exclusion: 333,372,41,59; density and MeanCount excluded by unchanged Count-only gate. Remaining metadata-eligible nonavian Count study IDs: 10,18,211,402,431,466 (43,399 raw rows before strict sample audit; row counts summed from source log).
- Full synthetic fixture: passed; exact and grid each two nonoverlapping windows, cancellation and ranked-TV checks passed.
- Actual source scoring STOP: no eligible four-year windows with >=2 species under unchanged sample rules. Not evidence for or against a biological hypothesis. Source-audit artifact 15000851081, run https://github.com/zuizui0223/knk/actions/runs/45567795227, job 139046748069.
- Teaching subset is not globally representative and is not a second independent BioTIME study source. Do not relax the eligibility rules on this source post hoc or pool it into AHIF. Full release still unavailable.

## Additional synthetic observation-only diagnostic
Executed locally with NumPy and seed 20261125; script benchmarks/simulate_biotime_stationary_noise.py. 4,000 independent synthetic communities per noise model, 30 taxa each, four years, two independent halves, two independent events per half. True species rates fixed through time; no biological temporal change, temporal memory, interaction, or succession. Poisson or independent gamma-Poisson counts (constant common dispersion 0.6 or fixed heterogeneous species dispersion). Same frozen M1-M4 definitions. TV lower is better; reported improvement is 100*(TV_M2-TV_M4) in percentage points, MC interval reflects only simulated-community sampling, not empirical uncertainty.

| Noise process | M1 TV | M2 TV | M3 TV | M4 TV | M4 improvement over M2 (pp) | Monte Carlo 95% interval (pp) |
|---|---:|---:|---:|---:|---:|---:|
| Poisson | .383644 | .301899 | .353216 | .332405 | -3.050577 | [-3.164975,-2.936178] |
| common overdispersion | .328172 | .371788 | .350683 | .360602 | +1.118651 | [.910136,1.327166] |
| species-specific overdispersion | .410284 | .434358 | .411100 | .414114 | +2.024410 | [1.782297,2.266522] |

M4 beats M2 in the latter two cases with no biological temporal memory. In those same cases M3 does not beat abundance-only M1. Approximate oracle sqrt-variance TV respectively .299536, .303655, .366552 (uses true parameters unavailable in observed data). This is a diagnostic counterexample, not a calibrated empirical null or a correction to earlier results. Temporal and observation splits alone do not rule out persistent abundance/dispersion differences. Before interpreting any full-release result, separately preregister and evaluate matched stationary-count simulations or within-site event permutations preserving sampling effort, time structure constraints, and species-specific dispersion, plus compare against training-only variance-aware baselines. Do not infer biological memory solely from M4 > M2. All actual ecological effects remain unscored.
