# Matched replication-scale results and directional overlap decomposition — 2026-11-25

Verified run https://github.com/zuizui0223/knk/actions/runs/45589332743, job 139119047796, source-locked scorer 86a75c6a71f13cea51cd34e8a7f6cbf689130998. Artifact 15009539604 SHA256 ad47b487ae798b0b98d189fd7329574228824acae199bf4d9f629ceb4b541bc6. Original source hashes, original five-source one-lag replay, author Hayashi complete-well export and synthetic identities passed. Seven existing study families, 43 analysis panels (Saine ten assignment strata, not the original five site-average panels), 53 eligible source-size configurations, 72,192 dependent fold/date cells per common 2+2 comparison across sources. Every reported config has full common support. No new publication, inferential p-value, sampling SE or formal meta-analysis. Saine secondary site-average results kept separately. Fixed shared-control partitions and paired block membership; physical count does not guarantee equivalent independent replication across designs. Kansas original linkage conditional; Hayashi density-position alias retained.

## Primary matched 2+2 results

All error gains in percentage points. G = TV(signblind) minus TV(hard signed forecast). H = expected TV of independent fair random hard signs at the same predicted absolute allocation minus TV(hard signed forecast). Concurrent uses train-side target-date data and is not a forecast. Fixed-mass concurrent advantage changes sign only, using the lag absolute allocation throughout.

| Source | G lag | H lag | G concurrent | H concurrent | Fixed-mass concurrent advantage | Increase training 1 to 2, held=2 (G) | Increase held 1 to 2, train=2 (G) |
|---|---:|---:|---:|---:|---:|---:|---:|
| pond | 8.5408 | 16.0553 | 9.4965 | 17.2886 | .4201 | 3.5409 | 1.6177 |
| Kansas | -3.0314 | 4.1066 | -3.0844 | 4.5263 | -1.2320 | 1.5393 | -.7028 |
| POEM | -6.4436 | .3551 | -3.4937 | 5.0554 | 3.8259 | .3963 | .3467 |
| field | 19.8563 | 26.4973 | 22.5118 | 29.8231 | -.6829 | 2.6279 | 2.7029 |
| Mancuso | 17.6485 | 23.6205 | 21.9590 | 27.3134 | 1.0152 | 1.4070 | 1.1213 |
| Saine | -3.0116 | .6797 | -3.4029 | .7386 | .4052 | .2152 | -.6782 |
| Hayashi | 5.4838 | 12.7082 | 8.5844 | 16.1421 | 1.6656 | 2.7877 | 1.8716 |

Four positive G, three negative G; seven positive H point estimates (POEM and Saine close to zero). H is a fair hard-sign reference expectation, not a treatment-assignment null distribution or significance test. Same mass uniform sign splitting is a convex-loss hedge, not the expected score of a hard random sign. The difference between H and G is mathematically required to be nonnegative, not a biological effect.

## Exact decomposition derived from frozen summary, not another fitted model

For each cell, p_i and q_i are normalized absolute train-lag and held-target response masses. Both signed profiles place each taxon's positive mass in one of two sign coordinates. O=sum_i min(p_i,q_i)=1−TV(p,q). M=sum_i min(p_i,q_i) I(sign agrees). Disjoint/zero-support taxa contribute zero to both. L_signed=1−M, L_flipped=1−(O−M), expected fair-hard loss=(L_signed+L_flipped)/2=1−O/2. Consequently H=M−O/2. Define A=M/O when O>0, then H=O(A−1/2). Also J=E_random_hard TV−TV(signblind)>=0 by convexity, and G=H−J. A ratio computed from aggregated H and O, A*=1/2+mean(H)/mean(O), is the overlap-mass-weighted agreement over the fixed study averaging scheme, NOT an unweighted average of taxon sign matches or a gross held response fraction. Absolute allocations and read/biomass lanes remain distinct. Derived after scoring, no new independent evidence.

| Source | Allocation overlap O (%) | Overlap-weighted agreement A* (%) | Hedging advantage J (pp) |
|---|---:|---:|---:|
| field | 64.46 | 91.11 | 6.64 |
| Mancuso | 66.96 | 85.28 | 5.97 |
| pond | 55.60 | 78.87 | 7.51 |
| Hayashi | 63.24 | 70.10 | 7.22 |
| Kansas | 49.59 | 58.28 | 7.14 |
| Saine | 20.41 | 53.33 | 3.69 |
| POEM | 33.95 | 51.05 | 6.80 |

Low overlap and weak directional agreement jointly limit Saine and POEM; Kansas has more overlap and weak positive agreement. POEM concurrent H +5.06 vs lag +.36 and fixed-mass concurrent advantage +3.83 pp suggest temporal mismatch may contribute, but concurrent G is still negative. Saine concurrent H +.74 vs lag +.68: no meaningful descriptive temporal gap, poor concurrent reproducibility as well. These are diagnostic differences, not identified mechanisms or proof that a community reorganization occurred. Concurrent absolute allocation differs from lag allocation except in the explicit fixed-mass diagnostic. The sign advantage in field is slightly negative at fixed lag mass despite positive concurrent G increment; do not equate concurrent-minus-lag G with pure directional temporal loss.

## Replication-scale curves

Hayashi (n+n), n=1/2/3/4/8/16/32: G +1.71/+5.48/+8.47/+10.73/+15.45/+19.31/+22.60 pp; A* 63.74/70.10/74.64/77.66/83.22/87.53/90.39%; O 63.28/63.24/64.54/65.38/68.20/70.28/72.63%. Both train and held increases improve G in this system: 32+2 +12.47 pp, 2+32 +10.13 pp, 32+32 +22.60 pp. Same wells and group contrasts, not independent experiments or a causal replication effect. Saine stratum-level (n=1/2/3/4), G -2.52/-3.01/-3.25/-3.42 pp, A* 52.55/53.33/53.95/54.19%, O 15.80/20.41/23.14/25.06%. More overlap but little directional agreement; G becomes more negative because sign-blind hedging benefit also changes. Full-size range limited by smallest arm. Site-average secondary Saine 2+2 per stratum G -3.34 pp, H +.89 pp; not the same total sample size as the primary stratum-level comparison.

Exact JSON summary validation/ahif_v1/replication_scale_summary_20261125.json; 53-row separate source-size registry database/ahif_mechanism_v0/replication_scale_effects_20261125.csv. Original v0.2 7-source study rows and 38 primary panels remain unchanged. All 53 rows are meta-ineligible: no calibrated SE, dependent folds/configurations, differing time intervals/measurement units and incomplete study independence. Do not use configurations as independent study effects. No manuscript change. Next priority recruit genuinely new publications with verified repeated treatment/control abundance data and source-level holdouts; avoid optimizing the same seven families further.
