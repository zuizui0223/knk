# Saine one-lag conditional sharp-null result

Verified run https://github.com/zuizui0223/knk/actions/runs/45587192721, job 139111682196, commit 9342ee84d1686044e9f186a30635f0bb98b71a05. Artifact 15008573588 sha256:a49067c1501c25fdf2a35f0a52958aff8cb79e3338769658ab4c3f353d35b3b7. Source checks, original assignment coverage and unit tests passed. All 20 fixed stratum-half cells had both treatment/control; 446 complete logs and 5 sites. No split redrawing. 999 permutations of labels within fixed stratum-half counts, same assignment across years. Exact summary validation/ahif_v1/saine_one_lag_sharp_null_summary_20261125.json; full null draws in artifact.

| Score | Fixed-split observed pp | Null median pp | Null 2.5–97.5 percentile pp | Conditional upper p | Holm p |
|---|---:|---:|---:|---:|---:|
| One-lag G | −3.4596 | −4.8959 | [−6.8085,−2.8934] | .086 | .172 |
| Lag vs pre signs at fixed 2020 mass | +.8909 | −.6144 | [−2.9733,+1.6243] | .100 | .172 |

Neither co-primary score rejects the conditional sharp null at .05. The first observed score is negative versus signblind, yet slightly better than the null median. The test preserves temporal covariance and cross-sectional measurement heterogeneity under label permutation; it is not a null of no directional memory, and conditioning on post-assignment complete-case QC makes exchangeability conditional. Failure to reject is not evidence of equivalence or absence of fungal introduction effects. The fixed partition score differs from the repeated-split descriptive estimate; both retained. No calibrated meta SE, formal pooling, or new independent publication. No original ledger/manuscript modification.