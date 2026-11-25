# Portal observational allocation pilot: verified results

2026-11-25 JST. Successful source-locked run https://github.com/zuizui0223/knk/actions/runs/45569020696, job 139050795088, commit c201b31ff22a770004c9137c0f5b21c942e7d5ea. Artifact 15000454275 sha256:ae96bbe6ece860851873dbf7d10c0ef01095aeadf115b2862ff64eaf814b2529. Five original hashes and synthetic checks passed. 79,828 raw capture rows, 10 author-confirmed constant control plots, 112 moons in scope, 73 complete full-effort QC periods, 21 target species codes, 6,477 eligible rows, 1 ambiguous individual ID (2 rows) and 2 duplicate capture rows removed, 6,473 counted unique plot-period-ID rows. Eight nonoverlapping four-period forecast windows (older training windows may overlap), four overlapping plot splits × two orientations; one ecological site. Captures, not detection-corrected population sizes.

| Split | Orientation | M2 | V | H | AH | VH | HH | M2−AH pp | V−VH pp | H−HH pp |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0 | .430614 | .450587 | .387025 | .417486 | .423208 | .387927 | 1.3129 | 2.7379 | -.0903 |
| 0 | 1 | .339584 | .359400 | .287374 | .325736 | .330899 | .282805 | 1.3849 | 2.8501 | .4568 |
| 1 | 0 | .318685 | .329984 | .303346 | .331177 | .335081 | .319184 | -1.2492 | -.5098 | -1.5838 |
| 1 | 1 | .453540 | .482899 | .433825 | .431118 | .442905 | .403597 | 2.2422 | 3.9994 | 3.0228 |
| 2 | 0 | .403646 | .422665 | .376948 | .355568 | .353241 | .347978 | 4.8078 | 6.9424 | 2.8970 |
| 2 | 1 | .365846 | .386856 | .341474 | .381808 | .378412 | .359096 | -1.5962 | .8445 | -1.7623 |
| 3 | 0 | .432628 | .474632 | .380442 | .437839 | .455404 | .403096 | -.5210 | 1.9227 | -2.2654 |
| 3 | 1 | .388458 | .419048 | .377097 | .337427 | .347855 | .343384 | 5.1031 | 7.1193 | 3.3714 |

Mean across eight dependent orientations: M2−AH +1.436 pp (5/8 positive), V−VH +3.238 pp (7/8 positive), H−HH +.506 pp (4/8 positive). Primary H−HH effectively near zero; no reliable incremental recency gain. H beats M2 in 8/8 dependent orientations and V in 8/8, descriptive; a long-run response distribution is useful within this site. Primary orientation dominant overlap .375, lag overlap .375, cancellation .466, signed gap .1375; swapped .6875/.4375/.502/.1375. No ecological-memory claim, no independent time-series inference, no p values, no causal assembly effect or meta-analysis. H_oldpair_mix reuses a pair included in H and thus is NOT an equal-information comparator. New paired comparison must use the same baseline excluding both candidate periods. Keep original outcome frozen. Do not select a favourable split. Full BioTIME still source STOP; teaching subset zero eligible windows. No manuscript/AHIF changes.
