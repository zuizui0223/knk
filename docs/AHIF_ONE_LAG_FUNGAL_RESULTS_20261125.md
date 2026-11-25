# One-lag sign forecast and fungal introduction — results 2026-11-25

Verified run https://github.com/zuizui0223/knk/actions/runs/45586981989, job 139110965349, source-locked scorer f7c3035bdd87172706d81a7affe4cb40a573c7c6. Artifact 15008623139 SHA256 477f01b4a65e048068e539eb9476339e81498c74692098d6a58023079d702d1c. All pinned source checks, original prior signed-mixture scores, original fungal metadata coverage and synthetic tests passed. Exact source summaries and all directional endpoints: validation/ahif_v1/one_lag_fungal_summary_20261125.json. No original effects ledger, prior results or manuscript modified.

## Harmonized adjacent-date forecast

G = TV(signblind) − TV(signed) with identical predicted absolute mass, both using only the immediately preceding training-unit contrast. Equal fold/date/panel weighting, original five-source design, 27 panels. The primary target dates remain those of previous PR8/9 (t>=2); all-adjacent t>=1 sensitivity adds the earliest forecastable date. Values are pp, not percentages of accuracy. All cells defined.

| Study | Prior target dates cells | G pp | All adjacent G pp | Same-direction gross primary |
|---|---:|---:|---:|---:|
| Narwani | 1872 | +8.6065 | +7.5087 | .678845 |
| Kansas | 78 | −2.9590 | −2.1588 | .552362 |
| POEM | 1200 | −7.0608 | −3.8301 | .506058 |
| van Steijn field | 320 | +21.1337 | +17.3947 | .841092 |
| Mancuso E1 | 1152 | +17.9880 | +14.4425 | .799933 |

Five-source sign pattern unchanged. Do not conflate original old+recent mixture estimand with one-lag G. All cells remain dependent. Kansas source linkage conditional. Mancuso compositional measurement lane.

## Saine fungal log experiment

Original 460 assigned physical logs, 446 complete post-QC triplets, 347 introduced and 99 control, five sites, ten allocation strata, 4653 non-target OTUs. Original RData/author script hashes, complete case and prior source coverage replayed. Relative read profiles only. Within-site: 256 seeded disjoint half partitions × 2 orientations × 5 sites (2560 dependent split-site cells). Half sizes from np.array_split (ceil/floor, both orientations). Binary introduction mixture is the primary contrast, balanced seven-identity sensitivity excludes two identities with singleton strata, uses 371 logs and equal label weights. This sensitivity changes estimand.

| Comparison | Within-site primary G pp | Held-site G pp | Balanced7 within G pp |
|---|---:|---:|---:|
| 2020 response vs 2021 | −2.4635 | −2.1969 | −2.6064 |
| 2019 pre-label contrast vs 2021 | −3.4053 | −4.1365 | −3.5008 |
| 2021 other-unit concurrent benchmark | −2.9174 | −.1892 | −3.2802 |
| Baseline-adjusted 2020 vs 2021 | −5.9317 | −2.3252 | −5.7840 |
| 2020 vs 2019 sign at fixed 2020 mass | +1.8786 | +1.6635 | +1.6929 |

Primary one-lag signed TV .794013, signblind .769378, absolute allocation TV .651650. Same-direction gross .468865, opposite .390605, abstain .140530. Within-site primary G negative all 5 sites (Kesijarvi −1.4991, Lapinjarvi −1.9558, Luukki −3.6971, Saajarvi −2.7173, Seitseminen −2.4482 pp). Held-site 1/5 positive (Saajarvi +.0736 pp), remaining negative. Baseline-adjusted within-site G all negative. Importantly concurrent signed forecast also negative all 5 sites, so weak directional transfer does NOT establish faster temporal reversal or ecological memory loss; low contrast strength, among-log heterogeneity and compositional observation can explain it. Post-treatment signs modestly outperform pre-treatment signs at fixed magnitude, but do not beat sign-blind forecast. Negative mean G is descriptive, not proof of a significant adverse effect or absence of directional persistence in all fungi.

Six existing publication-level sources now comparable in one-lag endpoint (three positive, three negative), but measurement lanes, lab dependence, QC and study selection differ. Saine already in AHIF source ledger; no newly recruited publication. No formal meta-analysis, inverse variance weighting, causal interaction claim or statistical inference from repeated folds. Fixed-split assignment null preregistered separately before these results and not yet used for this interpretation.
