# Narwani 2019 independent aquatic forecast results (2026-11-25 JST)

## Verified source and execution
- Paper: Narwani et al. (2019), DOI 10.1098/rspb.2019.1857; data DOI 10.5061/dryad.qv9s4mw99. Original public workbook from Zenodo 4941856, SHA256 6242f9c3a92cde7f451066ceebf6ae6c9cc6f924d5ed9bc0a7f11de7b9beb3c3.
- Source audit run 45556699625, prescore commit d26cbbe3d802c18696dfbd39f547d69d32feb4f0, forecast run 45557051796, job 139011521262: success including tests. Full derived forecast cells, scores and week-level summaries in artifact 14995009548 (sha256:f36c1bb71c430cb64200b0b665cdbaf64fdcbd59cf736fe7bc5ed394d719f023); this artifact expires according to Actions policy, source and scorer remain pinned.
- Twenty ponds, four per treatment, 660 raw pond-date rows, 17 microscopy times, 15 common complete times / 300 pond-date rows. Weeks 7 and 10 excluded (one incomplete biovolume row per week). 66 morphotypes excluding six published group totals; conservative 8 cyanobacteria / 43 other named algal morphotypes / 15 unresolved or other. All volume estimates expressed as mm^3/L, not dry mass or measured cell counts.
- Source metadata and treatment indicators disagree on the text coding of Nutrients; follow the original article and the source treatment labels/presence assignments. Spatial blocks are unbalanced, so use separate within-arm pond halves; no matched block interpretation. 36 overlapping 2+2 train/test selections per contrast, 10 forecast dates after 5 warmup observations; no p-values, no independent replication count of 360. Retain all ponds including 2C. 

## Response persistence, without forecasting
| Contrast | Cyan positive / other named negative (of 15) | Robust to unresolved assignment (of 15) | Source class totals, cyan positive / all five other negative (of 15) |
|---|---:|---:|---:|
| MD-C | 10 | 6 | 9 |
| M-C | 8 | 1 | 7 |
| D-C | 1 | 0 | 0 |
| C-O | 0 | 0 | 0 |

For MD-C, weeks 12,13,14,15,16,19 agree on opposition and remain robust to allocating every unresolved morphotype to either side. Weeks 24,28,32 are not robust under that rule. Published group totals corroborate opposition at 9/15 overall, including all six robust dates. The conservative mapping and source class totals are not identical: cyan subtotal agrees in 228/300 rows at absolute tolerance 1e-8, whole sum in 294/300. Keep these as distinct representations rather than force agreement. Pooled whole-arm signs and overlapping half-sample prevalences are descriptive only. The counts by contrast do not establish a factorial interaction or causal facilitation.

## Prospective held-pond and future-time forecasts
| Contrast | Best conservative-three forecast | SSE gain vs zero | Best original two-class forecast (simple baselines only) | SSE gain vs zero | Best morphotype forecast | SSE gain vs zero |
|---|---|---:|---|---:|---|---:|
| MD-C | cumulative | 10.17% | cumulative | 5.06% | cumulative | 1.67% |
| M-C | cumulative | 16.31% | cumulative | 15.35% | cumulative | 10.45% |
| D-C | cumulative | 10.62% | cumulative | 10.11% | cumulative | 11.49% |
| C-O | previous | 46.51% | previous | 46.76% | previous | 29.33% |

Primary MD-C, identical held three-group targets:
- zero: SSE 133431.82568, gain vs zero 0.000%, whole recorded-morphotype total SSE 88513.74526; opposition balanced accuracy 0.50000.
- previous: SSE 147435.03514, gain vs zero -10.495%, whole recorded-morphotype total SSE 156014.62236; opposition balanced accuracy 0.54554.
- cumulative: SSE 119865.86568, gain vs zero 10.167%, whole recorded-morphotype total SSE 105907.30650; opposition balanced accuracy 0.58036.
- bottom_history: SSE 239517.72179, gain vs zero -79.506%, whole recorded-morphotype total SSE 234226.16912; opposition balanced accuracy 0.55179.
- bottom_state: SSE 383220.45079, gain vs zero -187.203%, whole recorded-morphotype total SSE 419458.94656; opposition balanced accuracy 0.54107.
- direct_history: SSE 203189.23397, gain vs zero -52.279%, whole recorded-morphotype total SSE 216639.39844; opposition balanced accuracy 0.56339.
- direct_state: SSE 233004.53730, gain vs zero -74.624%, whole recorded-morphotype total SSE 257377.60758; opposition balanced accuracy 0.55536.
- direct_history_scaled: SSE 255383.63084, gain vs zero -91.396%, whole recorded-morphotype total SSE 253381.59471; opposition balanced accuracy 0.56875.
- direct_state_scaled: SSE 375890.20620, gain vs zero -181.710%, whole recorded-morphotype total SSE 420741.52605; opposition balanced accuracy 0.51964.
Cumulative mean is best group forecast for MD-C, M-C, D-C; last observed contrast wins C-O. Cumulative MD-C improves over zero at all three resolutions tested (morphotypes 1.67%, conservative groups 10.17%, published two-group totals 5.06%); these are relative gains within each target, not comparisons of absolute SSE across targets. Direct ridge is better than bottom-up ridge at penalty 10 in this study, but is worse than simple history for every contrast. Penalty sensitivity weakens the advantage. Do not infer that aggregation itself causes predictability, or that an improved historical mean beats a prospective forecasting baseline as a general mechanism. MD-C cumulative balanced accuracy 0.5804 vs constant classifier 0.5, despite positive opposition in 280/360 overlapping target cells. Whole recorded-morphotype total MD-C forecast best zero. 

## Interpretation boundary
Independent aquatic evidence for sustained opposite signed group response, not an independent confirmation of the Kansas intervention mechanism or causal assembly-history effect. Original authors already describe increased cyanobacterial dominance; our added diagnostic is out-of-pond/out-of-time forecast and fixed-group sign persistence with unresolved bounds. No new causal bridge, factorial interaction claim, formal meta-analysis, or ecological compensation of equal magnitude. No selective late-week scoring, no p-value using dependent cells, no tuning. Whole sum can include non-algal organisms. Freeze this diagnostic; prospective further source acquisition can target an actual randomized assembly-order experiment (Hayashi founder-size is not order manipulation). 
