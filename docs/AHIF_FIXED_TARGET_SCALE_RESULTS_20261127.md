# Fixed-training target aggregation — verified 2026-11-27

Frozen prescore 3192eeb889fde48e57ab6ceef8ab992591a9c447, parent PR26 9b0aa8e43e98e2686e594c9ab257d670e66c13ee. This is a retrospective diagnostic, not an independent publication or a geographic grain experiment. All prior results preserved.

## Existing temporal studies
Exact extraction of PR14 frozen run 45589332743, code 86a75c6a71f13cea51cd34e8a7f6cbf689130998, archived common-support results. 23 fixed-two-training source-size rows. Signed-TV error decreases from h1 to h2 in all seven sources; absolute mass-TV likewise. G direction gain increases in pond, POEM, field, Mancuso, Hayashi; decreases in Kansas and Saine. H increases in six, decreases Kansas. Source counts and common support unchanged. These normalized contrast forecasts change the target distribution with held replication; squared-error cancellation below does not apply to nonlinear TV normalization or absolute error. No new source rescore.

## Gibbs fixed predictions
Run 45951734930 SUCCESS, code a4af377e2ce65c45500235068771d5227e8c5970, job 140348915051, artifact 15176882763 SHA256 746fd6def131ffa51b65eff851eb2f3ad1c6291290873d9005059d4f28b6315b. Exact source audit, per-label and original variable-held group scores reproduced. Same 256 fixed 10-train/5-held label partitions. All 30 focal-pair groups retained; exact-five support 6325/7680 dependent fold-groups (37–256 eligible folds per group); at-least-three sensitivity 7643/7680. No missing-data imputation. Equal group weights then equal eligible folds. Enumerate all held subsets at each h, using identical predictions and support. 
| h | Zero MAE | Anonymous MAE | Held-pair identity MAE | Zero MSE | Identity MSE | Common target variance |
|---:|---:|---:|---:|---:|---:|---:|
|1|.509101|.534004|.532050|.593153|.531964|.413861|
|2|.405216|.393412|.387668|.334490|.273301|.155198|
|3|.356870|.331233|.322626|.248269|.187080|.068977|
|4|.328611|.293211|.283280|.205159|.143970|.025866|
|5|.309698|.267238|.256324|.179292|.118103|0|

Fixed prediction p, held pool k values, population variance V: mean over all size-h subsets of squared error = (p−full held mean)^2 + (k−h)/(h(k−1)) V. Exact cancellation of the target-variance term for absolute MSE contrasts between fixed predictors. All cells, models, measurement and support lanes pass (max error 2.23e−15). Held-pair additive identity absolute MSE gain vs zero .0611891194 for every h; vs anonymous .00550872703 for every h. Relative improvement vs zero increases as common denominator decreases. MAE crossover h1 to h2 is genuine under this loss and target but does not imply information growth. At least-three, hypothetical SA x10, AC/UR-only, density-standardized sensitivities retained; identity increment not uniformly positive vs anonymous (AC/UR-only standardized exact-five −.0019703). Note common target variance is finite observed subset variation, not exclusively measurement noise or stochastic observation error. Mathematical identity holds even with deterministic ecological variation within pool. Cannot infer spatial emergence, biological interaction strength, or relative importance of noise sources.

## Shen held inputs
Run 45951771253 SUCCESS, code ba54729558f5b17ca0e0aac64e10bade9b2c396d, job 140349031048, artifact 15176569226 SHA256 6403654124bb5575026eeb942e476b7406ad21aab97a17bd16f14489bd99fcea. Same source audit and earlier h2 scores reproduced. 512 folds ×8 dependent panels each lane, fixed two training cultures/arm. h1 average all four cross-products of held target/background culture, no physical pairing inferred; h2 means of same two cultures/arm. Linear predictors use held background, hence prediction variance is model-specific. 
| Window | MAE zero h1/h2 | MAE multiply h1/h2 | MSE gain multiply vs zero h1/h2 |
|---|---|---|---|
| common days3–17 | .164838/.150096 | .183525/.178554 | .004327/.003097 |
| strict days3–9 | .163715/.135498 | .233869/.220834 | −.038220/−.040247 |
| days9–17 | .228163/.224529 | .137448/.136475 | .071320/.069961 |
For linear predictors MSE h1−h2 equals variance of predictions plus variance of targets (after train baseline normalization). Maximum identity residual 4.27e−16. In days9–17 multiplication predicts better under both h1 and h2; absolute MSE advantage decreases slightly as aggregation reduces variability of no-change predictions more than multiplication. Clipped additive does not generally commute with averaging; no forced decomposition. Previous day3 baseline counterexample and initial imbalance limitations unchanged. No causal interpretation of predator effect.
 
## Decision
Treat common replicate averaging as an exact loss-geometry diagnostic, not evidence for a macroecological law. A real spatial scale test must alter spatially defined ecological context at fixed sample effort and use independent locations or held clusters. Compare same number of units contiguous vs dispersed, with spatial matching and blocked transfer; record composition/interaction opportunity differences, not just lower absolute error. Existing study-specific forecast targets are different; no pooled meta-analysis, no new publications, no signed interaction coefficients. 23 archived temporal summaries, 240 Gibbs and 18 Shen derived rows (281 total dependent records) archived separately from SQLite v0.8. No manuscript change.
