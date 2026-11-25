# Training-only contrast reliability: cross-source prediction results

Verified source-locked run https://github.com/zuizui0223/knk/actions/runs/45584955172, job 139104085121; artifact 15008132164 SHA256 d88e3c70a98de71a9aa4b5a4e0e4efbd982082f9524fa657e5577786afa39059. Exact summary validation/ahif_v1/signal_transfer_summary_20261125.json. All source checks and archived direction-score replays passed. Four absolute-measurement sources, 19 panels and 3470 dependent fold-date cells; held microbial compositional source, 8 panels and 1152 dependent cells. No new publications. All cells eligible. Kansas source linkage conditional; POEM and field potentially related teams.

Primary endpoint signed vs signblind forecast gain, secondary same-direction gross share. Training-only features at prior date, fixed ridge and equal source/panel/time weighting. U=(Q-W)/(Q+W); Q observed squared contrast, W estimated contrast variance trace. U is noisy, not identified true ecological signal. Positive improvement = relative held MSE reduction.

| Validation | U addition over strength: sign gain / same share | U model vs intercept: sign gain / same share |
|---|---:|---:|
| Four-source leave-one-out | +7.71% / +9.65% | -227.11% / -150.19% |
| POEM+field jointly held | +2.90% / +1.95% | -213.27% / -123.52% |
| Kansas excluded | +4.19% / +18.38% | +3.41% / +10.96% |
| Held microbial source, full panel | +20.46% / +21.38% | +18.09% / +12.23% |
| Microbial 3 depth-eligible panels (unrarefied) | +36.74% / +45.04% | +42.31% / +29.02% |
| Held microbial source, Kansas excluded from training | +11.49% / +17.59% | +3.35% / +34.49% |

Crucial qualification: four-source signed-gain MSE intercept .026202, structure .082392, strength .092869, U .085711, noise .087425. Same-share MSE intercept .076531, structure .194200, strength .211924, U .191472, noise .201904. Intercept far better than complex models. U addition reduces error from weak strength model but does not establish useful general cross-source predictor. Excluding Kansas structure remains better than U for both outcomes. Held microbial transfer is positive on this one previously examined publication and measurement-lane transfer is descriptive, not independent prospective validation. Alternative scaled variance (noise) also improves over strength in microbial transfer (+20.00%, +22.92%), and is stronger than U on three depth-eligible panels. Do not attribute gain uniquely to U. Correlations between U and signed gain within source (dependent weighted cells): pond .2594, Kansas -.0495, POEM -.2545, field .4910, microbial .5342; same share .2250, .0316, .2259, .4294, .4557. Heterogeneity explicit, no inferential p-values.

Conclusion: cross-source generalization not established; narrower observation that relative variance information can improve prediction conditional on strength model, and transfers to one microbial panel set. Static persistent treatment effects and detection noise can generate directional predictability; not evidence of dynamical memory or a biological interaction mechanism. Do not pool fold scores as meta-analysis, change weighting or tune models after results. Next priority: effect sizes at study-level with physical-unit uncertainty and independent publications, plus simple baselines, rather than further model optimization on these four sources. Prior manuscripts, effect ledger, and source branches unchanged.
