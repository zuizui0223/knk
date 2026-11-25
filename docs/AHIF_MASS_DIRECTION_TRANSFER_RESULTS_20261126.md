# Response allocation versus direction history: cross-source transfer results — 2026-11-26 JST

Frozen prescore: docs/AHIF_MASS_DIRECTION_TRANSFER_PRESCORE_20261125.md (1f3e011ecaa6f4a0f62460b4eb1b9703bdba78a9). Initial monolithic run 45599288300 stopped at official van Steijn Zenodo acquisition HTTP 504; no ecological output. Source-isolated retry run https://github.com/zuizui0223/knk/actions/runs/45659798054, code commit 3b6c36bdabc7df6c88c680ab3f85350706025d8c: six jobs SUCCESS, original source hash/roster and prior-score replays passed. Source-specific job and artifact hashes archived in validation/ahif_v1/mass_direction_transfer_summary_20261126.json. Scoring function, 2+2 physical splits, target dates, shared control and source weighting unchanged. No redefinition of outcome or selective source omission. Saine remains STOP: only two post-introduction dates. Portal is the motivating separate pilot, not included in this transfer.

## Primary equal one-observation result

33 contrast panels, 67,072 dependent split-date-orientation cells (all defined), six existing source families. Gains in percentage points of TV prediction error. Mass/sign contributions are symmetric two-order averages; exact sum equals old-full minus recent-full loss. They are not separately identified causal effects.

| Source | Mass symmetric | Sign symmetric | Total | Recent mass + old sign vs recent full |
|---|---:|---:|---:|---:|
| pond | +3.4718 | +1.1013 | +4.5731 | -1.1047 |
| Kansas | +2.0178 | -.3388 | +1.6790 | -.4912 |
| POEM | +13.3539 | +1.7225 | +15.0764 | -2.7147 |
| field | +16.9085 | +2.6125 | +19.5210 | -4.4477 |
| Mancuso | +10.5167 | +2.0262 | +12.5429 | -2.8047 |
| Hayashi | +4.7950 | +2.8140 | +7.6090 | -3.0692 |
| Equal-source descriptive mean | +8.5106 | +1.6563 | +10.1669 | |

All six mass contributions positive and larger than their sign contribution; all six total gains positive. Pure absolute-mass TV gain positive in all six as well. The ratio of equal-source mean mass contribution to equal-source mean total gain is about 83.7%, an algebraic descriptive ratio, not a pooled causal fraction or a percentage of independently explained variation. For every source, the full recent profile beats the proposed recent-mass/old-sign hybrid at the primary equal-one comparison. Thus Portal's particular hybrid superiority is not a general equal-information result. Kansas's negative symmetric sign contribution does not imply the hybrid wins: the sign update on recent mass is positive even though the update on old mass is negative.

## Cumulative and equal-two sensitivities

Cumulative older contrasts vs recent one date: hybrid beats both original full models in Kansas only (versus old +.4685 pp and versus recent +1.2219 pp). Restricting to targets where cumulative includes at least two dates leaves +.5885 and +.8590 pp respectively. Pure absolute-mass TV comparison in Kansas instead favors cumulative by 1.2491 pp on all primary dates, 1.0567 pp on longer-history dates. Therefore this hybrid result depends on conditioning signs on the mass distribution, not evidence of purely more local response allocation in Kansas.

Equal two nonoverlapping dates vs two more recent dates (pond, Kansas, Mancuso, only eligible three sources): mass symmetric +11.8495/+4.4291/+9.2818 pp, sign symmetric +1.7572/-.2749/+.0584 pp; hybrid over recent -1.8316/+.3025/-.1946 pp. Positive mass effect persists on these three source means; older-sign hybrid superiority only small and Kansas-specific. POEM has just three dates; field/Hayashi four, ineligible for the two-vs-two forecast. Cumulative POEM result duplicates its single-old result and is not independent.

## Biological limits

The robust descriptive pattern here is a larger time-local contribution of normalized taxon response allocation than of sign updating, on these fixed endpoints and source designs. It does NOT identify distinct biological memory timescales or rule out observation-process or estimation effects. Equal date count is not equal elapsed time or information quality; temporal covariance, changes in abundance/detection, unequal error, and loss geometry remain. Four-way mixing changes covariance of estimated signs and mass; exact decomposition is a predictive loss identity. The fixed illustrative observation-null simulation in docs/AHIF_MASS_DIRECTION_OBSERVATION_NULL_20261125.md was rerun with original seed and its rounded results reproduced, not calibrated to these experiments. Kansas source linkage conditional; Hayashi density-position aliased; compositional and absolute measures differ; repeated cells dependent. No calibrated study SE, p values, formal meta-analysis, causal claim, new publication or submission. All 25 defined source-comparison rows archived separately in database/ahif_mechanism_v0/mass_direction_transfer_effects_20261126.csv; STOP rows retained in JSON. Original AHIF v1 and previous registry/source branches frozen. Next priority is independent source recruitment and a sampling/observation-model calibrated validation, not repeated tuning of these same families.
