# ReplayEcology: initial profile persistence and representation, verified follow-up

2026-11-25 JST. Original hash-verified processed count and metadata files, 275 parent communities, 14 origins, 1100 final endpoints averaged within parent, 1454 ASVs. Same experiment and source matrices as prior state result. Run https://github.com/zuizui0223/knk/actions/runs/45525145757, job 138906509794, scorer commit af6090908d3e3b94cf8eb3841498b74d210a02bb; all four synthetic tests, exact source file checks, real-data scoring and artifact upload passed. Prespecified exploratory models in docs/AHIF_REPLAY2025_PERSISTENCE_PRESCORE_20261125.md. No p-values for this follow-up; descriptive held-origin errors, no independent new experiment.

## Residualized profile persistence

| Endpoint | Train mean SSE | Unchanged initial SSE | Depth SSE | Summary SSE | Depth+persistence SSE | Summary+persistence SSE | Summary gain over depth+persistence | Improved locations |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| named | 32.89770 | 52.71615 | 32.80070 | 32.48385 | 23.62522 | 23.47163 | 0.650% | 7/14 |
| rank | 6.96288 | 18.07805 | 7.07103 | 6.25437 | 6.56257 | 6.20905 | 5.387% | 9/14 |

Unchanged initial profile predicts worse than training mean for both representations. Train-only residualized persistence improves over depth alone for named (27.97%) and rank (7.19%). Adding baseline richness/Simpson/top1 over depth+persistence improves named by 0.65% (7/14 origins) and rank by 5.39% (9/14). This is an incremental score against one fitted linear persistence model, not proof of nonlinear dynamics or an inferential test. Original 11.55% rank gain was versus depth alone; the stricter estimate is smaller.

## Fixed-penalty representation diagnostics

| Target | Lambda | Initial rank over summary | Initial named over summary | Initial named over summary+rank | Initial rank over summary+named |
|---|---:|---:|---:|---:|---:|
| named | 0.1 | -6.110% | 2.978% | 4.693% | -4.235% |
| named | 1 | -1.726% | 17.451% | 17.249% | -1.975% |
| named | 10 | 0.033% | 13.120% | 12.513% | -0.665% |
| rank | 0.1 | 3.847% | -17.933% | -27.433% | -3.899% |
| rank | 1 | 6.761% | 8.080% | 0.221% | -1.210% |
| rank | 10 | 5.690% | 7.568% | 5.067% | 3.138% |

For named target, initial named kernel improves the summary+rank baseline at every penalty (4.69–17.25%); added rank to summary+named is negative at every penalty. For rank target, initial rank improves summary at every penalty (3.85–6.76%); adding initial named to summary+rank ranges −27.43% to +5.07%, including +0.22% at primary lambda=1. Named includes rank information nonlinearly; these are fixed model-specific comparisons, NOT pure information partition or significant differences. No heldout tuning. Do not claim taxonomic identity irrelevant for rank. The original article (Nature Communications 2025 DOI 10.1038/s41467-025-57591-2) already linked initial labeled composition to trajectories; initial-state dependence not novel.

## Next gate
Stop expanding same Replay data with selected hypotheses until independent intervention source is available. Current synthesis: identity-free initial summaries and rank profiles can improve held-origin prediction in this spontaneous-succession experiment; Saine predictive R-loss leads concern an intervention and are a different estimand. No causal mechanism, general effect-size comparison, statistical test of representational difference or new independent study here. All fold losses archived at validation/ahif_v1/replay2025_persistence_v1.json, separate from DB's 161 effect records and 6 prediction records. No main/manuscript changes.