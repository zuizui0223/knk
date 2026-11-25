# Ishizawa 2024: equal-information controls overturn the explicit-extrapolation interpretation

2026-11-25 JST. Source 10.5281/zenodo.8433832; original article https://doi.org/10.1073/pnas.2312396121; scoring run https://github.com/zuizui0223/knk/actions/runs/45526767902, job 138911956589, commit b9a81fbcf31bacce8b2db518bec41f09c2f7a708. Seven original MD5s, four unit tests, real-data scorer and artifact upload succeeded. Original author point predictions reproduced within 8e-15 log CFU. Same seven-strain experiment as prior analysis; not independent replication. All 64 fixed held-target subsets and nine prespecified models retained. No target outcome used in fitting or selecting models. Prescore docs/AHIF_ISHIZAWA2024_CONTROLS_PRESCORE_20261125.md.

## Equal-information model comparison

| Model | Named TV | Rank TV | Assignment residual TV | Simpson abs error | Top1 abs error | Log summed geometric CFU abs error |
|---|---:|---:|---:|---:|---:|---:|
| mono | 0.291917 | 0.165258 | 0.126659 | 0.097423 | 0.117439 | 0.206935 |
| pair_mean | 0.160250 | 0.135834 | 0.024416 | 0.079993 | 0.120612 | 0.132510 |
| pair | 0.213103 | 0.179880 | 0.033223 | 0.152278 | 0.157436 | 0.206937 |
| trio_mean | 0.093644 | 0.084737 | 0.008907 | 0.051968 | 0.064257 | 0.074386 |
| pair_on_trio_rhs | 0.172297 | 0.151052 | 0.021245 | 0.126882 | 0.140461 | 0.123795 |
| global_trio_coeff | 0.114524 | 0.106457 | 0.008067 | 0.073984 | 0.093059 | 0.082067 |
| trio | 0.116690 | 0.108842 | 0.007848 | 0.077680 | 0.093634 | 0.087324 |
| additive_log | 0.129903 | 0.118286 | 0.011617 | 0.090978 | 0.107379 | 0.074781 |
| inclusion_exclusion | 0.264057 | 0.144991 | 0.119066 | 0.095987 | 0.117066 | 0.152097 |

| Baseline → enhanced | Named relative gain | Rank relative gain | Improved rank targets |
|---|---:|---:|---:|
| mono → pair | 26.999% | -8.848% | 26/64 |
| pair → trio | 45.242% | 39.491% | 52/64 |
| pair_mean → trio_mean | 41.564% | 37.617% | 54/64 |
| trio_mean → trio | -24.610% | -28.447% | 25/64 |
| pair_on_trio_rhs → trio | 32.274% | 27.944% | 46/64 |
| global_trio_coeff → trio | -1.891% | -2.240% | 23/64 |
| additive_log → trio | 10.171% | 7.984% | 34/64 |
| inclusion_exclusion → trio | 55.809% | 24.932% | 36/64 |

**Important falsification:** the direct trio mean has the lowest mean error of all nine prespecified models on named/rank/Simpson/top1 and log summed geometric CFU (assignment residual alone is slightly lower in original trio). It is better than the original trio model on named TV by 19.750% and rank TV by 22.147%; do not use the reverse-denominator degradation as reduction. Direct mean pair→trio already improves named TV 41.56% and rank TV 37.62% (54/64), comparable to original pair→trio (45.24%,39.49%). The simple trio mean beats original trio in both named and rank mean error in all four target-richness groups. Replacing pair coefficients with contextual trio coefficients at fixed trio RHS/extrapolation factor does improve rank 27.94% and named 32.27%, but the zero-extrapolation control is better than either. Target-specific coefficients do not improve over all-third-strain mean coefficients (rank -2.24%, named -1.89%). Full seven-member target global and local agree, as required. The previous interpretation that explicit higher-order extrapolation was needed for rank prediction is not supported by these stronger controls. Three-strain observations contain contextual information, but data quantity, saturation and model misspecification are not separately identified; direct mean is not an interaction-free assay. No claim of a biological absence of higher-order interactions or proof of saturation. The finite-panel comparisons are descriptive; no p-values, no independent 64-experiment count.

## Matched error decomposition

Mean pair→trio direct-mean improvement named 0.066606 = rank 0.051097 + assignment 0.015509. This mathematical identity is not a causal mechanism partition. Original pair vs mono worsens rank, but direct pair mean vs mono improves rank; pairwise observations are therefore not inherently detrimental to rank, rather one pair extrapolation model mispredicts concentration. 

## Source grain audit

692 system × rep label groups; zero duplicated system-rep-strain rows; zero incomplete/duplicated roster groups. This supports complete label-level linkage, but physical flask pairing and independence have not been separately verified. No replicate-level CI or raw data redistribution. Cell geometric means, not mean within-flask composition.

## Stop rule

Archive all nine models and eight contrasts, all size-stratified results and 64 target rows at validation/ahif_v1/ishizawa2024_equal_information_v1.json. Do not optimize coefficients against these same held targets. Register finite-panel prediction lane separately from intervention effects and held-location CV. Next empirical priority is a new independent time-ordered study where partner effects are independently measured; not another reparameterization of this dataset. No causal history bridge; no meta-analysis eligibility. Main/manuscript unchanged.