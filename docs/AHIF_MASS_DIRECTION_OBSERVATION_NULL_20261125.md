# Illustrative observation-null control (2026-11-25)

Executed independently of the six-source run, without viewing its new scores. Script benchmarks/simulate_ahif_mass_direction_observation_null.py, seed 2026112512, 6,000 simulated independent communities per scenario, 40 taxa, 6 observation times, independent Poisson counts for two training and two held replicates per arm, constant latent ecological abundances and treatment contrasts. Second scenario adds shared species-specific time-correlated observation multiplier (AR coefficient .85, log SD .55). These parameters are illustrative, not fit to any empirical system; no claim of generality or exact empirical null calibration.

| Scenario and old history | Recent mass alone gain pp | Recent mass + old signs gain vs old full pp | Hybrid gain vs recent full pp | Recent full gain vs old full pp |
|---|---:|---:|---:|---:|
| Constant observation, equal-one | .0326 | -2.2219 | -2.2802 | .0584 |
| Constant observation, older four-time mean | -5.4260 | -5.5783 | +1.6321 | -7.2104 |
| Time-correlated observation multiplier, equal-one | +1.5243 | -.6774 | -2.1824 | +1.5049 |
| Time-correlated observation multiplier, older four-time mean | -1.8648 | -2.1482 | +1.4756 | -3.6238 |

All simulated comparisons defined. The Monte Carlo SE of the equal-one recent-full gain is .0790/.0829 pp respectively; SE of the cumulative hybrid-vs-recent gain is .0234/.0242 pp. These quantify simulation precision, not empirical SE. Constant ecological differences can produce an advantage of history signs over recent signs and correlated observation noise can produce apparent recent magnitude predictability. In these two fixed scenarios, no hybrid beats BOTH complete-history and recent profiles. Mixing magnitude/sign estimates can break their estimation covariance, and averaging different numbers of dates changes noise. The simulation does not reproduce the full Portal pattern, and does not identify biological timescales from a forecast gain. No parameter tuning, no independent empirical evidence or meta-analysis.