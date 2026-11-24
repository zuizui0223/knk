# Source methodology guard — Zou & Rudolf (2026) paired interaction and rank

**Updated 2026-11-24; post hoc research question selected after previewing the source author repo's group means. Not preregistered confirmatory.** Paper `Ecology` DOI https://doi.org/10.1002/ecy.70454, public data/code DOI https://doi.org/10.6084/m9.figshare.32118556.v2. Authors' original GitHub `hengxingzou/Zou2026Ecology` pinned at commit `5f28e5da78853bb31de1e1fa86f20f4892d7f2bb`. Original source Git blob hashes:
- `Counts.csv`: `8fef4ec5818319452e1febb4c51042530b39a676`
- `Coefficients/CompCoefficients.csv`: `5071d3d4bfdfb831326e5eb49d3960c0902c29fc`
- `Coefficients/Lambdas.csv`: `3380ce79795a8ca2419524480cf6a70e0acdfbd3`

**Original biology already established:** arrival timing modified two `Tribolium` competitors' fitted Lotka–Volterra competition parameters and adult outcomes, and generation number modulated stabilization/relative fitness differences. No claim that AHIF discovered these facts.

## Raw grain / program stops
- `Counts.csv`: 414 unique whole **vials** (3 season-generation conditions ×138 per generation). Across each season number, six mixed starting density combinations × 3 vials × seven arrival times =126 mixed; single-species controls (12) at simultaneous -> **378 mixed, 36 single controls**.
- Each mixture type includes initial `T:F` pairs (10:20, 15:15, 20:10, 20:40, 30:30, 40:20). 3 independent **vials** per time × density × generation. Other columns are dates and counts of adults **at end** of generation-specific season; not independent repeated surveys.
- `Relative_Arriv_Time`: −18, −12, −6 means *T. castaneum* (T) earlier; +6,+12,+18 means *T. confusum* (F) earlier; 0 simultaneous.
- Author-produced `CompCoefficients.csv`: 2 fitted **nonnegative interspecific coefficients** alpha12 and alpha21 × 7 relative days × 3 generation settings = 42 rows. They are **model posterior summary estimates FIT TO THE SAME Counts.csv DATA**. They are not individual vial-level experimental measurements and their posterior medians/10–90% intervals do not supply independent causal predictive validation.
- If any source hash mismatches, number of vials/coefficients/missing values wrong, independent vial names duplicate, or grouping incomplete, STOP. Do not silently discard vials or fill zero adults with pseudocounts. Each mixed vial must have positive total adults for proportion analyses.

## Focused, exploratory representation audit
The intended ecological interpretation differs from the article's original priority-coefficient claim. Restrict to **balanced initial density** (T=F=15 or T=F=30), yielding **126 whole vials**, 3 original replicates ×2 density strata×7 arrival times×3 generation treatments. This holds *relative starting species frequency constant* while varying total density.

For each generation × lead-duration `d∈{6,12,18}`, compare **T-first (−d)** with **F-first (+d)**, within total density (15:15 and 30:30):
- `T_share=Adults_T/(Adults_T+Adults_F)`: **identity-specific** dominance of T;
- `Top1_share=max(Adults_T,Adults_F)/(Adults_T+Adults_F)`: the **identity-free two-species rank concentration**. These are from the **same raw 2-element species vector**.
- Compute observed `ΔT_share=mean(T_share|-d)−mean(T_share|+d)` and `ΔTop1=mean(Top1|-d)−mean(Top1|+d)`. Each mean weights six real vials equally with equal density representation.
- For the 3-vs-3 assignment within each density level, enumerate all \(\binom{6}{3}^2=400\) exact combinations; report **two-sided permutation reference p** for both descriptors separately. Interpretation as a source-experiment randomized p requires verification of vial assignment and timing-block exchangeability. It is a **conditional label-permutation diagnostic**.
- Report **all nine 3-gen×3-lead contrasts**, not just the strongest and no "winner switch" selection. Compare descriptive magnitudes within the nine paired descriptors; note that T_share has a fixed orientation and Top1_share is a nonlinear fold under label erasure. They are different descriptive scales, not an information-theoretic entropy loss.
- Independently report 3×7 `T_share,Top1_share` means for all balanced vials to show sign/asymmetry and nonmonotonicity. `p>.05` never proves history-free architecture or equivalence.

## Coefficient and endpoint comparison — risk of information leakage
For each generation×lead day, align **author-published fitted coefficient medians** with the mean adult T_share across six balanced vials. Report the signed interaction asymmetry `log(alpha21/alpha12)` (alpha21 T effect on F and alpha12 F effect on T, from original two-species Stan matrix). Check coefficient positivity and exact 3×7 coverage. Report **Pearson correlation descriptively** across 21 condition-level cells and within-generation 7 cells, **NO p value** (coefficients fitted from same adult-count outcome, not independent observations, coefficient uncertainties correlated). Do NOT claim causal mediation or held-out predictive accuracy.

**Interpretation gate:** This is a **2-species** study, so "rank" means Top1 among 2 and **cannot** validate the proposed >=11 species generalization. The independent source demonstrates that even when history effects on the species label are large, unlabeled dominance may be relatively similar for opposite first-arriver histories, but these comparisons may still change absolute concentration. The prior study authors already demonstrated fitted interaction coefficients changing with arrival; this source is a powerful novelty challenger, not proof of a universally new discovery.

## AHIF contribution and limits
Before source scoring, the AHIF pilot records `ST_GENERATIONS2026` and `P_GENERATIONS` were candidate-only. If source-faithful results pass, mark as original-source-audited and numeric endpoint status `source_reanalysis`, **never as independent >=3 multispecies structure mediation**. This is one publication and one experiment with multiple repeated treatment cells, not 21 independent studies. Do not add it to the same meta effect pool as plant group, fungal wood, or soil chemistry.

One source study can be considered **joint interaction and community outcome at the conditional treatment level** but the alpha estimates are fitted from the same count outcomes. Explicitly code `same_panel_causal_interaction_to_structure_identified=no` to prevent circular validation.
