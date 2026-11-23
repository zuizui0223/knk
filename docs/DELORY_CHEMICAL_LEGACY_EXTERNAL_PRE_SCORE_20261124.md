# Delory et al. chemical legacy: independent post-submission process-channel validation

**Protocol dated 2026-11-24 and committed before KNK-specific outcomes were calculated.** Delory et al. (2021), *Plant, Cell & Environment*, DOI https://doi.org/10.1111/pce.13999, source DOI https://doi.org/10.5281/zenodo.4010162, immutable GitHub tree `BenjaminDelory/Rhizometabolome_2020` tag `v1.0.0`.

## Why this source

This is a two-phase environmental conditioning experiment: forb/grass/unsown soil communities create distinct soil solutions, then two **later-arriving plant species grow alone** but receive these solution legacies. It isolates the environmental-conditioning *route* from direct contemporaneous interspecific competition; it does not quantify the causal contributions of individual metabolites, microbes, nutrients, or root-exudation pathways. Original authors ALREADY showed reduced root foraging in Dianthus but not Festuca, with little direct relation to N uptake. KNK's reanalysis is an independent source-backed *process-versus-function* comparison, not a novel discovery of that mechanism.

## Fixed source and sample audit

Original source `data/Data_plants.csv` Git blob SHA `d11b4b45acb13bd9c240b004332e76b09eb6c3e3`, 60 plant observations, 44 columns, species `Dd` (Dianthus deltoides), `Fr` (Festuca rubra), conditioning `Soil` (unsown control), `Forbs`, `Grasses`. Expect exactly 10 replicates for every one of the 2×3 cells, unique `SampleID` and `Species×Treatment×Replicate`, all nonnegative finite values for the selected endpoints. Missingness is not zero. R source analyses and README `data/README.txt` are available in the original tag and must be matched.

## Primary process-vs-output targets (frozen)

Per species, compare both legacy treatments to `Soil` for:

1. **RSD**: root surface density (cm²/cm³), original proxy for root foraging (primary mechanistic process state).
2. **PDW**: total plant dry biomass (g), aggregate plant output.
3. **Nuptake**: total N uptake (g), nutrient acquisition (aggregate output).

For each species×legacy×response report mean in treated and soil conditions, absolute mean difference, ratio of means, log ratio of means. Independently resample *whole individual plants* within each Species×Treatment cell **20,000 times**, source-blind seed 20261124, and report bootstrap percentile 95% intervals for absolute and log ratio. Include the cross-species contrast in log ratios for each legacy (Dianthus minus Festuca), using independent plant cell bootstrap. **There are six responses per species across two legacy treatments (12 planned species×legacy×outcome contrasts)**; adjust any exploratory unblocked treatment-label permutation p-values by Holm across the 12. Do not call observational differences proof of equivalence or no impact.

## Secondary sensitivity

`RDW, SDW, SRL, SRA, RTD, RLD, RMF, TRL` are additional ecological states, but **secondary** only: report descriptive contrasts and compatibility intervals, do not select significant outcomes to sell a mechanism. `Nuptake` is the source's derived uptake index, not direct isotope-tracer flux. Root exploration and total uptake are different variables; lack of significant response in uptake does NOT prove constancy.

## Statistical and identification guardrails

Soil solutions were created in phase 1 and applied in phase 2; delivery was grouped by trays (the paper describes pots within a tray receiving the same solution). Whether independent plant pots are a valid treatment-randomization unit requires explicit source-method confirmation; source has **no tray ID**, so a plant-wise bootstrap is **conditional/descriptive** and may underestimate uncertainty if tray or solution-batch effects are shared. No claims of exact randomized p-values or resolved mediation without tray allocation records.

This is **one original study** with correlated multiple endpoints and species, **not 12 independent meta effects**. It tests chemical/soil-solution legacy *in absence of a resident plant at the second-stage assay*, not live competition; therefore do not pool its effects with the van Steijn live-neighbor timing contrast as though they share the same causal intervention.

## Decision

A large RSD response with weaker biomass/N responses would support **process–output decoupling** in this independent manipulation but it is a replication of the original study. Test a universal storage rule only after multiple published systems contain independently measured interaction-process and community-output responses. Leave KNK v28–v38 and current Ecology Letters manuscript untouched. Store only derived aggregate statistics, never original raw source, on this separate PR branch.
