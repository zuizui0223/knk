# Post hoc, source-corrected Delory process–output contrast, 2026-11-24

**This is explicitly post hoc**, motivated by the original-author-filtered 58-plant results and **not a preregistered confirmatory new study**. The original publication already discusses root foraging and N uptake dissociation. Reuse immutable GitHub source `BenjaminDelory/Rhizometabolome_2020@v1.0.0` `data/Data_plants.csv`, Git blob `d11b4b45acb13bd9c240b004332e76b09eb6c3e3`; follow author R notebook's `data<-data[-c(3,22),]` filtering. Preserve full 60-plant unfiltered sensitivity.

### Exact question and estimand

Instead of treating one significant RSD result and one non-significant total-biomass/Nuptake result as evidence that responses differ, **directly test a scale-free difference in treatment log ratios** for the SAME individual plants:

```
D(outcome) = log[mean(RSD|legacy) / mean(RSD|Soil)]
           - log[mean(outcome|legacy) / mean(outcome|Soil)]
outcome = PDW or Nuptake.
```

A negative D says root-surface-density has a **more negative proportional change** than the output, not that total biomass or nitrogen uptake is identical to the control.

For all **2 species ×2 legacy sources ×2 output endpoints =8 contrasts** (no subset selection based on results), resample within source Species×Treatment cell **the SAME plant row indices across RSD/PDW/Nuptake**, 20,000 resamples per cell (seed 20261124). Report observed D, 95% paired-plant bootstrap percentile CI; zero-inclusion cautiously, plus source-filtering sensitivity (58 vs 60). This is not an independent new experimental replicate.

**Main caveat:** the source plants share soil-solution source/tray conditions and the original archive has no tray assignment IDs. Plant-level bootstrap ignores possible cluster variance; CIs are conditional and cannot establish a causal mediation chain. RSD is derived from RDW×SRA, and N uptake from tissue nutrient mass, so the metrics share component measures. Process–output contrast is a measured phenotype relation, not a mechanistic root-resource flux measurement.

### Interpretation

Only if the direct D contrast differs from zero under the conditional bootstrap can we state a **descriptive larger proportional RSD reduction than PDW or Nuptake in the analyzed plants**; do not claim "biomass unaffected". Treat any effect sensitive to the two author-selected exclusions as fragility rather than universal biological regularity. Keep current KNK v28–v38, existing Ecology Letters manuscript, and separate experimental PRs unchanged. Documented as a new secondary exploratory extension of Draft PR #5.
