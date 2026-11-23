# Delory original-source analysis correction: 60 measured plants vs 58 analyzed

**Date: 2026-11-24.** Correction discovered by directly reopening the original `R/Analysis_PlantData_RhizoMetabolome_2020.Rmd` at immutable `v1.0.0` from https://github.com/BenjaminDelory/Rhizometabolome_2020. This audit was undertaken **after** our initial prescored 60-plant result; preserve that first result as an unfiltered sensitivity, not as a reproduction of the published paper.

## Exact original author rule

In the original R notebook's **Data filtering** section (source file lines 156–166), the authors explicitly stated:

- Plant **3** received **50 mL** rather than **25 mL** soil solution on day 2 (treatment manipulation error).
- Plant **22** had a very different phenotype from other plants (manuscript-identified outlier).

The original R code is `data<-data[-c(3,22),]`; the `SampleID` column in released `data/Data_plants.csv` increases sequentially from 1 to 60. Thus row 3 and row 22 coincide with SampleID 3 and SampleID 22. Since both are *Dianthus*, exclusion leaves a **9,10,9** balance across the Forbs, Grasses and Soil solutions, and *Festuca* remains 10/condition. Prior 60-plant reanalysis used an explicit 10/cell gate and therefore cannot be called source-method replication.

## Corrected analysis and retained comparison

1. **PRIMARY SOURCE-FAITHFUL**: recompute all 12 prespecified Species×Legacy×(RSD,PDW,Nuptake) contrasts on 58 plants (exclude precisely SampleID 3, 22). Keep each observed individual as source-defined unit but note solution tray/batch IDs are missing.
2. **SENSITIVITY ALL MEASURED**: preserve the initial 60-plant analysis unchanged, with exact archive SHA and counts, and compare effect direction and interval.
3. Both recomputations use the original precommitted bootstrap 20,000, label permutations 9,999 (conditional diagnostic ONLY), 12-effect Holm family and cross-species log-ratio contrasts; maintain seeds and response definitions. No new endpoint selected based on result.
4. Full report must explicitly identify whether exclusion changes root-foraging inference. The original paper uses bootstrap compatibility intervals, SEM, and other traits; this narrower mean-difference analysis does NOT automatically invalidate/replicate the original structural equation model.
5. Plant 22 is an **original-author designated outlier**, not an arbitrary decision made by KNK after observing the desired effect. Still, report both analyses because source-author selection can matter.
6. Plant pots were clustered in soil-solution trays. Because tray IDs are absent, bootstrap/permutation at plant grain is conditional and may be anti-conservative relative to true randomized solution batches. Do not claim population-level or universal mechanism, and do not pool with KNK direct neighbor experiment as the same causal intervention.

## Source content and hashes

Source plant CSV Git blob `d11b4b45acb13bd9c240b004332e76b09eb6c3e3`, published R script at `v1.0.0`. User-facing data rights: original source remains outside KNK; derived aggregates only in the independent unmerged Draft PR #5. Existing frozen v28–v38 and manuscript are unchanged.
