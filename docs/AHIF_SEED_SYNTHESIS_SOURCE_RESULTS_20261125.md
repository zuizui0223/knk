# Ladouceur 2020 secondary synthesis: official source audit complete

2026-11-25 JST. Source article https://doi.org/10.1111/ele.13566; data https://doi.org/10.6084/m9.figshare.12319682.v1, CC BY 4.0. Author code repository emma-ladouceur/SeedAdditionSynthesis commit ba33e2b75d54fa37de8d09067c11978be7c5108b. Verified source run https://github.com/zuizui0223/knk/actions/runs/45539188760, job 138953217266, commit a32a06f9afd5651150e5922dddabbe28d18ea7f0; first attempt had manifest download_url naming mismatch, second had strict UTF8 decoding failure, third succeeded. Two original archived source files match repository-declared MD5/size; both reversibly decoded cp1252 (some preexisting mojibake in reference strings, do not automatically rewrite taxa). 
- SeedAdd_Sp_level.csv: file 22714958, 23301361 bytes, MD5 420281357484766c25df281a4cd2ed45, SHA256 d9f09ea90b3bac6e9913000d93acf24e1c8b0293bebac563477c9e90a9fd590e
- SeedAdd_Plot_Level.csv: file 22714961, 962570 bytes, MD5 715c3ed27e6c11d41185a02a2ab07b90, SHA256 0455c64afd914af5c7f85273ae950104ece31bf94fe9737f1a172aa6d9f437a3
3393 plot-state rows, 58426 species-state records, exact unique.id_ join, no duplicate state/species keys or plot/species-only IDs, required metadata fields agree. 12 experiment identifiers, 987 controls and 2406 seeded states. Source table records plot, site, block, yr.trt; repeated-state IDs not independent plots. Some species-table seeded labels use Seed vs Seeds; seed.rich matches across files. Species-level biomass.measure direct 3212 records, indirect 55214 records. Actual state/experiment measurement-mode counts to be computed next. Source describes indirect biomass allocated from cover and community mass from neighboring plots/strips, so closure/allocation can mechanically create species cancellation even without directly measured compensatory masses. No individual seed membership column in either table. Species addition cannot be equated with successful immigration, and rising species cannot be called colonists. 
| Experiment | States | Matched zero/positive seed site/block/year strata | Species rows | Mass sum within 1% of reported total |
|---|---:|---:|---:|---:|
| ASGA_Michigan | 48 | 24 | 511 | 48 |
| CCR_04 | 28 | 2 | 197 | 14 |
| CCR_093 | 960 | 222 | 15695 | 925 |
| California_Invade | 80 | 40 | 720 | 80 |
| California_Prop_Limi | 531 | 36 | 5385 | 530 |
| Germany_Montane | 232 | 58 | 5306 | 232 |
| Halle | 120 | 20 | 1782 | 120 |
| Jena | 1020 | 6 | 22408 | 1020 |
| Jena2 | 36 | 6 | 688 | 36 |
| Kansas_KUFS_LTER_Hay_Meadow_Exp_2 | 120 | 60 | 1572 | 120 |
| Kansas_KUFS_LTER_Old_Field_Exp_1 | 168 | 84 | 3849 | 168 |
| Texas_Temple_Prarie | 50 | 1 | 313 | 38 |
Source original Table 1 distinguishes study durations, sites, doses; some control/treatment imbalance and raw species-mass totals disagree (especially CCR_04, Texas, selected CCR_093 and California_Prop_Limi). No rescaling/imputing to close totals. These are curated secondary data, not 12 independently audited original experiments. Original article already models diversity and biomass responses and compositional dissimilarity; new diagnostic addresses opposed signed species mass responses. Jena and CCR original study family overlap must be resolved before meta-analysis. All original data ephemeral; aggregate audit validation/ahif_v1/ladouceur2020_synthesis_source_audit_v1.json. No causal bridge or DB effect-count increment.