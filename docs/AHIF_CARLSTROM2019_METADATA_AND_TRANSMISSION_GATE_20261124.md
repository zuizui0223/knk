# AHIF v0.3 — Carlström 2019: metadata source audit and a resident transmission challenge

Date: 2026-11-24. Independent Draft PR #6. No changes to frozen KNK manuscript or main.

## Primary source and actual evidence

Carlström et al., 2019, Nature Ecology & Evolution, DOI 10.1038/s41559-019-0994-z.
Original full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC6774761/
Original code: https://github.com/cmfield/carlstrom2019
ENA accession PRJEB32997.

The authors experimentally removed and later introduced phylogenetic groups in a 62-strain Arabidopsis leaf microbiota and conducted separate single-strain dropout assays. The original paper already reports history-dependent assembly, a dropout-derived causal response network, and relative resistance of established residents. AHIF cannot claim novelty for any of those concepts. Published prose reports PERMANOVA R² ranges 0.19–0.62 for the introduced group, 0.12–0.28 for whole communities except Gammaproteobacteria, and 0.05–0.06 for founding taxa under early-vs-late contrasts. These are published summary ranges, NOT independent AHIF effects and NOT directly comparable to null-centered E. Late addition versus mock spray had no significant effect on the resident subset in the original report. Non-significance is not equivalence.

## Metadata and public script directly inspected via GitHub

| Original file | Git blob SHA-1 | Rows excluding header | Rows after author filter |
|---|---|---:|---:|
| metadata53.csv | a143701e52ac3324e5e1765d4252d2e00f52caf1 | 342 | 299 |
| metadata61.csv | 1056e7b5148b3d3825483178a29eb73f0ced4cc3 | 373 | 311 |
| metadata62.csv | 83cd360eb0c8428fad8c3586c74ec7d02c6bc840 | 373 | 311 |

All sample IDs unique within each file. Original analysis_publication.r blob 1bf958f12baacc790ac29c04120b266b61270c8f. Experiment 53 excludes Treatment=Ax and Spray=Inoc. Experiment 61 excludes Treatment=Axenic and Time=t0; 62 excludes Treatment=Ax and Time=t0. Total 1,088 metadata rows, 921 retained records; not a verified independent plant/box-level N. Experiment 53 t2 counts: ALL/Mg 24; No A/Mg 16 and No A/A 16; No B/Mg 15 and No B/B 16; No G/Mg 16 and No G/G 15; No P/Mg 16 and No P/P 15. Experiments 61 and 62 each have 26 conditions (control + 25 dropout targets); 12 records per condition except R1 -86 and R2 -88 have 11.

### Metadata alignment risk actually found

The original script assigns `meta62$initial <- droplevels(meta61$initial)` after preparing meta62. The correct source for that operation would be meta62 itself. Both filtered vectors have 311 entries, so the assignment can silently succeed. At 1-based filtered position 72, R1 sample 81 has treatment -88 and R2 sample 400 has treatment -86. Copying R1 labels into R2 would mislabel the latter record. All other filtered treatment positions match in the inspected CSVs. The affected source files and analysis script are blob-pinned above. This is a code/metadata discrepancy; its effect on published results cannot be quantified without the abundance matrices and exact analysis run. Independently join counts to sample IDs by explicit keys; never align by row position.

The repository's R script expects expt53, expt61 and expt62 otutab_final_classified.txt files. They are not in the displayed source repository file list. Original sequence data are archived under ENA PRJEB32997; sample-level outcome table has not been obtained in this audit. Do not substitute the publisher's control-only supplementary matrix for intervention outcomes or infer zero from absent raw files. No source-outcome audit or rank scoring is claimed.

## Narrow ecological hypothesis: resident transmission vs newcomer confinement

Initial manipulations can alter newcomer establishment without reorganizing incumbent dominance, or restructure the community beyond the manipulated group. Define the resident roster independently of observations, excluding the manipulated group from both arms. Compare early vs late arrival within a common final-time and fixed roster; compare late introduction vs time-matched mock spray among communities with identical starting exclusion. Close relative abundance within the resident roster only; separately record zero resident totals and, if available, absolute biomass/load. No claim of absolute effects from compositional counts. Do not include omitted/added taxa in the primary comparison, because structural zeros and species-number changes can mechanically produce apparent rank effects.

Pre-outcome analysis contract:
1. Obtain authenticated treatment count tables or reconstruct ENA reads with original reference sequences, documented mapping and QC. Resolve plant-level units, pooled plants, inoculation blocks, destructive versus repeated harvest and batch. Verify treatment counts against metadata and join by sample ID.
2. For each of the four group manipulations in experiment 53, score same-time control versus late arrival and late versus mock, separately for full common roster and fixed resident complement. Primary contrast is resident complement at t2, named vs within-sample sorted composition, with Bray-Curtis pseudo-R² and E=R² minus median design-respecting permutation null. Stop for zero totals or unknown exchangeability; record PERMDISP as a sensitivity. Preserve shared control dependence and nested Proteobacteria group overlap; do not pool as four independent papers.
3. Separate dropout experiments 61/62 could provide cross-batch, held-out predictions of net resident responses. Learn single-strain removal effects in 61, validate in 62 using corrected metadata and independent controls; repeat reversed as a secondary sensitivity. Predefine response breadth and signed effects using fixed resident roster. Compare against control dominance and richness/closure baselines. Class-level held-out transfer is exploratory because combined removals need not be additive, and batches/interventions differ.
4. Test whether cross-validated dropout response breadth predicts group-level resident architecture disruption. With only four partly nested group contrasts, no general inferential slope or meta-analysis. Claim only design-conditional predictive agreement or failure; independent manipulation of dynamic interaction coefficients remains absent.

Prediction: resident shielding gives large newcomer signal but weak resident rank signal; resident transmission gives nonzero resident rank signal, ideally validated against an independently assayed perturbation. Compare with the Swedish field experiment only after aligning intervention families and endpoint definitions; non-target plants are not proven later immigrants.

## DB gate

Seven substantive prior-art checks after this paper-level promotion; 30 study/candidate records, 32 panels, 118 aggregate outcomes, 10 outcome-source-audited studies unchanged. 921 metadata records are not new effect sizes. Same-panel independent causal mechanism-to-storage path remains unverified. No pooled effect. See reproducible audit script below.
