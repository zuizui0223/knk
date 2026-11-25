# Lear 2024: resident count response cancellation — source-verified boundary condition

Source article https://doi.org/10.1111/ele.14493; archive https://doi.org/10.5281/zenodo.12528023; successful run https://github.com/zuizui0223/knk/actions/runs/45527948842, job 138915847327, commit 2267688c5d6fd7468ad501bdac969278ba1c4b4b. Three synthetic tests, exact three original source MD5 checks, endpoint schema/coverage audit, scoring and artifact upload passed. This is an independent source experiment, not a replication of previous model-based effect estimands. Original paper already establishes effects of event order on resident diversity and taxon losses. Only resident endpoints analyzed; no analysis of newcomer performance or growth. Raw counts are common-scale colony observations, not biomass or true extinction. No experimental protocols reproduced.

Source rows 197 of 198 planned, 33 cells. One missing observation matches original publication. Two cells have duplicate rep labels with distinct rows; no guessed reassignment or matching on rep. Affected contrast is descriptive only. Other contrasts use approximate independent-group intervals conditional on the original physical-unit interpretation, not exact randomized inference.

| Comparison | Gross absolute response | Absolute net response | Pooled cancellation fraction | Opposing signs simultaneously supported | Cross-fitted pooled signed gap |
|---|---:|---:|---:|---|---:|
| primary | 35.1667 | 35.1333 | 0.0948% | false | -9.7931 |
| secondary 1 | 74.7500 | 74.7500 | 0.0000% | false | -0.0446 |
| secondary 2 | 29.6250 | 29.2917 | 1.1252% | not scored | not scored |
| secondary 3 | 39.2083 | 39.2083 | 0.0000% | false | -1.7396 |

All pooled cancellation fractions are small (0–1.1252%); the only comparison exceeding 1% has ambiguous replicate labels and is descriptive. No simultaneously supported positive and negative resident mean responses in any of the three inferentially scored comparisons, after Bonferroni approximate Welch intervals (five primary tests and ten secondary tests). All three pooled cross-fitted signed gaps are nonpositive. Within-level cancellation can exceed pooled cancellation, and one primary within-level cross-fitted gap is positive, but this is descriptive and not evidence for a general compensatory mechanism. Pooled metrics can obscure within-level heterogeneity. Repeated splits not independent, and nonpositive crossfit is not a test of equivalence or absence of antagonistic responses. Aggregate arithmetic cannot identify direct competition, facilitation, temporal interactions, or causal compensation. No exact randomization claim.

Prespecified exploratory cancellation test does not support a general opposing-response mechanism in this source. Preserve as a boundary condition rather than selecting favorable subgroup responses. All resident-only outcomes, source labels, approximate intervals and descriptive crossfit recorded in validation/ahif_v1/lear2024_resident_cancellation_v1.json. Separate source experiment, existing DB v0.7 effect and prediction tables unchanged until independent ledger integration/validation. No main or frozen manuscript changes.
