# AHIF-DB: original-review 268-study audit and scalable-interaction novelty gate

**2026-11-24.** This is a direct, hash-verified re-analysis of **the existing review's bibliographic metadata**, NOT a new ecological experiment, systematic re-screen of the original papers, or meta-analysis. Source DOI: https://doi.org/10.1016/j.tree.2023.08.001; original GitHub repository [hengxingzou/Zou2023TrendsEcolEvol](https://github.com/hengxingzou/Zou2023TrendsEcolEvol), `SupportingData.xlsx` source Git blob SHA `be2777d4e4bb1454db393b94a2118ada79118dfa` and original `LitReview.R` lines 45–53. Source has **268 study rows ×14 columns**; all 268 retained as `not_independently_screened` in AHIF's separate bibliography table.

## Source-code reproducibility: one mis-binned paper

The review's original R `case_when` bins the `# Manipulated Species` column by string membership in [1,2], [3..10] or literal `">10"`, and **all other codes are put into the no-temporal-treatment bin**. One original record (`ZOU2023_R0032`, Fukami & Morin 2003, *Nature*, *Productivity-biodiversity relationships depend on the history of community assembly*) uses the numeric code `18`, which is therefore incorrectly routed to `No temporal treatment` even though its `Temporal Gradient?` is `Yes`.

| Source design bin | Exact author's R code | Recode numbers as quantities |
|---|---:|---:|
| No temporal treatment | 46 | **45** |
| 1–2 manipulated taxa | 107 | 107 |
| 3–10 manipulated taxa | 72 | 72 |
| 11+ manipulated taxa | 43 (literal `>10` only) | **44** (plus numeric 18) |
| All | **268** | **268** |

This is a **metadata visualization/classification correction**. It does not establish any change in the original scientific claims or license us to replace author classifications in the original dataset. AHIF stores both original and corrected bins, source ID and audit explanation.

## Relevant metadata-only patterns

- Total review studies: **268**, years 1983–2022. These are **not independent AHIF-curated paper reviews**.
- Temporal gradient flagged by original authors: **58/268**. This is not the fraction with any arrival-time treatment: some categorical treatment designs have no continuous gradient.
- By original `Fitted Competition Model` field, **13/268** articles are assigned a specific non-`-` model code (e.g. Beverton–Holt, consumer-resource, niche-based, co-infection, plant–soil feedback). **255 have literal `-`**, which must NOT be interpreted as an absence of measured interactions or competition.
- Of those **13 model-coded reports**, three are author-coded **no temporal treatment**; the other **10 are in 1–2 taxa (6) or 3–10 taxa (4)**. Among **44 recoded 11+ manipulated-species reports**, **0 have a model code in this particular review field**. This is a credible *review-screening workload signal*, not a global proof that quantitative models are unknown for richer communities.
- Original generation-class labels: **<=1 generations 109**, **2–5 generations 27**, **>5 generations 121**, Mixed 10, Not Applicable 1. These bins are author-supplied categorical descriptors, not a re-estimated generational scale.
- There are 21 **nearby prior-art bibliography candidates** generated reproducibly as the union of 11 original-title high-risk cases (community history × function/resource/spatial preemption) and all 13 source-model-coded cases, with three overlaps. All 21 remain unreviewed for AHIF-specific effect eligibility.
- A non-exhaustive 2024–2026 research update is **mandatory** before any global rarity or first-in-literature claim. An old 2023 source review is not enough.

### Why this changes the novelty target

The **established** contributions of prior art already include:
1. Arrival history affecting community composition and ecosystem function, e.g. Fukami et al. 2010 (DOI 10.1111/j.1461-0248.2010.01465.x) and Leopold et al. 2017 (DOI 10.1111/ele.12803).
2. Measured species competition predicting complex community composition and historical priority patterns, e.g. Uricchio et al. 2019 (DOI 10.1086/701434).
3. Species arrival-dependent interspecific effects, coexistence and phenology/generation changes, e.g. Grainger et al. 2019 (DOI 10.1073/pnas.1803122116), Fragata et al. 2022, van Steijn et al. 2025 and Zou & Rudolf 2026 (DOI 10.1002/ecy.70454).
4. Multispecies priority assembly graphs already described theoretically and empirically in Song et al. 2021 (DOI 10.1111/ele.13870).

AHIF must **not claim any of these concepts as new**. The testable extension is **scale-and-representation aware interaction-to-structure transmission**: *Does experimentally measured history-dependent variation in partner effects predict named-species versus identity-erased abundance-rank memory in rich multispecies communities, and generalize to unseen communities?*

This proposal requires source-level support, and currently no AHIF-curated experimental panel identifies the whole same-unit causal chain.

## Highest priority raw-data test: Leopold et al. 2017

The 2017 factorial fungal experiment is scientifically close enough to **falsify a weak novelty claim**. Original full methods report 10 introduced fungal taxa, 4 initial fungal species plus an uninoculated control, **2 nitrogen ×2 fungivore ×5 histories ×5 replicate jars ×2 destructive harvests = 200 assembled/control microcosms**. An additional **400 single-species control microcosms** were conducted. The original *Dryad* data DOI https://doi.org/10.5061/dryad.7p2cv publicly lists three input files:

- `sample.data.csv` — jar metadata and wood decomposition
- `species.prevalence.csv` — prevalence of each fungal taxon in **nine wood-disk subsamples per jar** (do NOT treat these nine subsamples as nine independent jars)
- `collembola.csv` — grazer abundance category per jar

**Source audit attempt 2026-11-24:** DOI metadata and 3 original file IDs were verified through Dryad REST v2 (version ID 13946; file IDs 52325, 52326, 52327), but the original CSV binary endpoints gave **HTTP 401 Unauthorized / HTTP 403 Forbidden** from the available GitHub Actions environment. Source run: https://github.com/zuizui0223/knk/actions/runs/45316446032. **Zero original CSV files were retrieved**; do not claim the KNK named-versus-rank rescoring is done. There is no defensible biological effect estimate from the source yet.

Once original source data become available through normal authorized access, prespecify the design-respecting, within-nitrogen×fungivore×harvest label permutation, using the **microcosm** as the resampling unit, and score named fungal prevalence vectors and the same 10 counts sorted without species names. Separate direct-grazing intervention effects on resulting rank storage from fungal interaction coefficients (which would still require different controls or dynamics).

## Evidence hierarchy / reproducibility

- `review_candidates`: 268 **unreviewed** author-source bibliography entries.
- `review_priority`: 21 title/model-coded **unreviewed** high-priority candidates, with FK to the 268 source records.
- `review268_design_corrected_bins.csv`: per-candidate faithful original R bin + numerically corrected bin + source classification discrepancy.
- `review268_design_audit.json`: aggregate audited review descriptors, not ecological treatment effects.
- `studies`: 26 targeted AHIF-curated publication records including unreviewed paper-only screening leads.
- `panels`: 28 experimental panel records; 19 unscreened panels are **unknown**, not documented absence; 9 source-audited panels have no *identified* mechanism-to-structure causal path.
- `effects`: 45 numerical source-linked endpoints at original experimental grains (many post hoc) without random-effects meta pooling.
- `hypotheses`: 8 falsifiable research hypotheses, with H08 focused on transfer from 11+ taxa interactions to named/rank community memory.

**Do not call model coding counts the prevalence of ecological interactions**, original review R error a novel ecological finding, or the zero causal-bridge count proof that no existing paper has this design. This remains a targeted pilot and conservative prior-art challenge.
