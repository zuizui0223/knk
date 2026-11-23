# Jing/Yang/Da 2023 Zenodo soil-history × drought × mixture — analysis contract
**Frozen 2026-11-24, before KNK-specific biological scoring.** This is related to the Jing research group; potential overlap with 2024 Jing Dryad source 10.5061/dryad.jq2bvq8jj is unresolved. Do NOT count both as independent papers or experiments.

## Exact provenance and experimental grain
Zenodo original: https://doi.org/10.5281/zenodo.7856898. Workbook published 2023-04-23, MD5 5afb3a33ad49573aa7e255d2541ab0ed, SHA256 636ad7ebbe1f8c452b40fd8fde0e03cf611d21baf1db959fcce42385bc23b899. README MD5 ae5c5b22981990bdf3a701cd5c7d7393. GitHub Actions https://github.com/zuizui0223/knk/actions/runs/45245492910 verified source but did not score outcomes. Raw XLSX is never committed.

Original two sheets: monoculture_aboveground biomass = 40 rows and five species columns SK, AF, AC, PT, ES; mixture_aboveground biomass = 42 rows including two header rows, 24 columns with 6 assembly groups SK-AF, AF-AC, SK-AC, SK-AF-AC, SK-AF-AC-PT, SK-AF-AC-ES. Both contain Wet/Dry and inocula SK, AF, AC, MIX, S. Four reported replicates per water × soil × assembly. Original source README: aboveground biomass per species; not root biomass or actual competition coefficients. The row number is a reconstructed within-context replicate index, NOT a verified randomized greenhouse block. Mixture species in one pot must not be treated as separate experimental units. Potential count: 200 separate monoculture species pots and 240 mixture pots (subject to source methods), never only 40 total pots.

## Source-locked primary contrast
For each of three original TWO-species mixtures, both member taxa, five inocula, and Wet/Dry:
C(i,m,s,w) = log[ mean_aboveground_biomass_mixture(i,m,s,w) / mean_aboveground_biomass_monoculture(i,s,w) ].
C is SOLO-REFERENCED MIXTURE OUTPUT, not a per-capita interaction coefficient unless initial planting density is comparable (not yet verified).
D(i,m,s,w) = C(i,m,s,w) - C(i,m,sterile_S,w), for four inocula SK,AF,AC,MIX relative to S. Compare D(dry) vs D(wet) descriptively.

Report 3 pairs × 2 directed roles × 2 watering treatments × 5 inocula = 60 C contexts; 48 D comparisons. Pairwise order is only mean raw mixture biomass of its two species within a context, not general fitness rank. Explicitly report if that observed order flips across 10 inoculum×water contexts.

## Biological and statistical guardrails
- Validate full source MD5 and SHA, both original sheets, exact two-level mixture header, 2 watering levels ×5 soil types ×4 measured rows/cell, all expected species × mixture variables and finite nonnegative biomass.
- Zero individual biomass is allowed if group mean positive; a zero group mean triggers STOP instead of unannounced pseudocount.
- 10,000 BOOTSTRAP replicates, seed 20261124, resampling whole mixture pots with all component species together, and each species' independent monoculture pots within same water×soil. The same control draw is shared across all corresponding mixtures in a replicate. Report C, D, water interactions and 95% descriptive percentile intervals. Four pots/group = weak uncertainty; no pseudo-replication of species rows or 48 soil differences as independent studies.
- Soil inoculum sterile vs live contrasts can include nutrient/chemical as well as microbial changes. No isolation of the responsible mediator.
- No formal randomized permutation p values or causal alpha claims unless source methods give treatment-randomization blocks.
- Neither exact source paper nor experiment independence from 2024 Dryad is confirmed. This dataset is source-verified but held as a POSSIBLY OVERLAPPING research-group cluster.
- Not an experiment on direct species arrival order: it is soil-history and biotic-feedback mechanism lane, separate from KNK original priority sequence lane.
- KNK frozen v28–v38 and current Ecology Letters manuscript unchanged; outputs stored in unmerged Draft PR #5.
