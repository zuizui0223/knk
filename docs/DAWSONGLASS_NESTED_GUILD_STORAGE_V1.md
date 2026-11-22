# Hierarchical assembly-history signal — post hoc guild audit

**Status: exploratory, same Dawson-Glass 24-pots × 15-species experiment, NOT a new independent experiment.** The original KNK v28–v38 manuscript and source fingerprints must remain frozen.

## Source / taxonomy

Dawson-Glass et al. (2025), *Ecology*, DOI https://doi.org/10.1002/ecy.4504; Dryad https://doi.org/10.5061/dryad.sbcc2frh4. User-supplied original assembly CSV SHA256: `1514d4501a0bdf7a61ace6c250566c5a0b8e4a24a0f4737646d182ebad7b8230`. 24 whole pots, 12 per treatment, 15 species, 360 observations.

Functional grouping is botanically derived from **5 Poaceae grasses** (Andropogon, Bromus, Dactylis, Sorgastrum, Tridens) versus **10 broadleaf forbs**. The raw CSV uses the spelling `Sorgastrum nutans`; retain source spellings for all data joins. The separate Dryad `species_info.csv` Functional column was **not available for independent exact-value comparison**. The original paper uses the same grass/forb distinction, and botanical family membership is externally verifiable. Do not claim official trait-file-level verification.

## New nested representation result

Each taxon-by-pot biomass vector was evaluated with the same 9,999 whole-pot treatment-label permutations (seed 20261122), Bray–Curtis distance-based R², and `E=R²observed−median(R²permuted)`. In guild-specific analyses, each pot's species biomasses were closed **within that guild** before rank sorting. Holm-corrected p-values across **all eleven post hoc targets** are shown.

| Target representation | E | p | post hoc Holm p |
|---|---:|---:|---:|
| All 15 species, labeled | +0.3959 | .0001 | .0011 |
| All 15 species, unlabeled sorted | −0.0018 | .5256 | 1.0000 |
| Grass vs forb relative biomass shares | +0.5004 | .0001 | .0011 |
| All 10 forbs, unlabeled within-forb ranks | +0.6177 | .0001 | .0011 |
| All 5 grasses, unlabeled within-grass ranks | +0.0855 | .1135 | .5675 |
| Forbs excluding Rumex: 9 species, labeled within-forb | +0.3870 | .0001 | .0011 |
| Same nine forbs, unlabeled within-forb ranks | +0.0566 | .1142 | .5675 |
| Grass vs forb shares after deleting Rumex | −0.0219 | .9037 | 1.0000 |

**Biological result:** the full-community sorted rank profile contains a weak measured history signal, whereas sorting the within-forb 10-species profile produces a strong, dominance-centered signal. Removing Rumex weakens that within-forb rank signal, yet the identities of the other nine forbs still differ by treatment. The grass-only relative internal profile does not show similarly strong evidence. This is **nested, scale-dependent visibility**, not an information increase when names are removed.

Forb share: ambient **26.18%**, warmed-order **52.15%**, mean pot-level difference **+25.97 percentage points** (stratified whole-pot bootstrap 95% interval **+15.82 to +35.79 pp**, 20,000 replicates). Within-forb maximum share: ambient **49.67%**, warm **80.33%**. After removing Rumex and re-closing remaining biomass, forb9 share **17.74% vs 17.27%**. Grass mean mass **42.55→25.28 g/pot** and forb mean mass **13.90→26.72 g/pot**.

## Mapping placebo (essential falsification guard)

Since v33 previously showed that seemingly biological mapping results can be artifacts, every one of **C(15,5)=3,003** equal-sized partitions was enumerated. For each, recompute treatment R² of the **sorted, internally normalized 10-species side**. Real forb10 observed R²=**0.6415095**, ranking at or above **22 of 3,003 partitions (0.73%)**.

However, **conditional on placing Bromus and Dactylis in the five-species group, and Rumex in the other ten**, the real grouping ranks at or above **22 of 220 (10%)**. Thus the apparent functional-group specificity is **largely traceable to separating major dominant species**. The actual two-guild share R² is not unusual among arbitrary partitions (upper fraction **18.9%**). These mapping fractions are descriptive diagnostics, not treatment-randomization p-values.

## What was already known vs what is new

The original source paper already demonstrated species-dependent biomass responses and the connection to advanced germination dates. Its functional-type covariate on *species-level treatment response* was not significant (p=.127); this is **not** a test of the two-guild *community total biomass shares*, which we analyze. No contradiction is implied.

The added KNK discovery is a **representation-specific decomposition** across nested levels of the same experiment, not a newly demonstrated causal effect of Poaceae-vs-forb physiology. All tests are post hoc and reuse the same 24 pots. Distance-based R²/E values after renormalization cannot be added or read as percentages of ecological information. Non-significant rank tests are not equivalence tests.

The source does **not** experimentally separate earlier planting/initial size and genuine biotic interactions. Guild biomass is not a direct measurement of ecosystem functioning or trait diversity. Mechanism remains unidentified; priority-effect explanations are hypotheses.

## Reproducibility

```bash
python benchmarks/audit_dawsonglass_nested_storage.py \
  --assembly-csv /path/to/Dawson-Glass_ECY24-0464_assembly_data.csv \
  --output /private/path/nested_guild_results.json
```

The script enforces the source hash, uses only the approved 24×15 pot table, checks guild mass, returns the same whole-pot permutation stream for all representations, adjusts 11 exploratory p-values, and enumerates 3,003 mappings. It does not commit/redistribute raw CSV.

Detailed aggregate values: `validation/dawsonglass_external_postsubmission/nested_guild_storage_summary_v1.json`. This development stays in unmerged **Draft PR #3**, separate from the frozen journal-submission version.
