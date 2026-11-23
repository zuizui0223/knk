# Source-data missingness amendment BEFORE outcome scoring, Jing Zenodo 2023
**2026-11-24.** Authenticated original Zenodo source MD5 5afb3a33ad49573aa7e255d2541ab0ed; this is a data-quality AMENDMENT to docs/JING_ZENODO_MIXTURE_SOIL_HISTORY_PRE_SCORE_20261124.md, made after inspecting missingness but **before inspecting any treatment effect scores**. Do not silently rewrite the locked initial contract.

## Verified original source structure
The XLSX mixture worksheet has two header rows, exactly **40 true experimental condition rows**, and **one trailing entirely blank padded row**. The latter is Excel formatting residue and is ignored, documented as one removed blank row.

Source-only missingness audit identified exactly four missing plant species measurements:
1. Wet, SK soil, replicate position 3, four-species SK-AF-AC-ES mixture: AC weight missing, other components present.
2. Wet, AF soil, replicate position 3, same four-species mixture: SK weight missing, other components present.
3. Dry, S sterile soil, replicate position 1, two-species AF-AC mixture: both AF and AC missing, **entire separate mixture pot unobserved**.
4. Dry, S sterile soil, replicate position 1, four-species SK-AF-AC-ES mixture: ES weight missing, other components present.

**Locked scoring rule:** primary analysis uses only three two-species mixtures (SK-AF, AF-AC, SK-AC) and five monoculture species controls. The missing **whole pot** in AF-AC/Dry/S contributes **n=3 independently measured mixture pots** instead of 4; it is excluded without replacement or zero imputation. The paired species within that pot are jointly absent and cannot be counted as two zero-biomass observations. All other primary two-species mixture cells n=4. Full five monoculture species×water×soil cells are n=4 each.

The three partially missing species observations only occur in a FOUR-species mixture. They are **not included in the scored primary analyses** and are not set to zero. If a later multispecies score is attempted, a new source-level missing-data plan will be required before looking at outcomes. Such mixtures must not be used in the current primary result.

Use 10,000 bootstrap draws with each source factor-cell's ACTUAL observed pot count (3 or 4), resampling both partner species together within a mixture pot and independently resampling mono-control species pot units. Keep the original shared-control and sterile-baseline bootstrap covariance logic.

Counts expected: original 2×5×5×4 =200 source monoculture pots (if one species per separate pot); 3×2×5×4=120 nominal primary pairwise mixture pots, with **119 observed**, one missing. There are 3 other 3-/4-species mixture types, not scored; total physically independent pot IDs remain inferred from factor × combination × row order and have not been verified against author records. No claim of 440 independently randomized physical units.

This amendment is solely source-grain repair, **not a post hoc effect selection or biological mechanism claim**. Treat this 2023 Zenodo and Jing 2024 Dryad as **potentially overlapping** sources, not two new independent meta studies. Main KNK v28–v38 remains frozen.
