# Dawson-Glass follow-up: winner replacement with nearly unchanged dominance architecture
**Status:** Exploratory extension after the frozen primary Dawson-Glass test. Same 24 pots and original source CSV, **not a new independent experimental replication**. No changes to KNK v28–v38 and not a new claim for the original Ecology Letters submission.

## Data and statistical scope
Original `Dawson-Glass_ECY24-0464_assembly_data.csv`, SHA-256 `1514d4501a0bdf7a61ace6c250566c5a0b8e4a24a0f4737646d182ebad7b8230`, 24 pots × 15 species, ambient/warm order 12 each. Raw user-supplied input is not redistributed.

The original publication **already** established that earlier warming-shifted germination correlates with higher biomass and that *Rumex crispus* advanced by 19 days. Therefore KNK does not claim the species response or germination-to-biomass association as a novel result. Independent original source: https://doi.org/10.1002/ecy.4504.

## New decomposition of the previously scored structure
- Labeled-composition R² = **0.428959**; identity-free sorted-rank R² = **0.032899**.
- Paired `R²_labeled−R²_sorted` = **+0.396060**. Under the **same 9,999 whole-pot label permutations** (seed 20261122) null median = −0.001380, central 95% = [−0.1013, +0.1151]; upper-tail **p=0.0001**.
- This is a treatment-label permutation contrast under the *global no-treatment association* null, not a conditional test that rank response is exactly zero. Nor is R² loss a fraction of mutual information.
- Bray–Curtis distance between treatment **mean species-labeled proportional profiles** = **0.3841**, but between treatment **mean sorted profiles** = **0.0515**. These distances describe a substantial winner change without comparable change to rank geometry, not a causal partition.

## Who wins the community?
A winner is the **unique most-biomass species in each whole pot** (no ties). Winner identity shifted:

| Most-abundant species | Ambient | Warm-order |
|---|---:|---:|
| *Bromus inermis* | 6 | 2 |
| *Dactylis glomerata* | 6 | 2 |
| *Rumex crispus* | 0 | 8 |

Global winner × treatment association: likelihood ratio G²=**15.27634**. **Exact conditional p=0.003883** by enumerating 61 possible 2×3 winner tables weighted by the 2,704,156 possible 12-of-24 treatment assignments. Monte Carlo whole-pot permutation gives p≈0.0034. This is a **post hoc global analysis** after seeing the primary results; it demonstrates an association, not an independent test of mechanism.

## Absolute biomass: potential compensation, not proven competition
Changes in treatment mean **grams per pot**, warm-order minus ambient:
- *R. crispus*: **+15.93 g** (20,000 whole-pot bootstrap 95% **+11.92 to +19.94 g**).
- *B. inermis* and *D. glomerata* together: **−18.95 g** (95% **−28.10 to −9.57 g**).
- Other species combined except *R. crispus*: **−20.37 g**.
- Total aboveground community biomass: **−4.44 g** (95% **−11.82 to +3.14 g**); the primary treatment test for total biomass was p=0.3284.
- Mean largest relative share: ambient **0.486**, warm order **0.466** (bootstrap CI for difference **−0.084 to +0.043**). This cannot establish equivalence or a precisely zero effect.

Thus there is an **observed group-level biomass compensation pattern** (one growing species and two shrinking grasses with similar overall community output). This is a candidate *ecological explanation* for why named community composition changes greatly without a comparable unlabeled-rank or total-biomass change. But the design cannot demonstrate a within-pot competitive tradeoff or niche pre-emption.

## Additional explanatory constraints
- **Initial height (measured immediately before planting):** *R. crispus* average **0.78 cm** in ambient order vs **2.21 cm** in warm order. This variable may be affected by schedule/age and is not an unconfounded pre-treatment covariate; adjustment cannot by itself identify competition.
- Original authors did not directly separate species arrival-time advantages from changed interaction strength. Greenhouse mesocosms were under common growth temperatures; “warm order” was based on separate warmed-germination trials.
- Individual-species contrasts using 9,999 whole-pot permutations and Benjamini–Hochberg across 15 taxa show share differences for *R. crispus*, *Solidago canadensis*, *D. glomerata*, *H. punctatum*, and *R. pennsylvanicus* (exploratory q<0.05). These unadjusted pot-level comparisons are **not** substitutes for the original authors' phylogenetic mixed model, so do not assert discrepancies represent new findings.
- No raw data or cross-system meta-analytic synthesis changed. Correct next step for **mechanistic identification** is an arrival-time × initial-size × competitor-presence factorial design, preferably with repeated absolute-abundance measurements.

## Reproduce
From the KNK repository root on the separate post-submission branch:

```bash
python -m benchmarks.run_dawsonglass_winner_identity_followup \
  --assembly-csv /local/path/Dawson-Glass_ECY24-0464_assembly_data.csv \
  --output /local/private_results/winner_followup.json
```

The script checks the original file SHA and source panel, and writes only aggregate results. See `validation/dawsonglass_external_postsubmission/winner_followup_summary_v1.json` for headline numbers. No release or merge is requested before source-method and manuscript author review.

