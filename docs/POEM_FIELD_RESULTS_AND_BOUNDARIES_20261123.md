# POEM 2020/2021: real-data history-storage result and causal bounds (2026-11-23)

**Status:** independent real-data field study, SOURCE VERIFIED, exploration beyond the original publication, **not** a pooled meta-analysis. The existing KNK v28–v38 and Ecology Letters manuscript are unchanged. Original paper Alonso-Crespo et al. (2025), *Journal of Vegetation Science*, DOI https://doi.org/10.1111/jvs.70026; data Zenodo https://doi.org/10.5281/zenodo.10119982, version 1.0 (ZIP MD5 `3faa5f7d1619a3beae5941ba763b7a26`; SHA256 `7541ce47b61244f34a6b3b7f1187a2ef16f79bc173d73e6abed206203939ef8f`).

## Data restoration and quality

All **six original Excel biomass tables** were read in an ephemeral GitHub Actions runner. Two independent establishment cohorts (2020 / 2021), each **25 mixture plots**: five replicates of S (simultaneous), F (forbs first), G (grasses first), L (legumes first), or B (free succession). Ages 1–3 are repeated measures of the **same** 25 field plots within each cohort, not six experiments.

- 5,210 raw biomass rows in total; **62 repeated records for a plot×quadrat×identified-species cell** (mostly separately weighed fragments), summed in the same quadrat.
- 2,121 explicit zeros in the first 2020 harvest; later harvests mainly list species that were collected. Missing species within an otherwise sampled quadrat treated as zero harvested biomass (not as a missing study plot).
- Equal-area quadrats averaged within plot (four in 2020 age 1, two thereafter). `Unknown` excluded exactly as the published R script does; **49 named taxa**, of which **14 intentionally sown**.
- Primary F/G/L/S test excludes the B free-succession control, which is not a comparable seeded 14-species treatment: **20 whole plots per cohort per season**, 5 plots/treatment. Treatment labels permuted whole-plot within cohort (9,999, same permutation stream across all 3 ages).
- Protocol was committed **before experimental history-retention results were scored**: `docs/POEM_FIELD_PRE_SCORE_CONTRACT_20261123.md`. Script `benchmarks/score_poem_field_storage.py`. Original source data were not redistributed.

## Primary null-calibrated structural-retention estimates

The same 49 taxon biomass **proportions** were analyzed with and without identity names (the latter sorted by abundance per plot), using one-way Gower-centered Bray–Curtis R² and `E=R²observed−median(R²permuted)`. Both targets use the identical plots, dimensions, distance, and 9,999 within-cohort permutations. Six cell estimates are repeated states, **not** six independent meta-analysis effects.

| Cohort | Growing season | Named E | Named p | Rank-only E | Rank p | ΔE named−rank |
|---|---:|---:|---:|---:|---:|---:|
| POEM2020 | 1 | +0.0370 | .3063 | +0.1216 | .0677 | −0.0846 |
| POEM2020 | **2** | **+0.1872** | **.0038** | **+0.4270** | **.0022** | **−0.2398** |
| POEM2020 | 3 | +0.0749 | .0669 | +0.1269 | .0980 | −0.0520 |
| POEM2021 | 1 | +0.1465 | .0262 | +0.1519 | .0988 | −0.0054 |
| POEM2021 | **2** | **+0.3226** | **.0003** | +0.0230 | .3955 | **+0.2996** |
| POEM2021 | 3 | +0.0950 | .0430 | +0.1799 | .0463 | −0.0849 |

The two key second-season tests survive Holm multiplicity control across all **12 primary named/rank tests**: POEM2020 rank-only Holm p=.0242 and POEM2021 named Holm p=.0036. A difference in whether two p-values pass a threshold is **not itself** proof of a cohort difference; paired `ΔR²` and its within-plot bootstrap below supply a more appropriate descriptive contrast.

### Robustness to representation and treatment choice

- Alternate distance: **Hellinger–Euclidean**. In second season POEM2020, named E=+.1616, p=.0025; rank-only **E=+.3806, p=.0015**. In POEM2021, named **E=+.2861, p=.0001**; rank-only E=+.0214, p=.3745.
- F/G/L-only (15 original plots per cohort; synchronous S removed) preserves the main direction: POEM2020 age2 sorted E=+.3823, p=.0105; POEM2021 age2 named E=+.3782, p=.0001 and sorted E=+.0758, p=.2235.
- Only the **14 originally sown species**: POEM2020 age2 named E=+.1396, p=.0071 and rank E=+.1560, p=.1121; POEM2021 age2 named E=+.3499, p=.0004 and sorted E=−.0012, p=.5060. Importantly, the architecture-rich POEM2020 signal is less convincing when unsown species are omitted; this motivates a mechanism-facing hypothesis but **does not identify unsown recruitment as the cause**.
- POEM2021 age1 sown-14 only has two completely zero-total sown-species plots, and is correctly STOPPED rather than arbitrarily removing them.

### A plot-level longitudinal sensitivity (prevents temporal pseudo-replication)

Each plot's three annual species-proportion vectors were concatenated into **one 147-dimensional plot-level trajectory**, versus the matched three within-year sorted vectors. This exploratory *one plot, one response* score is:

| Cohort | Named E (p) | Rank trajectory E (p) |
|---|---|---|
| POEM2020 | +.1286 (.0011) | **+.2372 (.0004)** |
| POEM2021 | **+.1928 (.0002)** | +.0888 (.0480) |

Thus the cohort-specific emphasis persists when three seasons are kept together as a repeated trajectory. POEM2021 **does have** a marginal rank-trajectory association; do not call it universally devoid of architecture history.

### Second-season source-year contrast (exploratory, not a climate experiment)

Bootstrap treatment-stratified plot resampling: 3,000 replicates, fixed source groups. In season two:
- POEM2020 `R²named−R²sorted=−.2264`, 95% within-cohort bootstrap interval `[−.3370,−.0324]`.
- POEM2021 `R²named−R²sorted=+.3044`, interval `[−.0080,+.4979]`.
- Difference `(2021−2020)=+.5308`, 95% interval `[+.1331,+.7336]`.

This supports a **descriptive difference between the two tested establishment cohorts**, not a general causal weather effect. There is **only one starting field experiment for each weather year** and the cohorts were not randomly assigned to climate years.

### Biological context: dominance intensity versus winning identity

For 2020 at age2, the mean maximal-species biomass share by S/F/G/L was **.815/.422/.792/.530**; group G and S were especially dominance-heavy, often dominated by *Bromus hordeaceus*. For 2021 at age2, the same shares were **.512/.488/.606/.475**, much less spread; forb-first F frequently had *Jasione montana* as the top species (3 of 5 plots), while G was dominated by *Bromus hordeaceus* in all 5 plots. This is consistent with a switch in **statistical storage signature**, not conclusive evidence of a competition mechanism.

The deliberately post-hoc unsown-35-only check was not decisive: in POEM2020 age2 sorted E=.1369 (p=.1128), and POEM2021 has a plot with zero unsown biomass, so unsown-only endpoint STOPPED. The full-community result could be driven by shifts in composition, the ratio between sown and unsown species and geometry; **none is uniquely identified here**.

## Meta-analysis / theoretical conclusion

POEM is an additional **independently conducted field experiment** (distinct from the earlier rhizobox), but all its cohorts and seasons belong to **one publication cluster**. We now have evidence from **four independent published datasets** (fungi1, rhizobox grassland1, Dawson-Glass grassland1, POEM field grassland1), but they have different factor definitions and source-specified endpoints; publication/lab dependence remains. A cross-study pooled effect would be premature. The strong empirical advance is a within-site illustration of **cohort- and observation-dependent history retention**; the proposed universal macroecological law remains an untested hypothesis.

### Files & verification

- Prescored protocol: `docs/POEM_FIELD_PRE_SCORE_CONTRACT_20261123.md`
- Real-data code: `benchmarks/score_poem_field_storage.py`
- Repeated-plot and alternative-metric sensitivity: `benchmarks/audit_poem_history_channel_switch.py`
- Source-scored aggregate: `validation/macro_meta/poem_field_cohort_age_results_v1.json`
- Post-hoc sensitivity: `validation/macro_meta/poem_temporal_channel_sensitivity_v1.json`
- Source+score Actions run https://github.com/zuizui0223/knk/actions/runs/45085174472 (success)
- Follow-up Actions run https://github.com/zuizui0223/knk/actions/runs/45085506739 (success)
- Source ZIP never committed; under Draft PR #4, separate from current frozen journal manuscript.

**Claim boundary:** Do not say stronger or weaker p directly measures more or less information; no causal competition, climate-year effect, or root function has been identified. Null-calibrated E is geometry and design dependent.

