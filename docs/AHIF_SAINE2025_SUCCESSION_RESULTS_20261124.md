# Saine 2025: independent fungal introduction reanalysis — source verified

Source article https://doi.org/10.1093/ismeco/ycaf166 and exact public data release https://doi.org/10.5281/zenodo.13819389. Successful run https://github.com/zuizui0223/knk/actions/runs/45385872568, job 138416492000, commit eb3f43a3f3b242705b192b2bd876940f34f90a50. Three synthetic unit tests passed. RData MD5 4dbc27f81ffdcf1925887c7e0805b031. Original S2 taxonomy labels verified exactly. Source 4662 OTUs minus 9 unique introduced OTU labels = 4653 fixed columns; 7 columns have zero observations after author read-depth QC, yielding exactly 4646 observed non-target OTUs. All sample and OTU names align one-to-one. 1380 sample states, 11 low-read and 5 zero-nontarget states (overlap possible). 460 independently assigned physical logs; 446 complete triplets (347 treated, 99 controls); 14 excluded logs. All eligibility tests retain original 460 assignments; both Holm p=.432. Missingness gate passed but does not prove treatment-independent censoring. Same source log experiment as Ecology 10.1002/ecy.70013 and earlier baseline paper; do not double count. No raw data redistribution.

## Full endpoint ledger (no selection)

| Lane | Year | Endpoint | Partial R2 | Null-centered E | raw p | Holm p | complete logs |
|---|---:|---|---:|---:|---:|---:|---:|
| ten_group | 2020 | named | 0.021471 | 0.000983 | 0.2575 | 1.0000 | 446 |
| ten_group | 2020 | rank | 0.010140 | -0.008842 | 0.9375 | 1.0000 | 446 |
| ten_group | 2021 | named | 0.021335 | 0.000841 | 0.3010 | 1.0000 | 446 |
| ten_group | 2021 | rank | 0.023385 | 0.004230 | 0.3145 | 1.0000 | 446 |
| any_introduction | 2020 | named | 0.003372 | 0.001133 | 0.0345 | 0.1035 | 446 |
| any_introduction | 2020 | rank | 0.001185 | -0.000133 | 0.5385 | 0.5385 | 446 |
| any_introduction | 2021 | named | 0.004603 | 0.002368 | 0.0035 | 0.0140 | 446 |
| any_introduction | 2021 | rank | 0.004223 | 0.002837 | 0.1600 | 0.3200 | 446 |
| identity_among_introduced | 2020 | named | 0.023071 | -0.000713 | 0.6545 | 1.0000 | 347 |
| identity_among_introduced | 2020 | rank | 0.011616 | -0.010629 | 0.9285 | 1.0000 | 347 |
| identity_among_introduced | 2021 | named | 0.022264 | -0.001441 | 0.7765 | 1.0000 | 347 |
| identity_among_introduced | 2021 | rank | 0.024836 | 0.002288 | 0.4070 | 1.0000 | 347 |
| ten_group_scalar_secondary | 2020 | richness | 0.018426 | -0.001092 | 0.5440 | 1.0000 | 446 |
| ten_group_scalar_secondary | 2020 | simpson | 0.013081 | -0.005774 | 0.7695 | 1.0000 | 446 |
| ten_group_scalar_secondary | 2020 | top1 | 0.007536 | -0.011182 | 0.9530 | 1.0000 | 446 |
| ten_group_scalar_secondary | 2021 | richness | 0.013024 | -0.006488 | 0.7865 | 1.0000 | 446 |
| ten_group_scalar_secondary | 2021 | simpson | 0.018021 | -0.000905 | 0.5375 | 1.0000 | 446 |
| ten_group_scalar_secondary | 2021 | top1 | 0.023213 | 0.004215 | 0.3360 | 1.0000 | 446 |
| baseline_placebo | 2019 | named | 0.017402 | -0.003222 | 0.9725 | 0.9725 | 446 |
| baseline_placebo | 2019 | rank | 0.025484 | 0.006302 | 0.2325 | 0.4650 | 446 |

Only two-year any-introduction labeled community change survives the separate four-test Holm family (E=.002368, p=.014). This is small and not a universal effect. Under this prespecified abundance-change statistic, no evidence for identity-free rank shifts, nine-way introduced identity differences, or richness/concentration/top1 omnibus responses. Failure to reject is not equivalence; original HMSC presence-absence and phylogenetic results address different estimands, not contradicted by these nulls. Our primary 10-group family was entirely non-significant. One significant binary secondary lane must not be presented as a primary omnibus success, and family-wise correction is not global over all analyses. Relative read proportions are not absolute biomass, occupancy changes, or direct causal interaction coefficients. Conditional QC sample may be selected by post-treatment process. Site and log context not randomized. All 20 results registered (18 post-introduction and 2 pre-treatment placebo), 2 coverage tests stay QA only.

## Prior-art and next predictive check

Original Table 1 assigns treatment x year 2020/2021 plus treatment x log type shares of explained variance: presence-absence 2.9%, OTU richness 15.5%, DNA amount 6.9%, community composition 2.2% (sum of published component percentages, not comparable overall R2 or AHIF E). Original study already differentiated these response facets. Depth accounts for 24% of explained richness variance in its model; do not conflate observed OTU richness with true richness. Need read-depth-standardized sensitivity and out-of-site predictive comparison, not stronger causal interpretation. Source publication: https://doi.org/10.1093/ismeco/ycaf166.

Future prescore priority: fixed depth sensitivity on non-target reads, report loss of eligible logs and treatment dependence; no outcome-driven depth choice. For cross-site prediction, train first-year taxon response features only on other sites and use baseline composition to construct predicted second-year rank change; nested site cross-fitting, compare to baseline richness/concentration/log context/treatment identity, report each site's held-out loss. Distinguish predictive response features from independently measured competitive coefficients. Do not use the same held-out outcomes to estimate the effects that predict them. Five sites only.

Current study cluster is one additional independent experiment relative to Swedish grassland. AHIF v0.5 counts 33 study records / 35 panels / 161 effects / 11 source-audited records; no formal pooling or causal bridge. Existing frozen main and manuscript unchanged.
