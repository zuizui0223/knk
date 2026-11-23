# Jing 2024 soil-inoculum competition source-access gate (2026-11-24)

**Scientific status: PAPER AND FILE METADATA VERIFIED, RAW XLSX NOT OBTAINED. NO KNK EFFECT SCORED.**

Paper: Dabu/Xilatu, Jingying Jing, Ji, Yang & Bezemer (2024), *Functional Ecology*, DOI https://doi.org/10.1111/1365-2435.14631. Public deposit DOI https://doi.org/10.5061/dryad.jq2bvq8jj, July 26 2024, `Fulldata.xlsx` (~89.85 kB) and README (~4.53 kB). The full 11-sheet schema is described in the public Dryad dataset landing page, but the cell-level workbook was not available through our environment.

**Empirical result from original paper, NOT new KNK research:** three focal grassland taxa `Stipa krylovii`, `Artemisia frigida`, `Agropyron cristatum` in monoculture/mixtures, 5 soil inocula × water treatments. Soil source and drought changed relative yield and biotic feedback, **but the rank of competitive abilities was unchanged** (A. cristatum > S. krylovii > A. frigida). Source: https://datadryad.org/dataset/doi%3A10.5061/dryad.jq2bvq8jj.

**Authentic access attempted with public GitHub Actions** using Dryad's REST v2 endpoint, DOI URL encoding, public version 307544. The API's dataset/version/files metadata GET requests succeeded, but binary `/api/v2/files/{id}/download` returned **401 Unauthorized** and browser-style `/downloads/file_stream/{id}` returned an HTML `Validating...` challenge, rather than XLSX bytes. The code correctly rejects non-XLSX. Run logs: https://github.com/zuizui0223/knk/actions/runs/45245012338 . No authorization bypass, browser evasion, invented data, or guessed effect size is permitted. The absence of raw data blocks reconstruction of pot units and independent source reanalysis.

**Fallback screening only:** Zenodo https://doi.org/10.5281/zenodo.7856898 ("Drought alters soil microbiome-mediated plant coexistence"), another 2023 record from Jing/Yang/Da, offers public original `aboveground biomass-Drought alters soil microbiome-mediated plant coexistence.xlsx` with MD5 `5afb3a33ad49573aa7e255d2541ab0ed`. **Potential source-study overlap with 2024 Jing dataset is unverified**; do not count it as an independent paper or combine corresponding effects in a meta-analysis until original study, provenance and experimental treatment levels are matched.

**Next strict gates if bytes are acquired:** source hash, sheet grain, independent experimental pot IDs, replica vs author-derived rows, species × soil × drought × alone/mixture design, and reciprocal independent design controls. Only then a predeclared competition-difference-in-differences should be computed.

The development stays in KNK Draft PR #5, and original frozen v28–v38 and Ecology Letters manuscript remain untouched.
