# Canonical pipeline order

| Order | Step | Purpose |
|---:|:---|:---|
| 1 | 10 | Build aligned 25-km genus-by-cell incidence matrices and mutually exclusive latitude bands. |
| 2 | 10T | Normalize the authoritative CSV/workbook `exclusive_tier` + `primary_C3_group` schema and validate all matrix genera. |
| 3 | 11 | Foundational richness, completeness, and beta-diversity diagnostics. |
| 4 | 11B | Equal-cell total richness and Jaccard/Simpson/nestedness. |
| 5 | 11F | Equal-cell mean expected genus rarefaction. |
| 6 | 11G | Paired C3-versus-fixed-N0 richness and turnover; D4 excluded. |
| 7 | 11E | Four-boundary Jaccard turnover map. |
| 8 | 10A-10B | Ecoregion layer audit and cell crosswalk. |
| 9 | 10C-10E | Ecoregion richness, boundary turnover, and a priori test latitudes. |
| 10 | 10G | Latitude-band C3:N0 ratios and adjacent-band tests. |
| 11 | 10I | Integrated centerpiece map, Step 10I v7. |
| 12 | 12A-12C | Environmental audit, raster installation, and C3/N0 cell table. |
| 13 | 12D-12E | GLM screening, reduced models, and leave-one-band-out validation. |
| 14 | 12F | Environmental score preparation independent of GDM. |
| 15 | 12K | Spatial+ deconfounding and validation figure. |
| 16 | 12K1 | Sampled-cell Spatial+ map. |
| optional | 12J | Exploratory supplementary GDM only. |
| optional | 12N | Exploratory OMI environmental-niche position and breadth analysis; C3 versus fixed N0 with occupancy and taxonomy controls. |
| optional | 11H | Separate C3 and fixed-N0 geographic grouping, clustering, PCoA, Mantel comparison, and Baja map. |

<!-- PHASE13_PIPELINE_START -->
## Supplementary Phase 13 extension

| Order | Step | Purpose |
|---:|:---|:---|
| 17 | 13A | Audit and lock retained incidence, trait, environmental, and boundary inputs. |
| 18 | 13B | Construct independently frozen historical-boundary signals without outcome testing. |
| 19 | 13C | Construct the frozen domain-balanced contemporary-environment signal. |
| 20 | 13D | Build paired C3 and N0 Jaccard and Simpson dissimilarities. |
| 21 | 13E | Join historical, environmental, geographic, and community signals. |
| 22 | 13F | Primary paired historical-versus-contemporary cell-label permutation test. |
| supplementary | 13G | Pre-specified robustness and sensitivity analyses. |
| exploratory | 13H1–13H3 | Post-hoc Mulegé/ecoregion-junction follow-up. |
<!-- PHASE13_PIPELINE_END -->

<!-- PHASE14_15_START -->
## Supplementary Phase 14 temporal H3 extension

| Order | Step | Purpose |
|---:|:---|:---|
| 23 | 14A | Audit dated occurrence coverage and freeze eligible repeated cell-period comparisons. |
| 24 | 14B | Estimate paired C3 and N0 temporal Jaccard and Simpson turnover with equal-event resampling. |
| 25 | 14C0–14C1 | Prepare eligible cells and extract annual ERA5-Land and MODIS stress anomalies. |
| 26 | 14C–14D | Join temporal stress changes and run the original exploratory H3 model. |
| supplementary | 14E0–14F | Audit ECOSTRESS coverage and test prespecified extended heat, drought, water-balance, and vegetation sensitivities. |
| supplementary | 14G–14H | Test sampling-continuity and C1/C2/C3 trait-threshold sensitivities. |
| 27 | 14I | Synthesize all Phase 14 evidence without selecting the strongest-looking predictor. |

## Supplementary Phase 15 Bayesian H3 extension

| Order | Step | Purpose |
|---:|:---|:---|
| 28 | 15A | Freeze the original Phase 14 pair-level response, stress predictor, and observation uncertainty. |
| 29 | 15B | Fit the primary robust Bayesian paired hierarchical measurement-error model. |
| supplementary | 15C | Audit zero-centered prior and observation-uncertainty sensitivity. |
| supplementary | 15D | Compare stress and null models by leave-one-cell-out predictive CRPS. |
| supplementary | 15E | Run design-matched Monte Carlo power simulations with exact cell sign-flip inference. |
| 30 | 15F | Apply frozen synthesis rules and manuscript guardrails. |
| 31 | 15G | Generate the frozen spatial-temporal synthesis figure and editable-text SVG. |
<!-- PHASE14_15_END -->
