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
