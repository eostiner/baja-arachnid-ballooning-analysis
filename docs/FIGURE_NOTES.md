# Figure-generation notes

## Centerpiece: biogeographic dispersal-group balance

Step 10I v7 integrates the independently mapped ecoregions, equal-cell C3/N0 richness, ecoregion B:N ratios, latitude-band B:N ratios and subset-resampling intervals, and the four a priori test latitudes. It writes editable vector output as well as raster versions:

```text
PROJECT_ROOT/04_analysis/C3_pipeline_rebuild/09_C3_biogeographic_concordance/
10I_final_map_and_ratios/publication_outputs/
├── Figure_3_Biogeographic_Dispersal_Balance.pdf
├── Figure_3_Biogeographic_Dispersal_Balance.svg
├── Figure_3_Biogeographic_Dispersal_Balance.png
└── Figure_3_Biogeographic_Dispersal_Balance.jpg
```

For manuscript layout, use the PDF or SVG. The dashed horizontal lines are **a priori test latitudes**, not boundaries demonstrated by the analysis. Donut area represents mean expected C3+N0 richness from eight occupied cells; D4 is excluded. The band intervals are percentile distributions among equal-cell subsets of the observed occupied cells, not confidence intervals for total regional richness.

## Turnover figure

Step 11E displays all four adjacent-band Jaccard comparisons. Step 11B writes the corresponding Simpson replacement and nestedness summaries, and Step 11G writes the paired C3-versus-N0 turnover contrasts.

## Spatial+ figure

Step 12K writes the coefficient and validation figure. The Step 12K1 map is a sampled-cell descriptive visualization and should not be described as a suitability or peninsula-wide prediction map.
