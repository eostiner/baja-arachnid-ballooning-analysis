# Step 12N — OMI environmental niche analysis

## Biological question

Do ballooning-capable Baja arachnid genera occupy broader realized environmental niches than fixed non-ballooning genera?

This analysis addresses environmental niche position and breadth. It does **not** test the immediate weather conditions that trigger an individual ballooning event.

## Inputs

Step 12N uses retained 25-km pipeline products:

- Step 10 genus-by-cell incidence matrix;
- Step 10 taxonomy-strict incidence matrix, when available;
- Step 12C environmentally eligible occupied cells;
- Step 12F environmental predictor scores;
- reviewed D1–D4/N0 genus trait table;
- Step 10 cell polygons and genus taxonomy lookup.

## Primary analysis

- Grid: 25 km.
- Trait definition: C3 = D1 + D2 + D3 versus fixed N0.
- D4: excluded from the primary comparison.
- Minimum occupancy: five analyzed cells per genus.
- Environmental dimensions:
  - vapor-pressure deficit;
  - wind seasonality;
  - vegetation/phenology PC1;
  - topography PC1.
- Method: Outlying Mean Index analysis using `ade4::dudi.pca()`, `ade4::niche()`, and `ade4::niche.param()`.
- Global niche separation: Monte Carlo randomization test.
- Main trait test: `log(1 + total OMI tolerance)` modeled against ballooning classification while controlling for `log(1 + occupied cells)` and taxonomic order when estimable.

Total tolerance is calculated as `Tol + Rtol`. It represents realized environmental niche breadth. OMI is marginality: the distance between a genus's mean occupied environmental conditions and the average conditions available among analyzed cells.

## Sensitivity analyses

- C3 genera occurring in at least ten cells;
- C1, C2, and C4 ballooning definitions against the same fixed N0 reference;
- taxonomy-strict incidence matrix when available.

## Interpretation hierarchy

1. Use the occupancy- and order-adjusted C3 coefficient as the main trait result.
2. Treat raw boxplots and Wilcoxon tests as descriptive because niche breadth generally increases with the number of occupied cells.
3. Retain a main-text result only when the C3 direction is consistent at the five- and ten-cell thresholds and the uncertainty intervals support the same conclusion.
4. Otherwise report Step 12N as exploratory or supplementary.

## Outputs

The output folder is:

```text
PROJECT_ROOT/04_analysis/12N_omi_environmental_niche/
```

Key files:

```text
12N_omi_genus_niche_parameters_primary_C3.csv
12N_trait_effect_models.csv
12N_global_omi_randomization_tests.csv
12N_retention_screen.txt
12N_validation.csv
Figure_12N_caption.txt
figures/Figure_12N_OMI_niche_analysis_combined.{png,pdf,svg}
```

## Run

```bash
Rscript scripts/environment/step_12N/12N_install_packages.R

Rscript scripts/environment/step_12N/12N_omi_environmental_niche_analysis.R \
  "/Users/estiner/Desktop/OLD BALLOONING/Baja_Ballooning_Pipeline" \
  5000 20260717
```

## Pipeline position

Step 12N runs after Step 12F inputs exist. In the canonical full workflow it is placed after Step 12K1 so the retained Spatial+ analysis remains the primary environmental model and OMI is evaluated as a complementary environmental-niche analysis.
