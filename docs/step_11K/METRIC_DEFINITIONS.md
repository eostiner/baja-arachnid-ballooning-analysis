# Metric definitions and interpretation

For two genus-incidence assemblages, let **a** be the number of shared genera, **b** the number unique to the first assemblage, and **c** the number unique to the second assemblage.

## Total Jaccard dissimilarity

`beta_jac = (b + c) / (a + b + c)`

This measures total compositional difference. It does not distinguish balanced replacement from richness imbalance.

## Baselga Jaccard replacement–nestedness partition

`beta_jtu = 2 min(b,c) / (a + 2 min(b,c))`

`beta_jne = beta_jac - beta_jtu`

- `beta_jtu` is the Jaccard-family turnover/replacement component.
- `beta_jne` is the Jaccard nestedness-resultant component.
- By construction, `beta_jac = beta_jtu + beta_jne`.

## Sørensen/Simpson replacement–nestedness partition

`beta_sor = (b + c) / (2a + b + c)`

`beta_sim = min(b,c) / (a + min(b,c))`

`beta_sne = beta_sor - beta_sim`

- `beta_sim` is Simpson replacement/turnover.
- `beta_sne` is Sørensen nestedness-resultant.
- By construction, `beta_sor = beta_sim + beta_sne`.

The Simpson replacement value is related to, but not numerically identical to, Baselga's Jaccard turnover component. Both are reported so that the manuscript can use precise terminology.

## Paired trait contrasts

Each Monte Carlo iteration draws the same occupied 25-km cells for ballooning-capable and non-ballooning genera. Trait contrasts are therefore computed as:

`ballooning-capable metric - non-ballooning metric`

A negative replacement contrast means less balanced genus replacement among ballooning-capable genera. A near-zero total-Jaccard contrast means that total compositional difference is similar between the groups, even when its replacement and nestedness components differ.

## iNEXT

The iNEXT stage is a separate coverage-standardized Hill-diversity diagnostic for q = 0, 1, and 2. It uses incidence-frequency vectors derived from the same frozen genus-by-cell matrix. It should not be described as the Baselga partition.
