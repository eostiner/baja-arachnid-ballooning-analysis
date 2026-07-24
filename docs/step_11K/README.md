# Step 11K–11N publication diversity and turnover extension

## Step 11K

Step 11K provides the publication quality-control extension for adjacent
latitude-band diversity and turnover. It uses the frozen genus-by-25-km-cell
matrix, 5,000 paired equal-cell iterations, and the reviewed D1–D4/N0 trait
classification. The primary comparison is D1–D3 versus N0, with D4 excluded.

The analysis reports total Jaccard dissimilarity, its additive Baselga
replacement and nestedness-resultant components, the corresponding
Sørensen-family sensitivity metrics, and coverage-standardized iNEXT Hill
diversity at q = 0, 1, and 2.

## Step 11N main figure

Step 11N generates the main map plus inset figure:

- Panel A: within-band q = 0 richness, q = 1 common-genera diversity, and q = 2
  dominant-genera diversity, shown as ballooning-capable / non-ballooning;
- Panel B: adjacent-boundary ballooning-minus-non-ballooning contrasts in
  Jaccard replacement and Jaccard nestedness-resultant.

Panel B is deliberately not q-specific. Boundary-level q values belong to Hill
beta diversity and are reported separately in the Step 11K supplementary
outputs. Baselga replacement/nestedness and Hill beta answer different
questions and must not be presented as interchangeable measures.

## Compact archived outputs

The `results/` directory contains the selected summaries needed to reproduce
reported values. The multi-million-row raw iteration tables and fitted iNEXT
RDS object remain locally generated rather than committed.
