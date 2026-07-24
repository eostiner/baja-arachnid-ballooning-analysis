# Baja Ballooning publication QC report

Run completed from frozen inputs on 2026-07-23T19:49:16.683900+00:00.
- Incidence genera: 267
- Occupied 25-km cells: 205
- Equal cells sampled per band: 22
- Paired Monte Carlo iterations: 5,000
- Trait evidence counts: {'D1': 34, 'D2': 34, 'D3': 19, 'D4': 40, 'N0': 140}

## Interpretation rule
Jaccard total dissimilarity is partitioned into Baselga Jaccard turnover/replacement plus Jaccard nestedness-resultant. Simpson replacement is the turnover component of the Sørensen-family partition; it is related to, but not numerically identical to, Baselga Jaccard turnover.

## Files to inspect first
1. `00_INPUT_PROVENANCE.json`
2. `PUBLICATION_NUMBERS_ADJACENT_BANDS.csv`
3. `06_TRAIT_CONTRAST_SUMMARY.csv`
4. `figures/Figure_QC_02_Baselga_Jaccard_partition.png`
5. `PUBLICATION_CAPTION_DRAFT.txt`
6. iNEXT outputs added by `02_run_iNEXT_hill.R` when R is available.
