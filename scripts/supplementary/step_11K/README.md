# Step 11K — publication nestedness/replacement QC and iNEXT

Step 11K reruns the latitude-band trait comparison from the frozen 25-km
incidence matrix and reviewed D1–D4/N0 trait table.

## Primary trait coding

- ballooning-capable: D1–D3;
- non-ballooning: N0;
- D4: excluded from the primary comparison.

## Core analyses

- 5,000 paired equal-cell Monte Carlo iterations by default;
- 22 occupied cells per latitude band;
- total Jaccard dissimilarity;
- Baselga Jaccard replacement and nestedness-resultant components;
- Sørensen-family Simpson replacement and nestedness-resultant sensitivity;
- ballooning-capable minus non-ballooning paired contrasts;
- incidence-based iNEXT at q = 0, 1, and 2, standardized to shared sample coverage.

The iNEXT Hill estimates are a diversity analysis and are **not equivalent** to
the Baselga replacement–nestedness partition.

## Run

```bash
bash scripts/supplementary/step_11K/run_step11K.sh   /path/to/Baja_Ballooning_Pipeline 5000 20260723 200
```

Outputs are written to:

```text
PROJECT_ROOT/04_analysis/11K_publication_nestedness_replacement_QC/
```

Large raw iteration tables are regenerated locally and are intentionally not
stored in GitHub. Compact summaries, provenance, iNEXT results, and QC figures
are archived in `docs/step_11K/`.
