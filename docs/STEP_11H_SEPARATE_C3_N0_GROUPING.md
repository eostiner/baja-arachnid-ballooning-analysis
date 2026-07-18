# Step 11H — Separate C3 and N0 geographic grouping

## Status

**Optional exploratory / supplementary analysis.**

Step 11H does not replace the paired Step 11G C3-versus-fixed-N0 turnover
analysis. It analyzes ballooning-capable and non-ballooning assemblages
independently to ask whether latitude bands and ecoregions group in the same way.

## Trait definition

- **C3 ballooning:** D1 + D2 + D3.
- **N0 non-ballooning:** fixed non-ballooning reference.
- **D4:** excluded.

## Methods

- 5,000 equal-cell Monte Carlo resamples.
- Jaccard dissimilarity as the primary composition metric.
- Simpson replacement as a sensitivity metric.
- Average-linkage hierarchical clustering.
- Classical PCoA.
- Resampling stability for closest regional pairs and most-distinct regions.
- Spearman Mantel permutation comparison of C3 and N0 distance structures.
- Two-panel Baja map using the validated mainland outline.

## July 2026 full-arachnid result

- C3: 87 genera; N0: 140 genera.
- C3 most often grouped 23–24°N with 24–26°N (56.4%).
- N0 most often grouped 24–26°N with 26–28°N (47.1%).
- C3 most often isolated 28–30°N (61.6%).
- N0 strongly isolated 30–32°N (86.6%).
- Latitude distance-structure concordance: Mantel rho = 0.430, p = 0.2110.
- Ecoregion distance-structure concordance: Mantel rho = 0.503, p = 0.0023.

Interpret cautiously: broad regional structure is shared, but the strongest
expression of compositional isolation differs between C3 and N0.

## Run

```bash
python run_pipeline.py \
  --project-root /path/to/Baja_Ballooning_Pipeline \
  --from-step 12K1 --to-step 12K1 --skip-step 12K1 \
  --include-grouping
```

Or directly:

```bash
python scripts/biogeography/step_11H/11H_separate_C3_N0_grouping.py \
  --project-root /path/to/Baja_Ballooning_Pipeline \
  --iterations 5000 \
  --mantel-permutations 9999 \
  --seed 20260717
```

## Interpretation boundary

Step 11H is exploratory. Differences in the most frequently isolated band do
not by themselves prove that the complete C3 and N0 spatial structures differ.
The Mantel comparison and primary paired Step 11G analysis remain essential.
