# Phase 13 — historical versus contemporary mechanisms

**Repository release:** 1.0.6
**Locked:** 2026-07-18

Phase 13 is a supplementary mechanistic extension of the retained C3 analysis. It asks whether lower compositional replacement of ballooning-capable assemblages is associated with stronger contemporary-environment tracking and/or weaker historical-boundary structure.

## Retained definition

- C3 = D1 + D2 + D3: **87 genera**
- Fixed N0 non-ballooning: **140 genera**
- D4 excluded: **40 genera**
- Total retained: **267 genera**
- Occupied 25-km cells: **205**
- Contemporary-environment eligible cells: **189**

The authoritative Phase 13D classification is the retained normalized Step 11G trait lookup. Phase 13 does not infer D1–D4 from older descriptive evidence fields.

## Inferential hierarchy

### Confirmatory: Step 13F

Primary response: paired C3 minus N0 Simpson replacement. Secondary response: paired Jaccard dissimilarity. The primary complete matrix contained 120 cells and 7,140 pairs. Inference used Freedman–Lane-style cell-label residual permutation with 4,999 permutations; pairwise rows were not treated as independent.

### Robustness: Step 13G

Pre-specified cell thresholds, boundary coding, environmental domains, and secondary contextual transitions were used to test robustness rather than search for significance.

### Exploratory/post hoc: Steps 13H1–13H3

The visually noticed Mulegé/central-Gulf feature and ecoregion-junction hypothesis were tested separately. These analyses cannot overturn 13F/13G.

## Locked results

C3 assemblages had lower raw mean Simpson replacement than N0 assemblages (approximately 0.381 versus 0.613), but the proposed historical-versus-contemporary mechanism was not supported.

| Primary Simpson effect | C3−N0 coefficient | Permutation P |
|---|---:|---:|
| Contemporary environment | −0.02365 | 0.291 |
| La Paz boundary | −0.00882 | 0.881 |
| Vizcaíno boundary | +0.00562 | 0.895 |
| Historical boundaries jointly | partial R² = 0.000062 | 0.975 |

Step 13G retained this conclusion across eight confirmatory sensitivity specifications. No specification supported the predicted stronger positive contemporary-environment effect, and none showed a significant joint historical-boundary difference in the predicted direction.

The exploratory Mulegé follow-up found no block-bootstrap-supported discrete local bump or general junction effect. The fixed 75-km effect and junction score were positive but uncertain, and the fitted Gulf-side curve indicated a slight local dip at the exact Mulegé latitude relative to adjacent flanks.

## Biological interpretation

> Ballooning-capable assemblages show lower overall compositional replacement than non-ballooning assemblages, but this difference is not explained by a simple contrast in which ballooning taxa track contemporary environmental conditions more strongly while non-ballooning taxa retain stronger signatures of historical biogeographic boundaries.

## Canonical order

1. 13A input audit
2. 13B frozen historical-boundary signals
3. 13C frozen contemporary-environment signal
4. 13D paired C3/N0 dissimilarities
5. 13E pairwise master table
6. 13F primary permutation test
7. 13G robustness/sensitivity
8. 13H1 input audit
9. 13H2 frozen Mulegé/ecotone predictors
10. 13H3 exploratory local inference

## Run

```bash
python scripts/supplementary/phase_13/run_phase13.py \
  --project-root /path/to/Baja_Ballooning_Pipeline
```

Large pairwise working tables are reconstructed locally and are not stored in this compact public package. Selected summaries, manifests, diagnostics, and figures are under `docs/phase_13/results/`.
