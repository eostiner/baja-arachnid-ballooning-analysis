# Release 1.0.5 — final post-rerun publication package

This release incorporates every compatibility patch identified during the fresh end-to-end manuscript rerun. The R-enabled workflow completed successfully through the Spatial+ analysis and final mapping.

## Fresh-run validation

- 205 occupied cells were represented in the Step 12C table.
- 189 cells met the primary model criteria; 90 met the >=5-genus sensitivity criterion; 43 met the >=10-genus criterion.
- Step 12C passed 13/13 validation checks.
- Step 12K Spatial+ passed 13/13 validation checks.
- All primary and sensitivity models converged and all five leave-one-band-out folds ran.
- The centerpiece Step 10I map and the publication figure collection were regenerated.

## Final compatibility fixes

- Step 10B is independent of the internal `sf` geometry-column name and safely handles zero-overlap cells.
- Step 12B can audit already installed raster products without the original Downloads export directory.
- Step 12C distinguishes documented occupied cells from eligible primary-model cells when validating completeness and denominators.
- Steps 12D and 12E no longer compare the current data with hard-coded historical cell counts. They validate that the >=5 and >=10 sets are nested correctly and satisfy their thresholds.
- Step 12K1 reports the current modeled-cell count rather than assuming 195 cells.
- The trait normalizer accepts both CSV and the authoritative Excel workbook.

## Scientific scope

The primary trait definition remains **C3 = D1 + D2 + D3 versus fixed N0**, with D4 excluded. No metric formula, randomization seed, equal-cell design, or Spatial+ model specification was changed by these compatibility fixes.

The ordinary GLM remains a screening step. Passing the reduced-model diagnostic does not authorize an extrapolative habitat-suitability or prediction heat map; the retained final environmental analysis is Spatial+.
