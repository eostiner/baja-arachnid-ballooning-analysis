# Changelog

## 1.0.5 — 2026-07-16

Final post-rerun publication release.

- Incorporated the Step 10B fixes for arbitrary `sf` geometry-column names and occupied cells with no positive-area ecoregion overlap.
- Step 12B now reuses and audits installed environmental rasters when the original export folder is absent.
- Step 12C retains all 205 documented occupied cells but applies completeness and positive-denominator validation to the 189 primary-model candidates.
- Replaced obsolete fixed candidate counts in Steps 12D and 12E with dynamic checks of set nesting and >=5/>=10 denominator thresholds.
- Removed stale 195-cell wording from Steps 12E, 12K, and 12K1; map captions now report the current modeled-cell count.
- Step 10T and preflight now accept the authoritative Excel workbook directly, including the `Genus_Trait_Master_267` sheet and `final_tier_for_current_build` field.
- Updated Step 12E language so model screening proceeds to Spatial+ deconfounding rather than authorizing an extrapolative prediction heat map.
- The complete R-enabled workflow was rerun successfully through Step 12K1 on 2026-07-16.

## 1.0.4 — 2026-07-16

Step 11F optional-confidence compatibility fix.

- Treats a missing trait-confidence column as `UNSPECIFIED` rather than aborting.
- Preserves all genera in the explicit low-confidence sensitivity when no LOW designations exist.
- Adds a regression test using the authoritative trait schema without a confidence field.
- No richness formula, incidence matrix, randomization design, seed, or trait definition changed.

## 1.0.3 — 2026-07-16

Step 11B full-reference export fix.

- Passed the explicit `n0_mask`, `d4_mask`, and `classified_mask` arguments to the full-cell `band_metrics()` call.
- This fixes the post-resampling crash that occurred only after all 5,000 iterations had completed.
- Added a regression test that verifies every `band_metrics()` call supplies the complete C3/N0/D4 mask set.

## 1.0.2 — 2026-07-16

Trait-schema compatibility fix.

- Added Step 10T to combine `exclusive_tier` with `primary_C3_group`.
- Correctly resolves fixed N0 genera whose evidence-tier cell is blank but whose primary group is explicitly non-ballooning.
- Writes a canonical normalized trait table used by Steps 11, 11B, 11G, 12C, and supplementary 12J.
- Updated preflight, documentation, runner order, tests, and manifests.
- Added a synthetic regression test for the exact two-column authoritative schema.

## 1.0.1 — 2026-07-16

Publication rerun release.

- Enforced the primary C3 definition (D1+D2+D3) versus fixed N0 in Steps 11G and 12C; D4 is excluded rather than recoded as N0.
- Restored the Step 10I v7 support module and retained one canonical script per step.
- Corrected Figure 1 to plot and label mean expected genus richness.
- Updated the Step 10I centerpiece map legend to use the eight-cell ecoregion threshold and the term “a priori test latitude.”
- Made Step 10I sample-size and iteration annotations derive from the current output tables instead of fixed text.
- Added publication figure aliases, a locked visual reference, preflight checks, manifest validation, and GitHub Actions syntax checks.
- Kept GDM optional and supplementary; Spatial+ no longer depends on it.

## 1.0.0 — 2026-07-16

Initial canonical public release assembled from the retained analysis workflow.
