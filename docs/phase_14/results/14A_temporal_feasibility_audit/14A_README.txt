PHASE 14A — TEMPORAL FEASIBILITY AUDIT
======================================
Status: CONDITIONAL_LIMITED_H3_TEST

Input occurrence table: $PROJECT_ROOT/02_data_clean/05_final_qc_flags/05_biodiversity_final_records.tsv
Input trait table: $PROJECT_ROOT/ANALYSIS_READY_INPUTS/03_trait_tables/07_reviewed_genus_trait_lookup_normalized.csv
Primary analysis window: 2001-2025
Retained C3/N0 records in primary window: 5,027
Eligible adjacent cell-period comparisons: 15
Eligible cells: 12
Eligible latitude bands: 23-24N, 24-26N, 30-32N

Eligibility per trait group within each cell-period:
  records >= 5
  unique events >= 3
  genera >= 3

Recommendation:
Proceed only as an explicitly limited/exploratory H3 test, with equal-event resampling, dataset-continuity sensitivity analyses, and cautious geographic scope.

Interpretation:
This audit determines whether the existing occurrence data can support a repeated,
within-cell temporal comparison. It does not itself test stress tracking and does not
convert opportunistic GBIF records into standardized surveys.
