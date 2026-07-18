PHASE 13C — FROZEN CONTEMPORARY ENVIRONMENT SIGNAL

Purpose:
Create an outcome-independent contemporary environmental dissimilarity signal
for the retained Baja 25-km cells.

Primary environmental distance:
1. Z-standardize frozen predictors across eligible cells.
2. Compute standardized Euclidean distance separately within thermal, moisture,
   wind, and vegetation domains.
3. Divide each domain distance by sqrt(number of variables in that domain).
4. Average the four domain distances with equal domain weight.

This prevents the moisture domain (which has more predictors) from dominating
the environmental signal simply because it contains more variables.

Predictor selection is not based on C3/N0 Phase 13 outcomes.
ECOSTRESS is not included because the retained Step 12L feasibility analysis
failed coverage thresholds. Existing ERA5-Land/MODIS variables provide broader
and already-audited peninsula-wide coverage.

Do not alter predictor membership after inspecting Phase 13 trait results.
Any alternative predictor sets must be labeled sensitivity analyses.
