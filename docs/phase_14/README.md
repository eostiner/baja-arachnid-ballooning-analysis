# Phase 14 — Recent environmental-stress tracking (H3)

## Scientific question

Phase 14 tests whether ballooning-capable C3 assemblages changed through time in
closer correspondence with recent environmental stress than fixed N0 assemblages.
The paired response is temporal turnover within the same retained 25-km cell:

`C3 temporal turnover − N0 temporal turnover`

A positive stress coefficient is in the H3-predicted direction.

## Current empirical scope

The real temporal audit found only 15 eligible adjacent-period comparisons in 12
cells, representing three latitude bands. Consequently every Phase 14 result is a
limited exploratory observational test. It is not a peninsula-wide confirmatory
test and cannot establish that environmental stress caused redistribution.

## Canonical steps

### 14A — Temporal feasibility audit

Recovers dates from final-QC occurrences, joins authoritative traits, assigns the
retained 25-km grid, applies frozen minimum records/events/genera thresholds, and
returns PASS, CONDITIONAL, or FAIL.

### 14B — Paired temporal turnover

Builds C3 and N0 assemblages in frozen five-year periods and calculates Jaccard and
Simpson turnover after equal-event resampling. The primary response is
`resampled_delta_simpson_C3_minus_N0_median`.

### 14C–14D — Original real-stress test

Uses ERA5-Land air temperature, VPD, precipitation and soil water plus MODIS EVI.
The completed real-data result was positive but highly uncertain:

- coefficient: +0.064450 per 1 SD worsening composite stress
- 95% cell-cluster bootstrap interval: -1.071152 to +0.307895
- wild-cluster p: 0.616793

This does not demonstrate H3.

### 14E0 — ECOSTRESS coverage audit

Queries `NASA/ECOSTRESS/L2T_LSTE/V2` for intersections with the eligible Baja cells.
ECOSTRESS is excluded unless spatial and temporal coverage passes. Scene intersection
alone is not sufficient; repeated high-quality observations in the eligible periods
would also be required.

### 14E1 — Extended annual stress extraction

Adds independent products and event-focused metrics:

- **ERA5-Land:** annual and upper-tail air temperature, hot days, mean and upper-tail
  VPD, high-VPD frequency, and root-zone soil water.
- **MODIS Terra + Aqua LST:** clear-sky daytime and nighttime land-surface temperature,
  upper-tail temperature, maximum temperature, hot-observation frequency, and counts.
- **CHIRPS:** annual rainfall, dry-day frequencies, driest-month rainfall, heavy-rain
  days, and maximum one-day rainfall.
- **TerraClimate:** climatic water deficit, PDSI, soil moisture, and AET/PET ratio.
- **MODIS EVI:** mean EVI and lower-tail EVI.

Each raw metric is standardized within cell against 2001–2020. Signs are oriented so
positive values always mean greater stress. Domain composites prevent individual
variables from silently changing the primary hypothesis.

### 14E2 — Prespecified temporal stress summaries

For every annual stress predictor and five-year period, calculates:

1. mean stress;
2. worst annual stress;
3. positive stress burden;
4. fraction of years at least one within-cell SD above baseline.

The primary expanded predictor is the period-2 minus period-1 change in mean extended
composite stress. Alternative summaries are frozen sensitivities, not candidates for
post-hoc selection.

### 14F — Extended environmental sensitivities

Runs a frozen set of models covering the composite and thermal, surface-heat, VPD,
rainfall, soil-water, water-balance, and vegetation domains. It also includes
independent single-product checks, Jaccard response sensitivity, and reduced-control
models.

Every model reports:

- effect per 1 SD worsening stress;
- cell-cluster bootstrap interval;
- exact wild-cluster sign-flip p-value when there are 15 or fewer cells;
- Benjamini–Hochberg false-discovery-rate value across nonprimary models;
- leave-one-cell-out coefficient range and sign stability.

No secondary predictor is selected after the results are seen.

### 14G — Sampling-continuity sensitivity

Reconstructs temporal turnover after restricting records to:

- datasets represented in both periods and both trait groups;
- calendar quarters represented in both periods and both trait groups;
- both dataset and quarter continuity simultaneously.

These checks address changes in data contributors and collection season. They normally
reduce sample size and therefore cannot prove the absence of sampling bias.

### 14H — Trait-threshold temporal sensitivity

Repeats temporal eligibility, turnover, and guarded H3 models for:

- C1: D1 only;
- C2: D1 + D2;
- C3: D1 + D2 + D3, the primary definition.

C4 is not reconstructed from the primary temporal index because D4 was intentionally
excluded upstream.

### 14I — Evidence synthesis

Combines the original primary model, expanded primary model, major environmental
domains, sampling-continuity checks, and trait-threshold checks into one table and
figure. The synthesis determines manuscript language so that the strongest-looking
individual result is not cherry-picked.

## Interpretation hierarchy

1. **Repeated exploratory directional support:** both primary models are positive and
   their intervals exclude zero, with no clear contradiction from robustness checks.
2. **Positive but uncertain across primary models:** both primary estimates are
   positive but at least one interval includes zero.
3. **Mixed:** estimates vary materially across domains or robustness restrictions.
4. **No consistent directional support:** the expanded analyses do not retain the
   H3-predicted direction.

Even the strongest category remains exploratory because biological replication is
limited and opportunistic.

## Required manuscript guardrails

- Do not call any Phase 14 result peninsula-wide.
- Do not describe a positive coefficient alone as support.
- Do not select a predictor because it has the smallest p-value.
- Do not interpret ECOSTRESS without passing the coverage and quality gates.
- Do not infer causation from temporal association.
- State the number of temporal comparisons and cells whenever H3 is discussed.
