## [1.2.0] - 2026-07-24

### Added
- Step 11K publication QC with 5,000 paired equal-cell iterations, Baselga Jaccard replacement/nestedness, Sørensen-family sensitivity metrics, and input provenance.
- Incidence-based iNEXT Hill diversity at q = 0, 1, and 2 with shared-coverage standardization, 200 bootstrap replicates, recorded seed, curves, and R session information.
- Step 11N map-plus-inset main figure separating within-band Hill diversity from adjacent-boundary Baselga effects.
- Compact selected summaries and figures under `docs/step_11K/`; raw multi-million-row iterations remain locally generated.

### Scientific guardrails
- iNEXT Hill diversity is not described as a replacement–nestedness partition.
- Jaccard replacement plus Jaccard nestedness-resultant equals total Jaccard dissimilarity.
- Simpson replacement remains a Sørensen-family sensitivity metric.
- Boundary-level Hill beta q values remain supplementary and are not substituted for Baselga effects.

# Changelog

## [1.1.0] - 2026-07-18

### Added
- Phase 14 scripts 14A–14I, frozen temporal/stress configurations, Earth Engine extraction, sampling-continuity and trait-threshold sensitivities.
- Phase 15 scripts 15A–15G for Bayesian measurement-error estimation, prior sensitivity, held-out prediction, power simulation, evidence synthesis, and the final editable-text SVG figure.
- Compact selected empirical outputs, figures, provenance, and output manifests under `docs/phase_14/` and `docs/phase_15/`.

### Scientific result
- Phase 14 retained a mostly positive but uncertain H3 direction across 15 comparisons in 12 cells.
- Phase 15 estimated a positive standardized stress effect (`beta = 0.200426`, 95% CrI `[0.022920, 0.358910]`, `Pr(beta > 0 | data) = 0.989100`).
- Leave-one-cell-out prediction was not decisive (`exact p = 0.292969`); the official conclusion is positive but uncertain, not causal or peninsula-wide confirmation.
- The final synthesis figure retains the frozen spatial Simpson contrast (`C3 - N0 = -0.141745`, 95% interval `[-0.275743, -0.013406]`) and treats the latitude-boundary profile as descriptive.

### Reproducibility
- Phase 15G is locked to the original accepted 5,000-resample spatial input by checksum.
- Posterior draws, raw Earth Engine tables, credentials, caches, and synthetic validation data remain uncommitted.
- `SCRIPT_MANIFEST.tsv` is regenerated and verified.

<!-- PHASE13_CHANGELOG_START -->
## [1.0.6] - 2026-07-18

### Added
- Phase 13 scripts 13A–13H3 and the two frozen configuration files.
- `run_phase13.py` reproducibility wrapper.
- Compact selected outputs, manifests, diagnostics, and figures in `docs/phase_13/results/`.
- Documentation separating confirmatory 13F, robustness 13G, and exploratory 13H.

### Scientific result
- Confirmed lower overall Simpson replacement for C3 than N0.
- Found no support for stronger contemporary-environment tracking by C3 or stronger historical-boundary retention by N0.
- Confirmed the negative mechanism conclusion was robust across the pre-specified 13G analyses.
- Found no block-bootstrap-supported discrete Mulegé bump or general ecoregion-junction effect.

### Reproducibility
- Large pairwise working tables remain locally generated rather than committed.
- Local absolute paths are removed from copied public text outputs.
- `SCRIPT_MANIFEST.tsv` is regenerated.
<!-- PHASE13_CHANGELOG_END -->
