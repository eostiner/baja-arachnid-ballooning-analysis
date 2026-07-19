# Phase 15 — Bayesian hierarchical evidence and predictive validation for H3

## Purpose

Phase 14 found a mostly positive but uncertain relationship between worsening recent environmental stress and the difference in temporal turnover between ballooning-capable (C3) and non-ballooning (N0) assemblages. Phase 15 asks whether that signal becomes more informative when pair-level resampling uncertainty is propagated, repeated observations within cells are partially pooled, prior dependence is audited, and prediction is evaluated in held-out cells.

Phase 15 does **not** treat the Phase 14 sensitivity models as independent replications and does not claim that additional Monte Carlo iterations create new biological information.

## Frozen primary question

Does worsening contemporary stress predict a positive C3-minus-N0 temporal Simpson-replacement response?

The primary predictor remains the original Phase 14 composite:

`delta_stress_composite_z_period2_minus_period1`

The extended Phase 14 environmental metrics are not substituted after seeing their results.

## Model

For temporal comparison `i` in cell `j`:

- The Phase 14B median turnover contrast is an observed estimate of a latent pair-level contrast.
- The Phase 14B 2.5% and 97.5% resampling quantiles are converted to an approximate observation-dispersion scale.
- Latent pair-level responses follow a robust Student-t regression with a cell random intercept.
- Stress change, common-event depth, and transition midpoint are included.
- The stress predictor and controls are standardized.

The primary stress-effect prior is zero-centered Normal with SD 0.25. Frozen sensitivity priors use SD 0.15 and 0.50. Process and cell standard deviations receive half-Normal regularization. The sampler is a transparent NumPy Metropolis-within-Gibbs implementation; it requires no external Bayesian package.

## Steps

- **15A** audits and freezes the original Phase 14 input.
- **15B** fits the primary Bayesian hierarchical measurement-error model and performs posterior predictive checks.
- **15C** repeats the model under skeptical, regular, and broad priors, halves and doubles the approximate observation-dispersion scale, and performs prior predictive checks.
- **15D** withholds each cell in turn and compares the stress model with a matched null model using CRPS.
- **15E** runs design-matched Monte Carlo simulations with exact cell sign-flip inference to quantify detectable effect sizes.
- **15F** applies frozen synthesis rules and writes manuscript-safe wording.

## Interpretation hierarchy

1. Posterior computation must pass convergence and effective-sample diagnostics.
2. Report the posterior median, credible interval, `P(beta > 0)`, `P(beta > 0.10)`, and the probability within the practical-equivalence region `[-0.05, 0.05]`.
3. Inspect whether direction and probability are robust to all frozen priors.
4. Treat held-out-cell prediction as corroboration only when the stress model improves CRPS over the matched null model.
5. Use the power simulation to distinguish low precision from evidence of a negligibly small effect.

## Scope

All inference remains restricted to the repeated cell-period comparisons retained by Phase 14A. The analysis remains observational and cannot establish that environmental stress caused redistribution.
