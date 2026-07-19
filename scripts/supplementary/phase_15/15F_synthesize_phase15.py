#!/usr/bin/env python3
"""Phase 15F — synthesize Bayesian estimation, prior robustness, prediction, and power."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

from phase15_common import phase15_root, read_csv, write_csv, write_json

SCRIPT_VERSION = "15F_v0.1.0_2026-07-18"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Synthesize Phase 15 evidence.")
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--output-dir", type=Path)
    return p.parse_args()


def one(path: Path):
    rows = read_csv(path)
    if not rows:
        raise RuntimeError(f"No rows in {path}")
    return rows[0]


def main() -> int:
    args = parse_args()
    root = args.project_root.expanduser().resolve()
    base = phase15_root(root)
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else base / "15F_evidence_synthesis"
    outdir.mkdir(parents=True, exist_ok=True)
    primary = one(base / "15B_primary_bayesian_model" / "15B_primary_result.csv")
    priors = read_csv(base / "15C_prior_sensitivity" / "15C_prior_sensitivity_results.csv")
    measurements = read_csv(base / "15C_prior_sensitivity" / "15C_measurement_error_sensitivity.csv")
    prediction = one(base / "15D_leave_one_cell_out_prediction" / "15D_predictive_comparison_summary.csv")
    power = one(base / "15E_design_power_simulation" / "15E_power_summary.csv")

    p_positive = float(primary["posterior_probability_beta_gt_0"])
    p_meaningful = float(primary["posterior_probability_beta_gt_0_10"])
    rope = float(primary["posterior_probability_abs_beta_le_0_05"])
    median = float(primary["stress_beta_median"])
    lo = float(primary["stress_beta_q025"])
    hi = float(primary["stress_beta_q975"])
    diagnostics = int(float(primary["diagnostics_pass"])) == 1
    prior_probs = [float(row["posterior_probability_beta_gt_0"]) for row in priors]
    prior_medians = [float(row["posterior_median"]) for row in priors]
    measurement_probs = [float(row["posterior_probability_beta_gt_0"]) for row in measurements]
    measurement_medians = [float(row["posterior_median"]) for row in measurements]
    sensitivity_min_probability = min(prior_probs + measurement_probs)
    prior_robust = (
        sensitivity_min_probability >= 0.90
        and all(value > 0 for value in prior_medians + measurement_medians)
    )
    predictive_gain = float(prediction["mean_crps_improvement_null_minus_stress"])
    predictive_p = float(prediction["exact_cell_sign_flip_p_two_sided"])

    if diagnostics and lo > 0 and prior_robust and predictive_gain > 0:
        status = "CONVERGENT_BAYESIAN_AND_PREDICTIVE_SUPPORT_FOR_H3"
        wording = (
            "The paired Bayesian model estimated a positive recent-stress effect whose 95% credible interval excluded zero, "
            "the direction was robust to prior choice, and the stress model improved held-out-cell prediction. This constitutes "
            "convergent exploratory support for H3 within the repeatedly sampled cells, but not peninsula-wide or causal proof."
        )
    elif diagnostics and p_positive >= 0.95 and prior_robust:
        status = "BAYESIAN_DIRECTIONAL_SUPPORT_BUT_NOT_DECISIVE"
        wording = (
            "The Bayesian analysis assigned high probability to a positive recent-stress effect and the direction was robust "
            "to prespecified priors, but the 95% credible interval or predictive comparison remained inconclusive. The result "
            "provides directional exploratory evidence rather than a definitive demonstration of H3."
        )
    elif median > 0:
        status = "POSITIVE_BUT_UNCERTAIN"
        wording = (
            "The Bayesian estimate remained positive, but posterior uncertainty, prior sensitivity, or out-of-cell prediction "
            "did not provide decisive corroboration. H3 remains plausible but unconfirmed."
        )
    else:
        status = "NO_DIRECTIONAL_SUPPORT_FOR_H3"
        wording = "The Phase 15 analysis did not support a positive recent-stress effect for ballooning assemblages."

    summary_rows = [
        {"item": "overall_status", "value": status},
        {"item": "posterior_median", "value": median},
        {"item": "posterior_95pct_interval", "value": f"[{lo}, {hi}]"},
        {"item": "posterior_probability_beta_gt_0", "value": p_positive},
        {"item": "posterior_probability_beta_gt_0_10", "value": p_meaningful},
        {"item": "posterior_probability_rope", "value": rope},
        {"item": "minimum_probability_positive_across_priors", "value": min(prior_probs)},
        {"item": "minimum_probability_positive_across_measurement_scales", "value": min(measurement_probs)},
        {"item": "minimum_probability_positive_across_all_sensitivities", "value": sensitivity_min_probability},
        {"item": "loco_mean_crps_improvement", "value": predictive_gain},
        {"item": "loco_exact_sign_flip_p", "value": predictive_p},
        {"item": "minimum_tested_effect_reaching_80pct_power", "value": power["minimum_grid_effect_reaching_80pct_power"]},
        {"item": "manuscript_wording", "value": wording},
        {"item": "scope_guardrail", "value": "Inference remains limited to 15 temporal comparisons in 12 repeatedly sampled cells."},
        {"item": "causal_guardrail", "value": "Observational temporal association cannot demonstrate stress-caused redistribution."},
    ]
    write_csv(outdir / "15F_conclusion_summary.csv", summary_rows)

    caption = (
        "Bayesian estimate of the difference in recent-stress-associated temporal turnover between ballooning-capable and "
        "non-ballooning assemblages. Positive values indicate a stronger turnover response among ballooning assemblages. "
        f"The posterior median was {median:.3f} with a 95% credible interval of [{lo:.3f}, {hi:.3f}], and "
        f"P(beta > 0) = {p_positive:.3f}. Results are restricted to repeatedly sampled cells and remain observational."
    )
    (outdir / "15F_figure_caption.txt").write_text(caption + "\n", encoding="utf-8")
    readme = f"""PHASE 15F — BAYESIAN H3 EVIDENCE SYNTHESIS
===========================================
Overall status: {status}

Primary posterior median: {median:.6f}
95% credible interval: [{lo:.6f}, {hi:.6f}]
P(beta > 0): {p_positive:.6f}
P(beta > 0.10): {p_meaningful:.6f}
P(|beta| <= 0.05): {rope:.6f}
Minimum P(beta > 0) across frozen priors: {min(prior_probs):.6f}
Minimum P(beta > 0) across observation-SD sensitivities: {min(measurement_probs):.6f}
LOCO CRPS improvement (null minus stress): {predictive_gain:.6f}
LOCO exact sign-flip p: {predictive_p:.6f}

Manuscript interpretation:
{wording}

Phase 15 is designed to improve estimation, propagate pair-level resampling
uncertainty, test prior dependence, and evaluate out-of-cell prediction. It does not
convert sensitivity models or Monte Carlo iterations into biological replicates.
"""
    (outdir / "15F_README.txt").write_text(readme, encoding="utf-8")
    write_json(outdir / "15F_run_status.json", {
        "phase": "15F", "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETED", "overall_status": status,
        "manuscript_wording": wording,
    })
    print(readme)
    print(f"OUTPUT_DIR={outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
