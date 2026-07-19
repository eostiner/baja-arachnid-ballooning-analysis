#!/usr/bin/env python3
"""Phase 15C — prior sensitivity, prior predictive audit, and posterior robustness."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from phase15_common import (
    build_design, load_prepared, phase15_root, prior_sd_vector, run_chain,
    summarize_chains, write_csv, write_json,
)

SCRIPT_VERSION = "15C_v0.1.0_2026-07-18"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run Phase 15 Bayesian prior sensitivity.")
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--input-file", type=Path)
    p.add_argument("--output-dir", type=Path)
    p.add_argument("--chains", type=int, default=4)
    p.add_argument("--warmup", type=int, default=2000)
    p.add_argument("--draws", type=int, default=3000)
    p.add_argument("--seed", type=int, default=20260718)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    root = args.project_root.expanduser().resolve()
    inp = args.input_file.expanduser().resolve() if args.input_file else phase15_root(root) / "15A_input_audit" / "15A_bayesian_model_input.csv"
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else phase15_root(root) / "15C_prior_sensitivity"
    outdir.mkdir(parents=True, exist_ok=True)
    data = load_prepared(inp)
    X, names, _ = build_design(data, include_stress=True)
    scales = [("skeptical", 0.15), ("regular_primary", 0.25), ("broad", 0.50)]
    rows = []
    posterior_by_scale = {}
    for index, (label, scale) in enumerate(scales):
        print(f"Running prior sensitivity: {label} (stress prior SD={scale}) ...", flush=True)
        prior_sd = prior_sd_vector(names, scale)
        chains = [
            run_chain(data.y, data.se, X, data.cell_index, prior_sd,
                      draws=args.draws, warmup=args.warmup, thin=1,
                      seed=args.seed + index * 100000 + chain * 10007)
            for chain in range(args.chains)
        ]
        summary, arrays = summarize_chains(chains, names)
        stress = arrays["stress_beta"].reshape(-1)
        posterior_by_scale[label] = stress
        stress_summary = next(row for row in summary if row["parameter"] == "stress_beta")
        rows.append({
            "prior_label": label,
            "stress_prior_sd": scale,
            "posterior_median": float(np.median(stress)),
            "posterior_q025": float(np.quantile(stress, 0.025)),
            "posterior_q975": float(np.quantile(stress, 0.975)),
            "posterior_probability_beta_gt_0": float(np.mean(stress > 0)),
            "posterior_probability_beta_gt_0_10": float(np.mean(stress > 0.10)),
            "posterior_probability_abs_beta_le_0_05": float(np.mean(np.abs(stress) <= 0.05)),
            "stress_beta_rhat": stress_summary["rhat"],
            "stress_beta_ess_bulk_approx": stress_summary["ess_bulk_approx"],
        })
    write_csv(outdir / "15C_prior_sensitivity_results.csv", rows)

    # The Phase 14 resampling interval is only an approximate observation-dispersion
    # scale. Test whether the H3 direction depends on halving or doubling that scale.
    measurement_rows = []
    measurement_posteriors = {}
    for index, multiplier in enumerate((0.5, 1.0, 2.0)):
        label = f"observation_sd_x{multiplier:g}"
        if multiplier == 1.0:
            stress = posterior_by_scale["regular_primary"]
        else:
            chains = [
                run_chain(data.y, data.se * multiplier, X, data.cell_index, prior_sd_vector(names, 0.25),
                          draws=args.draws, warmup=args.warmup, thin=1,
                          seed=args.seed + 500000 + index * 100000 + chain * 10007)
                for chain in range(args.chains)
            ]
            _, arrays = summarize_chains(chains, names)
            stress = arrays["stress_beta"].reshape(-1)
        measurement_posteriors[label] = stress
        measurement_rows.append({
            "measurement_label": label,
            "observation_sd_multiplier": multiplier,
            "posterior_median": float(np.median(stress)),
            "posterior_q025": float(np.quantile(stress, 0.025)),
            "posterior_q975": float(np.quantile(stress, 0.975)),
            "posterior_probability_beta_gt_0": float(np.mean(stress > 0)),
            "posterior_probability_beta_gt_0_10": float(np.mean(stress > 0.10)),
        })
    write_csv(outdir / "15C_measurement_error_sensitivity.csv", measurement_rows)

    # Prior predictive audit: use the frozen scales and the actual design. This asks whether
    # the priors regularly generate impossible C3-N0 turnover differences outside [-1, 1].
    rng = np.random.default_rng(args.seed + 777)
    prior_rows = []
    for label, scale in scales:
        prior_sd = prior_sd_vector(names, scale)
        simulations = 10000
        beta = rng.normal(0.0, prior_sd, size=(simulations, len(names)))
        sigma = np.abs(rng.normal(0.0, 0.25, size=simulations))
        tau = np.abs(rng.normal(0.0, 0.15, size=simulations))
        generated = []
        for d in range(simulations):
            u = rng.normal(0.0, tau[d], size=len(data.unique_cells))
            mu = X @ beta[d] + u[data.cell_index]
            ysim = mu + rng.standard_t(4.0, size=len(data.y)) * sigma[d] + rng.normal(0.0, data.se)
            generated.append(ysim)
        generated = np.asarray(generated)
        prior_rows.append({
            "prior_label": label,
            "stress_prior_sd": scale,
            "prior_predictive_fraction_outside_minus1_plus1": float(np.mean((generated < -1) | (generated > 1))),
            "prior_predictive_q005": float(np.quantile(generated, 0.005)),
            "prior_predictive_q995": float(np.quantile(generated, 0.995)),
        })
    write_csv(outdir / "15C_prior_predictive_audit.csv", prior_rows)

    min_prior_prob = min(float(row["posterior_probability_beta_gt_0"]) for row in rows)
    min_measurement_prob = min(float(row["posterior_probability_beta_gt_0"]) for row in measurement_rows)
    min_prob = min(min_prior_prob, min_measurement_prob)
    all_positive_medians = (
        all(float(row["posterior_median"]) > 0 for row in rows)
        and all(float(row["posterior_median"]) > 0 for row in measurement_rows)
    )
    status = "PRIOR_AND_MEASUREMENT_ROBUST_POSITIVE_DIRECTION" if all_positive_medians and min_prob >= 0.90 else (
        "POSITIVE_DIRECTION_BUT_SENSITIVITY_DEPENDENT" if all_positive_medians else "DIRECTION_NOT_ROBUST_TO_SENSITIVITY"
    )

    plot_outputs = []
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.4, 4.8))
        for label, _ in scales:
            values = posterior_by_scale[label]
            hist, edges = np.histogram(values, bins=60, density=True)
            centers = (edges[:-1] + edges[1:]) / 2
            ax.plot(centers, hist, label=label)
        ax.axvline(0, linewidth=1.0)
        ax.axvspan(-0.05, 0.05, alpha=0.10)
        ax.set_xlabel("Stress effect on C3 − N0 turnover")
        ax.set_ylabel("Posterior density")
        ax.set_title("Phase 15 prior sensitivity")
        ax.legend(frameon=False)
        fig.tight_layout()
        for ext in ("png", "svg"):
            path = outdir / f"15C_prior_sensitivity.{ext}"
            fig.savefig(path, dpi=300 if ext == "png" else None)
            plot_outputs.append(str(path))
        plt.close(fig)
    except ImportError:
        pass

    readme = f"""PHASE 15C — PRIOR SENSITIVITY AND PREDICTIVE AUDIT
=====================================================
Status: {status}
Minimum P(beta > 0) across all frozen sensitivities: {min_prob:.6f}
Minimum P(beta > 0) across priors: {min_prior_prob:.6f}
Minimum P(beta > 0) across observation-SD multipliers: {min_measurement_prob:.6f}
All posterior medians positive: {all_positive_medians}

Three zero-centered priors were specified before this Phase 15 fit: skeptical
(SD 0.15), regular (SD 0.25), and broad (SD 0.50). Because Phase 14 resampling
quantiles only approximate observation dispersion, the model is also repeated with
that scale halved and doubled. A sign change or sharp probability change across these
sensitivities must not be presented as robust Bayesian evidence.
"""
    (outdir / "15C_README.txt").write_text(readme, encoding="utf-8")
    write_json(outdir / "15C_run_status.json", {
        "phase": "15C", "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": status, "minimum_probability_positive": min_prob,
        "minimum_prior_probability_positive": min_prior_prob,
        "minimum_measurement_probability_positive": min_measurement_prob,
        "all_positive_medians": all_positive_medians,
        "plot_outputs": plot_outputs,
    })
    print(readme)
    print(f"OUTPUT_DIR={outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
