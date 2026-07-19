#!/usr/bin/env python3
"""Phase 15B — fit the primary Bayesian paired hierarchical measurement-error model."""
from __future__ import annotations

import argparse
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from phase15_common import (
    build_design, load_prepared, phase15_root, posterior_predictive,
    prior_sd_vector, run_chain, summarize_chains, write_csv, write_json,
)

SCRIPT_VERSION = "15B_v0.1.0_2026-07-18"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Fit Phase 15 Bayesian H3 model.")
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--input-file", type=Path)
    p.add_argument("--output-dir", type=Path)
    p.add_argument("--chains", type=int, default=4)
    p.add_argument("--warmup", type=int, default=3000)
    p.add_argument("--draws", type=int, default=5000)
    p.add_argument("--thin", type=int, default=1)
    p.add_argument("--seed", type=int, default=20260718)
    p.add_argument("--beta-prior-scale", type=float, default=0.25)
    p.add_argument("--nu", type=float, default=4.0)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if args.chains < 2 or args.draws < 1000 or args.warmup < 1000:
        raise ValueError("Use at least 2 chains, 1000 warmup iterations, and 1000 retained draws per chain.")
    root = args.project_root.expanduser().resolve()
    inp = args.input_file.expanduser().resolve() if args.input_file else phase15_root(root) / "15A_input_audit" / "15A_bayesian_model_input.csv"
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else phase15_root(root) / "15B_primary_bayesian_model"
    outdir.mkdir(parents=True, exist_ok=True)
    data = load_prepared(inp)
    X, names, scaling = build_design(data, include_stress=True)
    prior_sd = prior_sd_vector(names, args.beta_prior_scale)
    chains = []
    for chain in range(args.chains):
        print(f"Running primary Bayesian chain {chain + 1}/{args.chains} ...", flush=True)
        chains.append(run_chain(
            data.y, data.se, X, data.cell_index, prior_sd,
            draws=args.draws, warmup=args.warmup, thin=args.thin,
            seed=args.seed + chain * 10007, nu=args.nu,
        ))
    summary_rows, arrays = summarize_chains(chains, names)
    write_csv(outdir / "15B_posterior_parameter_summary.csv", summary_rows)

    stress = arrays["stress_beta"].reshape(-1)
    p_positive = float(np.mean(stress > 0))
    p_meaningful = float(np.mean(stress > 0.10))
    p_rope = float(np.mean(np.abs(stress) <= 0.05))
    median = float(np.median(stress))
    q025, q975 = [float(x) for x in np.quantile(stress, [0.025, 0.975])]
    q055, q945 = [float(x) for x in np.quantile(stress, [0.055, 0.945])]

    rng = np.random.default_rng(args.seed + 999)
    yrep = posterior_predictive(rng, chains, X, data.cell_index, data.se, nu=args.nu)
    pred_lo, pred_hi = np.quantile(yrep, [0.05, 0.95], axis=0)
    ppc_rows = []
    for i, row in enumerate(data.rows):
        ppc_rows.append({
            "pair_id": row["pair_id"], "grid_cell_id": row["grid_cell_id"],
            "observed": data.y[i], "posterior_predictive_median": float(np.median(yrep[:, i])),
            "posterior_predictive_q05": float(pred_lo[i]), "posterior_predictive_q95": float(pred_hi[i]),
            "inside_90pct_interval": int(pred_lo[i] <= data.y[i] <= pred_hi[i]),
        })
    write_csv(outdir / "15B_posterior_predictive_by_pair.csv", ppc_rows)
    ppc_coverage = float(np.mean((pred_lo <= data.y) & (data.y <= pred_hi)))
    ppc_rmse = float(np.sqrt(np.mean((np.median(yrep, axis=0) - data.y) ** 2)))
    out_of_bounds = float(np.mean((yrep < -1) | (yrep > 1)))

    max_rhat = max(float(row["rhat"]) for row in summary_rows if math.isfinite(float(row["rhat"])))
    min_ess = min(float(row["ess_bulk_approx"]) for row in summary_rows if math.isfinite(float(row["ess_bulk_approx"])))
    stress_diag = next(row for row in summary_rows if row["parameter"] == "stress_beta")
    stress_rhat = float(stress_diag["rhat"])
    stress_ess = float(stress_diag["ess_bulk_approx"])
    diagnostics_pass = max_rhat <= 1.05 and stress_rhat <= 1.01 and stress_ess >= 400
    if q025 > 0 and diagnostics_pass:
        classification = "POSTERIOR_INTERVAL_EXCLUDES_ZERO_POSITIVE"
    elif p_positive >= 0.95 and diagnostics_pass:
        classification = "BAYESIAN_DIRECTIONAL_EVIDENCE_POSITIVE"
    elif median > 0:
        classification = "POSITIVE_BUT_UNCERTAIN"
    else:
        classification = "NO_POSITIVE_DIRECTIONAL_EVIDENCE"

    result = {
        "model": "primary_original_composite_bayesian_measurement_error",
        "n_pairs": len(data.y), "n_cells": len(data.unique_cells),
        "stress_beta_median": median,
        "stress_beta_q025": q025, "stress_beta_q975": q975,
        "stress_beta_q055": q055, "stress_beta_q945": q945,
        "posterior_probability_beta_gt_0": p_positive,
        "posterior_probability_beta_gt_0_10": p_meaningful,
        "posterior_probability_abs_beta_le_0_05": p_rope,
        "prior_scale_stress_beta": args.beta_prior_scale,
        "student_t_df": args.nu,
        "max_rhat": max_rhat, "minimum_ess_bulk_approx": min_ess,
        "stress_beta_rhat": stress_rhat, "stress_beta_ess_bulk_approx": stress_ess,
        "diagnostics_pass": int(diagnostics_pass),
        "posterior_predictive_90pct_coverage": ppc_coverage,
        "posterior_predictive_rmse": ppc_rmse,
        "posterior_predictive_fraction_outside_response_bounds": out_of_bounds,
        "classification": classification,
    }
    write_csv(outdir / "15B_primary_result.csv", [result])
    draw_rows = []
    max_export = min(args.draws, 5000)
    for chain_id, chain in enumerate(chains, start=1):
        for d in range(min(max_export, len(chain["beta"]))):
            row = {"chain": chain_id, "draw": d + 1, "process_sd": chain["sigma"][d], "cell_sd": chain["tau"][d]}
            for k, name in enumerate(names):
                row[name] = chain["beta"][d, k]
            draw_rows.append(row)
    write_csv(outdir / "15B_posterior_draws.csv", draw_rows)

    plot_outputs = []
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.2, 4.8))
        ax.hist(stress, bins=50, density=True, alpha=0.75)
        ax.axvline(0, linewidth=1.0)
        ax.axvline(0.10, linewidth=1.0, linestyle="--")
        ax.axvspan(-0.05, 0.05, alpha=0.12)
        ax.set_xlabel("Stress effect on C3 − N0 turnover per 1 SD worsening stress")
        ax.set_ylabel("Posterior density")
        ax.set_title("Phase 15 primary Bayesian H3 effect")
        ax.text(0.02, 0.97, f"median={median:.3f}\n95% CrI [{q025:.3f}, {q975:.3f}]\nP(β>0)={p_positive:.3f}", transform=ax.transAxes, va="top")
        fig.tight_layout()
        for ext in ("png", "svg"):
            path = outdir / f"15B_primary_stress_posterior.{ext}"
            fig.savefig(path, dpi=300 if ext == "png" else None)
            plot_outputs.append(str(path))
        plt.close(fig)
    except ImportError:
        pass

    readme = f"""PHASE 15B — PRIMARY BAYESIAN H3 MODEL
======================================
Classification: {classification}
Pairs: {len(data.y)}
Cells: {len(data.unique_cells)}

Stress effect median: {median:.6f}
95% credible interval: [{q025:.6f}, {q975:.6f}]
89% credible interval: [{q055:.6f}, {q945:.6f}]
Posterior probability beta > 0: {p_positive:.6f}
Posterior probability beta > 0.10: {p_meaningful:.6f}
Posterior probability |beta| <= 0.05: {p_rope:.6f}

Maximum split R-hat: {max_rhat:.4f}
Minimum approximate bulk ESS: {min_ess:.1f}
Stress-effect split R-hat: {stress_rhat:.4f}
Stress-effect approximate bulk ESS: {stress_ess:.1f}
Posterior-predictive 90% interval coverage: {ppc_coverage:.3f}
Posterior-predictive fraction outside [-1,1]: {out_of_bounds:.4f}

The posterior probability is a graded measure of evidence, not a conversion of the
15 comparisons into additional biological replicates. The model propagates Phase 14
resampling dispersion and partially pools repeated observations within cells.
"""
    (outdir / "15B_README.txt").write_text(readme, encoding="utf-8")
    write_json(outdir / "15B_run_status.json", {
        "phase": "15B", "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETED", "result": result, "scaling": scaling,
        "parameter_names": names, "plot_outputs": plot_outputs,
        "chains": args.chains, "warmup": args.warmup, "draws_per_chain": args.draws,
        "seed": args.seed,
    })
    print(readme)
    print(f"OUTPUT_DIR={outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
