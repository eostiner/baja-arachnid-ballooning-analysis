#!/usr/bin/env python3
"""Phase 15E — design-matched Monte Carlo sensitivity and power analysis.

This step does not re-estimate the observed effect. It asks what effect sizes the
existing 15-pair/12-cell design could reliably detect using the exact cell sign-flip
procedure already used in Phase 14.
"""
from __future__ import annotations

import argparse
import itertools
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from phase15_common import build_design, load_prepared, phase15_root, read_csv, write_csv, write_json

SCRIPT_VERSION = "15E_v0.1.0_2026-07-18"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run Phase 15 design-matched power simulation.")
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--input-file", type=Path)
    p.add_argument("--posterior-summary", type=Path)
    p.add_argument("--output-dir", type=Path)
    p.add_argument("--simulations", type=int, default=2000)
    p.add_argument("--seed", type=int, default=20260718)
    p.add_argument("--effect-grid", default="0,0.05,0.10,0.15,0.20,0.25,0.30")
    return p.parse_args()


def exact_wild_pvalue(y: np.ndarray, X: np.ndarray, cell_index: np.ndarray, signs: np.ndarray) -> tuple[float, float]:
    stress_index = 1
    pinv = np.linalg.pinv(X)
    coefficient_weights = pinv[stress_index]
    beta = float(coefficient_weights @ y)
    X0 = np.delete(X, stress_index, axis=1)
    fitted0 = X0 @ np.linalg.lstsq(X0, y, rcond=None)[0]
    residual0 = y - fitted0
    cluster_scores = np.array([
        float(np.sum(coefficient_weights[cell_index == j] * residual0[cell_index == j]))
        for j in range(signs.shape[1])
    ])
    randomized = signs @ cluster_scores
    pvalue = float(np.mean(np.abs(randomized) >= abs(beta) - 1e-12))
    return beta, pvalue


def main() -> int:
    args = parse_args()
    if args.simulations < 200:
        raise ValueError("Use at least 200 simulations.")
    root = args.project_root.expanduser().resolve()
    inp = args.input_file.expanduser().resolve() if args.input_file else phase15_root(root) / "15A_input_audit" / "15A_bayesian_model_input.csv"
    posterior = args.posterior_summary.expanduser().resolve() if args.posterior_summary else phase15_root(root) / "15B_primary_bayesian_model" / "15B_posterior_parameter_summary.csv"
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else phase15_root(root) / "15E_design_power_simulation"
    outdir.mkdir(parents=True, exist_ok=True)
    data = load_prepared(inp)
    X, names, _ = build_design(data, include_stress=True)
    summaries = {row["parameter"]: row for row in read_csv(posterior)}
    beta_template = np.array([float(summaries[name]["median"]) for name in names])
    sigma = float(summaries["process_sd"]["median"])
    tau = float(summaries["cell_sd"]["median"])
    effects = [float(item) for item in args.effect_grid.split(",") if item.strip()]
    unique_cells = len(data.unique_cells)
    signs = np.asarray(list(itertools.product((-1.0, 1.0), repeat=unique_cells)), dtype=float)
    rng = np.random.default_rng(args.seed)
    rows = []
    for effect in effects:
        estimates = []
        significant_positive = 0
        positive = 0
        pvalues = []
        beta_true = beta_template.copy()
        beta_true[1] = effect
        for _ in range(args.simulations):
            u = rng.normal(0.0, tau, size=unique_cells)
            process = rng.standard_t(4.0, size=len(data.y)) * sigma
            observation = rng.normal(0.0, data.se)
            y = X @ beta_true + u[data.cell_index] + process + observation
            estimate, pvalue = exact_wild_pvalue(y, X, data.cell_index, signs)
            estimates.append(estimate)
            pvalues.append(pvalue)
            positive += int(estimate > 0)
            significant_positive += int(estimate > 0 and pvalue < 0.05)
        estimates = np.asarray(estimates)
        rows.append({
            "true_stress_effect": effect,
            "simulations": args.simulations,
            "probability_estimate_positive": positive / args.simulations,
            "exact_test_positive_power_p_lt_0_05": significant_positive / args.simulations,
            "median_estimated_effect": float(np.median(estimates)),
            "estimated_effect_q025": float(np.quantile(estimates, 0.025)),
            "estimated_effect_q975": float(np.quantile(estimates, 0.975)),
            "median_exact_pvalue": float(np.median(pvalues)),
        })
        print(f"Effect {effect:.2f}: power={rows[-1]['exact_test_positive_power_p_lt_0_05']:.3f}", flush=True)
    write_csv(outdir / "15E_design_power_curve.csv", rows)
    adequate = [row for row in rows if row["exact_test_positive_power_p_lt_0_05"] >= 0.80]
    min_80 = min((row["true_stress_effect"] for row in adequate), default=None)
    false_positive = next((row["exact_test_positive_power_p_lt_0_05"] for row in rows if abs(row["true_stress_effect"]) < 1e-12), None)
    summary = {
        "n_pairs": len(data.y), "n_cells": unique_cells,
        "process_sd_used": sigma, "cell_sd_used": tau,
        "simulations_per_effect": args.simulations,
        "exact_sign_patterns": len(signs),
        "estimated_false_positive_rate_at_zero": false_positive,
        "minimum_grid_effect_reaching_80pct_power": "not_reached" if min_80 is None else min_80,
    }
    write_csv(outdir / "15E_power_summary.csv", [summary])

    plot_outputs = []
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.2, 4.8))
        ax.plot([r["true_stress_effect"] for r in rows], [r["exact_test_positive_power_p_lt_0_05"] for r in rows], marker="o")
        ax.axhline(0.80, linewidth=1.0, linestyle="--")
        ax.set_ylim(0, 1)
        ax.set_xlabel("True C3 − N0 stress effect per 1 SD")
        ax.set_ylabel("Probability of positive exact-test detection")
        ax.set_title("Phase 15 design-matched power")
        fig.tight_layout()
        for ext in ("png", "svg"):
            path = outdir / f"15E_design_power_curve.{ext}"
            fig.savefig(path, dpi=300 if ext == "png" else None)
            plot_outputs.append(str(path))
        plt.close(fig)
    except ImportError:
        pass

    readme = f"""PHASE 15E — DESIGN-MATCHED MONTE CARLO POWER
================================================
Pairs: {len(data.y)}
Cells: {unique_cells}
Simulations per effect: {args.simulations}
Exact sign patterns per simulation: {len(signs)}
Estimated false-positive rate at beta=0: {false_positive}
Minimum tested effect reaching 80% power: {min_80 if min_80 is not None else 'not reached'}

This analysis does not strengthen the observed coefficient. It quantifies which true
effects the existing design could detect and distinguishes an inconclusive result from
evidence that the biological effect is near zero.
"""
    (outdir / "15E_README.txt").write_text(readme, encoding="utf-8")
    write_json(outdir / "15E_run_status.json", {
        "phase": "15E", "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETED", "summary": summary, "plot_outputs": plot_outputs,
        "effect_grid": effects, "seed": args.seed,
    })
    print(readme)
    print(f"OUTPUT_DIR={outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
