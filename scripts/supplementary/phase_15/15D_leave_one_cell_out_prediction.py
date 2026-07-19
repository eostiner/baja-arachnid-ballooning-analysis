#!/usr/bin/env python3
"""Phase 15D — leave-one-cell-out predictive comparison of stress versus null models."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from phase15_common import (
    build_design, empirical_crps, exact_sign_flip_pvalue, load_prepared,
    phase15_root, posterior_predictive, prior_sd_vector, run_chain,
    write_csv, write_json,
)

SCRIPT_VERSION = "15D_v0.1.0_2026-07-18"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run leave-one-cell-out Bayesian predictive validation.")
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--input-file", type=Path)
    p.add_argument("--output-dir", type=Path)
    p.add_argument("--chains", type=int, default=2)
    p.add_argument("--warmup", type=int, default=1200)
    p.add_argument("--draws", type=int, default=1800)
    p.add_argument("--seed", type=int, default=20260718)
    return p.parse_args()


def fit_fold(data, train_mask, include_stress, seed, args):
    X_all, names, _ = build_design(data, include_stress=include_stress, train_mask=train_mask)
    prior_sd = prior_sd_vector(names, 0.25)
    train_cells = [data.cells[i] for i in np.where(train_mask)[0]]
    unique_train = sorted(set(train_cells))
    mapping = {cell: j for j, cell in enumerate(unique_train)}
    train_cell_index = np.array([mapping[cell] for cell in train_cells], dtype=int)
    chains = [
        run_chain(data.y[train_mask], data.se[train_mask], X_all[train_mask], train_cell_index,
                  prior_sd, draws=args.draws, warmup=args.warmup, thin=1,
                  seed=seed + chain * 10007)
        for chain in range(args.chains)
    ]
    return X_all, chains


def new_cell_predictive(rng, chains, X_hold, se_hold, max_draws=2500):
    beta = np.concatenate([c["beta"] for c in chains], axis=0)
    sigma = np.concatenate([c["sigma"] for c in chains])
    tau = np.concatenate([c["tau"] for c in chains])
    if len(beta) > max_draws:
        idx = rng.choice(len(beta), max_draws, replace=False)
        beta, sigma, tau = beta[idx], sigma[idx], tau[idx]
    result = np.empty((len(beta), X_hold.shape[0]))
    for d in range(len(beta)):
        cell_effect = rng.normal(0.0, tau[d])
        result[d] = X_hold @ beta[d] + cell_effect + rng.standard_t(4.0, size=X_hold.shape[0]) * sigma[d] + rng.normal(0.0, se_hold)
    return result


def main() -> int:
    args = parse_args()
    root = args.project_root.expanduser().resolve()
    inp = args.input_file.expanduser().resolve() if args.input_file else phase15_root(root) / "15A_input_audit" / "15A_bayesian_model_input.csv"
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else phase15_root(root) / "15D_leave_one_cell_out_prediction"
    outdir.mkdir(parents=True, exist_ok=True)
    data = load_prepared(inp)
    rng = np.random.default_rng(args.seed + 555)
    pair_rows = []
    cell_rows = []
    for fold, cell in enumerate(data.unique_cells):
        print(f"LOCO fold {fold + 1}/{len(data.unique_cells)}: hold out {cell}", flush=True)
        hold = np.array([value == cell for value in data.cells], dtype=bool)
        train = ~hold
        X_stress, chains_stress = fit_fold(data, train, True, args.seed + fold * 200000, args)
        X_null, chains_null = fit_fold(data, train, False, args.seed + fold * 200000 + 100000, args)
        pred_stress = new_cell_predictive(rng, chains_stress, X_stress[hold], data.se[hold])
        pred_null = new_cell_predictive(rng, chains_null, X_null[hold], data.se[hold])
        cell_crps_stress = []
        cell_crps_null = []
        hold_indices = np.where(hold)[0]
        for local, global_i in enumerate(hold_indices):
            crps_s = empirical_crps(pred_stress[:, local], data.y[global_i])
            crps_n = empirical_crps(pred_null[:, local], data.y[global_i])
            cell_crps_stress.append(crps_s)
            cell_crps_null.append(crps_n)
            pair_rows.append({
                "held_out_cell": cell,
                "pair_id": data.pair_ids[global_i],
                "observed": data.y[global_i],
                "stress_model_predictive_median": float(np.median(pred_stress[:, local])),
                "null_model_predictive_median": float(np.median(pred_null[:, local])),
                "stress_model_crps": crps_s,
                "null_model_crps": crps_n,
                "crps_improvement_null_minus_stress": crps_n - crps_s,
            })
        cell_rows.append({
            "held_out_cell": cell,
            "n_pairs": int(np.sum(hold)),
            "mean_stress_model_crps": float(np.mean(cell_crps_stress)),
            "mean_null_model_crps": float(np.mean(cell_crps_null)),
            "mean_crps_improvement_null_minus_stress": float(np.mean(cell_crps_null) - np.mean(cell_crps_stress)),
        })
    write_csv(outdir / "15D_loco_pair_predictions.csv", pair_rows)
    write_csv(outdir / "15D_loco_cell_summary.csv", cell_rows)
    improvements = np.array([row["mean_crps_improvement_null_minus_stress"] for row in cell_rows], dtype=float)
    pvalue, patterns = exact_sign_flip_pvalue(improvements)
    mean_improvement = float(np.mean(improvements))
    cells_improved = int(np.sum(improvements > 0))
    if mean_improvement > 0 and pvalue < 0.05:
        classification = "STRESS_MODEL_PREDICTS_HELD_OUT_CELLS_BETTER"
    elif mean_improvement > 0:
        classification = "PREDICTIVE_IMPROVEMENT_POSITIVE_BUT_UNCERTAIN"
    else:
        classification = "NO_PREDICTIVE_IMPROVEMENT_OVER_NULL"
    summary = {
        "classification": classification,
        "n_held_out_cells": len(cell_rows),
        "cells_with_positive_crps_improvement": cells_improved,
        "mean_crps_improvement_null_minus_stress": mean_improvement,
        "exact_cell_sign_flip_p_two_sided": pvalue,
        "exact_sign_patterns": patterns,
    }
    write_csv(outdir / "15D_predictive_comparison_summary.csv", [summary])

    plot_outputs = []
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.6, 4.8))
        order = np.argsort(improvements)
        ax.barh(np.arange(len(order)), improvements[order])
        ax.axvline(0, linewidth=1.0)
        ax.set_yticks(np.arange(len(order)), [cell_rows[i]["held_out_cell"] for i in order], fontsize=7)
        ax.set_xlabel("CRPS improvement: null minus stress (positive favors stress model)")
        ax.set_title("Phase 15 leave-one-cell-out predictive validation")
        fig.tight_layout()
        for ext in ("png", "svg"):
            path = outdir / f"15D_loco_predictive_improvement.{ext}"
            fig.savefig(path, dpi=300 if ext == "png" else None)
            plot_outputs.append(str(path))
        plt.close(fig)
    except ImportError:
        pass

    readme = f"""PHASE 15D — LEAVE-ONE-CELL-OUT PREDICTIVE VALIDATION
========================================================
Classification: {classification}
Held-out cells: {len(cell_rows)}
Cells where stress model had lower CRPS: {cells_improved}
Mean CRPS improvement (null minus stress): {mean_improvement:.6f}
Exact cell sign-flip p: {pvalue:.6f} ({patterns} patterns)

This is a prediction check, not another test fitted to the same observations. Each
cell is withheld in turn, and the stress model is compared with a null model containing
the same sampling and time controls. Positive CRPS improvement means the stress model
predicted the held-out cell more accurately.
"""
    (outdir / "15D_README.txt").write_text(readme, encoding="utf-8")
    write_json(outdir / "15D_run_status.json", {
        "phase": "15D", "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETED", "summary": summary, "plot_outputs": plot_outputs,
        "chains_per_fold_model": args.chains, "warmup": args.warmup, "draws": args.draws,
    })
    print(readme)
    print(f"OUTPUT_DIR={outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
