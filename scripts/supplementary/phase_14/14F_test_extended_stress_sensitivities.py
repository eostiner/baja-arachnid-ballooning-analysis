#!/usr/bin/env python3
"""Phase 14F — prespecified extended stress model multiverse and influence diagnostics."""
from __future__ import annotations

import argparse
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from phase14_common import default_analysis_output_root, load_periods, read_delimited, write_csv, write_json
from phase14_model_utils import as_float, bh_fdr, classify, leave_one_cell_out, run_model

SCRIPT_VERSION = "14F_v0.3.0_2026-07-18"


def control_flags(control_set: str) -> tuple[bool, bool, bool]:
    if control_set == "full":
        return True, True, True
    if control_set == "effort_time":
        return False, True, True
    if control_set == "effort_only":
        return False, False, True
    if control_set == "unadjusted":
        return False, False, False
    raise ValueError(f"Unknown control_set: {control_set}")


def write_forest_plot(outdir: Path, rows: list[dict[str, Any]]) -> list[str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    plotted = [row for row in rows if row.get("status") == "OK" and row.get("control_set") == "full"]
    if not plotted:
        return []
    plotted = list(reversed(plotted))
    y = list(range(len(plotted)))
    estimates = [float(row["coefficient_per_predictor_sd"]) for row in plotted]
    lower = [float(row["cluster_bootstrap_q025"]) for row in plotted]
    upper = [float(row["cluster_bootstrap_q975"]) for row in plotted]
    labels = [str(row["model_id"]).replace("E_", "").replace("_", " ") for row in plotted]
    fig, ax = plt.subplots(figsize=(9.5, max(6.0, 0.42 * len(plotted) + 1.8)))
    ax.errorbar(estimates, y, xerr=[
        [estimate - lo for estimate, lo in zip(estimates, lower)],
        [hi - estimate for estimate, hi in zip(estimates, upper)],
    ], fmt="o", capsize=3)
    ax.axvline(0, linewidth=1)
    ax.set_yticks(y, labels)
    ax.set_xlabel("Change in C3 − N0 turnover per 1 SD worsening stress")
    ax.set_title("Phase 14F prespecified H3 stress sensitivities")
    fig.tight_layout()
    outputs = []
    for suffix in ("png", "svg"):
        path = outdir / f"14F_extended_stress_forest.{suffix}"
        fig.savefig(path, dpi=300 if suffix == "png" else None)
        outputs.append(str(path))
    plt.close(fig)
    return outputs


def write_loco_plot(outdir: Path, loco_rows: list[dict[str, Any]], primary_model_id: str) -> list[str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    selected = [row for row in loco_rows if row.get("model_id") == primary_model_id and row.get("status") == "OK"]
    if not selected:
        return []
    selected.sort(key=lambda row: float(row["coefficient"]))
    fig, ax = plt.subplots(figsize=(8.5, max(5.0, 0.35 * len(selected) + 1.5)))
    y = list(range(len(selected)))
    ax.scatter([float(row["coefficient"]) for row in selected], y)
    ax.axvline(0, linewidth=1)
    ax.set_yticks(y, [str(row["omitted_cell"]) for row in selected])
    ax.set_xlabel("Primary coefficient after omitting one cell")
    ax.set_ylabel("Omitted 25-km cell")
    ax.set_title("Phase 14F leave-one-cell-out influence")
    fig.tight_layout()
    outputs = []
    for suffix in ("png", "svg"):
        path = outdir / f"14F_primary_leave_one_cell_out.{suffix}"
        fig.savefig(path, dpi=300 if suffix == "png" else None)
        outputs.append(str(path))
    plt.close(fig)
    return outputs


def main() -> int:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Run prespecified extended H3 stress sensitivities.")
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--model-input", type=Path)
    parser.add_argument("--model-specs", type=Path, default=here / "configs" / "phase_14_extended_model_specs_frozen.csv")
    parser.add_argument("--period-config", type=Path, default=here / "configs" / "phase_14_temporal_windows_frozen.csv")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--bootstrap-iterations", type=int, default=5000)
    parser.add_argument("--permutation-iterations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260718)
    parser.add_argument("--min-pairs", type=int, default=10)
    parser.add_argument("--min-cells", type=int, default=8)
    args = parser.parse_args()

    root = args.project_root.expanduser().resolve()
    base = default_analysis_output_root(root)
    input_path = args.model_input.expanduser().resolve() if args.model_input else (
        base / "14E_extended_stress_join" / "14E2_H3_extended_model_input_all_pairs.csv"
    )
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else base / "14F_extended_stress_sensitivity"
    outdir.mkdir(parents=True, exist_ok=True)
    if not input_path.exists():
        raise FileNotFoundError(f"Extended model input not found: {input_path}")

    _, rows = read_delimited(input_path)
    _, specs = read_delimited(args.model_specs.expanduser().resolve())
    periods = {period["period_id"]: period for period in load_periods(args.period_config.expanduser().resolve())}

    results: list[dict[str, Any]] = []
    loco_rows: list[dict[str, Any]] = []
    for index, spec in enumerate(specs, start=1):
        model_id = str(spec["model_id"])
        predictor = str(spec["predictor"])
        response = str(spec["response"])
        control_set = str(spec.get("control_set", "full"))
        include_band, include_time, include_effort = control_flags(control_set)
        available = [
            row for row in rows
            if as_float(row.get(predictor)) is not None and as_float(row.get(response)) is not None
        ]
        cells = {str(row.get("grid_cell_id", "")) for row in available}
        base_result: dict[str, Any] = {
            **spec,
            "n_pairs": len(available),
            "n_cells": len(cells),
        }
        if len(available) < args.min_pairs or len(cells) < args.min_cells:
            results.append({**base_result, "status": "SKIPPED_INSUFFICIENT_COMPLETE_CASES"})
            continue
        try:
            model = run_model(
                rows, predictor, response, periods,
                args.bootstrap_iterations if spec.get("tier") == "primary" else max(1500, args.bootstrap_iterations // 2),
                args.permutation_iterations,
                args.seed + index * 1000,
                include_band=include_band,
                include_time=include_time,
                include_effort=include_effort,
            )
            model_cells = {item["cell"] for item in model["built"]["retained"]}
            results.append({
                **base_result,
                "status": "OK",
                "n_pairs": len(model["built"]["retained"]),
                "n_cells": len(model_cells),
                "coefficient_per_predictor_sd": model["beta"],
                "cluster_bootstrap_q025": model["ci_low"],
                "cluster_bootstrap_q975": model["ci_high"],
                "wild_cluster_p_two_sided": model["pvalue"],
                "wild_cluster_method": model["permutation_method"],
                "r_squared_descriptive": model["fit"]["r2"],
                "controls_retained": ";".join(model["built"]["names"][2:]),
                "controls_omitted": ";".join(model["built"]["omitted_controls"]),
            })
            for loco in leave_one_cell_out(
                rows, predictor, response, periods,
                include_band=include_band,
                include_time=include_time,
                include_effort=include_effort,
            ):
                loco_rows.append({"model_id": model_id, **loco})
        except Exception as exc:
            results.append({**base_result, "status": f"MODEL_FAILED: {exc}"})

    pvalues = [
        float(row["wild_cluster_p_two_sided"]) if row.get("status") == "OK" and row.get("tier") != "primary" else None
        for row in results
    ]
    qvalues = bh_fdr(pvalues)
    for row, qvalue in zip(results, qvalues):
        row["bh_fdr_q_across_nonprimary_models"] = qvalue if qvalue is not None else ""
        if row.get("status") == "OK":
            beta = float(row["coefficient_per_predictor_sd"])
            low = float(row["cluster_bootstrap_q025"])
            high = float(row["cluster_bootstrap_q975"])
            pvalue = float(row["wild_cluster_p_two_sided"])
            row["conclusion"] = classify(beta, low, high, pvalue, qvalue)
            matching_loco = [item for item in loco_rows if item["model_id"] == row["model_id"] and item.get("status") == "OK"]
            coefficients = [float(item["coefficient"]) for item in matching_loco]
            row["loco_valid_n"] = len(coefficients)
            row["loco_positive_fraction"] = sum(value > 0 for value in coefficients) / len(coefficients) if coefficients else ""
            row["loco_min_coefficient"] = min(coefficients) if coefficients else ""
            row["loco_max_coefficient"] = max(coefficients) if coefficients else ""

    write_csv(outdir / "14F_extended_model_results.csv", results)
    write_csv(outdir / "14F_leave_one_cell_out.csv", loco_rows)
    plot_outputs = write_forest_plot(outdir, results)
    plot_outputs += write_loco_plot(outdir, loco_rows, "E_PRIMARY_COMPOSITE_MEAN")

    primary = next((row for row in results if row.get("model_id") == "E_PRIMARY_COMPOSITE_MEAN"), None)
    if primary and primary.get("status") == "OK":
        primary_wording = (
            f"The expanded composite estimate was {float(primary['coefficient_per_predictor_sd']):.3f} "
            f"with a 95% cell-bootstrap interval from {float(primary['cluster_bootstrap_q025']):.3f} "
            f"to {float(primary['cluster_bootstrap_q975']):.3f}."
        )
    else:
        primary_wording = "The expanded primary model could not be fit under the frozen minimum-coverage rules."

    positive_full = [row for row in results if row.get("status") == "OK" and row.get("control_set") == "full" and float(row["coefficient_per_predictor_sd"]) > 0]
    significant_full = [row for row in positive_full if float(row["cluster_bootstrap_q025"]) > 0 and float(row["wild_cluster_p_two_sided"]) < 0.05]
    readme = f"""PHASE 14F — EXTENDED RECENT-STRESS SENSITIVITY ANALYSIS
========================================================
Status: COMPLETED_EXPLORATORY_MULTIVERSE
Models specified: {len(specs)}
Models fit successfully: {sum(row.get('status') == 'OK' for row in results)}
Full-control positive estimates: {len(positive_full)}
Full-control positive estimates with interval above zero and raw p<0.05: {len(significant_full)}

Primary expanded model:
{primary_wording}

Interpretation rule:
No individual secondary predictor is selected after seeing the result. Effect sizes,
cell-bootstrap intervals, exact wild-cluster p-values, false-discovery-rate values,
and leave-one-cell-out sign stability must be considered together. The biological
replication remains 15 comparisons in 12 cells, so every result remains exploratory.
"""
    (outdir / "14F_README.txt").write_text(readme, encoding="utf-8")
    status = {
        "phase": "14F",
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETED_EXPLORATORY_MULTIVERSE",
        "models_specified": len(specs),
        "models_fit": sum(row.get("status") == "OK" for row in results),
        "primary_model": primary or {},
        "plot_outputs": plot_outputs,
        "guardrail": "Expanded predictors improve measurement and robustness assessment but do not create new biological replicates.",
    }
    write_json(outdir / "14F_run_status.json", status)
    print(readme)
    print(f"OUTPUT_DIR={outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
