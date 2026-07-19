#!/usr/bin/env python3
"""
Phase 14D — limited exploratory test of H3 recent environmental-stress tracking.

Primary question:
Does worsening five-year environmental stress predict a larger temporal Simpson
replacement response for C3 than N0 assemblages?

Primary response:
  resampled_delta_simpson_C3_minus_N0_median
Primary predictor:
  delta_stress_composite_z_period2_minus_period1

Inference emphasizes effect size, cell-cluster bootstrap uncertainty, and an exact
(or Monte Carlo when necessary) wild-cluster sign-flip test. With the current Phase
14A status this analysis remains explicitly exploratory.
"""
from __future__ import annotations

import argparse
import itertools
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from phase14_common import default_analysis_output_root, load_periods, read_delimited, write_csv, write_json

SCRIPT_VERSION = "14D_v0.2.0_2026-07-18"
PRIMARY_RESPONSE = "resampled_delta_simpson_C3_minus_N0_median"
SECONDARY_RESPONSE = "resampled_delta_jaccard_C3_minus_N0_median"
PRIMARY_PREDICTOR = "delta_stress_composite_z_period2_minus_period1"
SECONDARY_PREDICTORS = (
    "abs_delta_stress_composite_z_period2_minus_period1",
    "delta_thermal_stress_z_period2_minus_period1",
    "delta_vpd_anomaly_z_period2_minus_period1",
    "delta_moisture_stress_z_period2_minus_period1",
    "delta_vegetation_stress_z_period2_minus_period1",
)


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Test limited exploratory temporal H3 recent-stress tracking.")
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--model-input", type=Path)
    parser.add_argument("--period-config", type=Path, default=here / "configs" / "phase_14_temporal_windows_frozen.csv")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--bootstrap-iterations", type=int, default=5000)
    parser.add_argument("--permutation-iterations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260718)
    parser.add_argument("--min-pairs", type=int, default=12)
    parser.add_argument("--min-cells", type=int, default=10)
    return parser.parse_args()


def as_float(value: Any) -> float | None:
    try:
        result = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def zscore(values: np.ndarray) -> tuple[np.ndarray, float, float]:
    mean_value = float(np.mean(values))
    sd_value = float(np.std(values, ddof=1))
    if not math.isfinite(sd_value) or sd_value <= 0:
        raise RuntimeError("Cannot standardize a constant or invalid predictor.")
    return (values - mean_value) / sd_value, mean_value, sd_value


def fit_ols(x: np.ndarray, y: np.ndarray) -> dict[str, Any]:
    beta, _, rank, _ = np.linalg.lstsq(x, y, rcond=None)
    if rank < x.shape[1]:
        raise RuntimeError("Model design matrix is rank deficient.")
    fitted = x @ beta
    residual = y - fitted
    n, p = x.shape
    rss = float(residual @ residual)
    df = n - p
    sigma2 = rss / df if df > 0 else float("nan")
    covariance = sigma2 * np.linalg.inv(x.T @ x) if df > 0 else np.full((p, p), np.nan)
    se = np.sqrt(np.diag(covariance))
    t_values = beta / se
    t_values[~np.isfinite(t_values)] = np.nan
    tss = float(((y - np.mean(y)) ** 2).sum())
    r2 = 1.0 - rss / tss if tss > 0 else float("nan")
    return {
        "beta": beta,
        "se": se,
        "t": t_values,
        "fitted": fitted,
        "residual": residual,
        "rss": rss,
        "df": df,
        "r2": r2,
    }


def build_design(rows: list[dict[str, str]], predictor: str, response: str, periods: dict[str, dict[str, Any]]):
    retained: list[dict[str, Any]] = []
    for row in rows:
        y = as_float(row.get(response))
        x = as_float(row.get(predictor))
        events = as_float(row.get("common_event_resample_n"))
        if y is None or x is None or events is None:
            continue
        p1 = periods.get(row.get("period_1", ""))
        p2 = periods.get(row.get("period_2", ""))
        if not p1 or not p2:
            continue
        midpoint = (float(p1["start_year"]) + float(p2["end_year"])) / 2.0
        retained.append(
            {
                "row": row,
                "y": y,
                "x": x,
                "log_events": math.log1p(events),
                "midpoint": midpoint,
                "band": str(row.get("latitude_band", "unknown")),
                "cell": str(row.get("grid_cell_id", "")),
            }
        )
    if not retained:
        raise RuntimeError(f"No complete rows for response={response}, predictor={predictor}")

    y = np.array([item["y"] for item in retained], dtype=float)
    x_raw = np.array([item["x"] for item in retained], dtype=float)
    x_z, x_mean, x_sd = zscore(x_raw)
    names = ["intercept", "stress_predictor_1SD"]
    design = np.column_stack([np.ones(len(retained)), x_z])
    omitted_controls: list[str] = []
    scaling: dict[str, float | str] = {
        "predictor_mean": x_mean,
        "predictor_sd": x_sd,
    }

    def try_add_numeric(name: str, raw_values: np.ndarray) -> None:
        nonlocal design
        sd = float(np.std(raw_values, ddof=1))
        if not math.isfinite(sd) or sd <= 0:
            omitted_controls.append(f"{name}:constant")
            return
        mean_value = float(np.mean(raw_values))
        column = (raw_values - mean_value) / sd
        candidate = np.column_stack([design, column])
        if np.linalg.matrix_rank(candidate) <= np.linalg.matrix_rank(design):
            omitted_controls.append(f"{name}:collinear")
            return
        design = candidate
        names.append(name)
        scaling[f"{name}_mean"] = mean_value
        scaling[f"{name}_sd"] = sd

    try_add_numeric(
        "log_common_events_1SD",
        np.array([item["log_events"] for item in retained], dtype=float),
    )
    try_add_numeric(
        "transition_midpoint_1SD",
        np.array([item["midpoint"] for item in retained], dtype=float),
    )

    bands = sorted({item["band"] for item in retained})
    reference_band = bands[0]
    for band in bands[1:]:
        name = f"latitude_band_{band}_vs_{reference_band}"
        column = np.array([1.0 if item["band"] == band else 0.0 for item in retained])
        candidate = np.column_stack([design, column])
        if np.linalg.matrix_rank(candidate) <= np.linalg.matrix_rank(design):
            omitted_controls.append(f"{name}:collinear")
            continue
        design = candidate
        names.append(name)

    return {
        "retained": retained,
        "y": y,
        "x_raw": x_raw,
        "x_z": x_z,
        "design": design,
        "names": names,
        "reference_band": reference_band,
        "scaling": scaling,
        "omitted_controls": omitted_controls,
    }


def cluster_bootstrap(design: np.ndarray, y: np.ndarray, cells: list[str], iterations: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    unique_cells = sorted(set(cells))
    indices_by_cell = {cell: np.array([i for i, value in enumerate(cells) if value == cell], dtype=int) for cell in unique_cells}
    draws: list[float] = []
    for _ in range(iterations):
        sampled = rng.choice(unique_cells, size=len(unique_cells), replace=True)
        indices = np.concatenate([indices_by_cell[cell] for cell in sampled])
        xb = design[indices]
        yb = y[indices]
        try:
            fit = fit_ols(xb, yb)
        except (RuntimeError, np.linalg.LinAlgError):
            continue
        coefficient = float(fit["beta"][1])
        if math.isfinite(coefficient):
            draws.append(coefficient)
    return np.array(draws, dtype=float)


def wild_cluster_pvalue(
    design: np.ndarray,
    y: np.ndarray,
    cells: list[str],
    observed: float,
    iterations: int,
    seed: int,
) -> tuple[float, int, str]:
    reduced = np.delete(design, 1, axis=1)
    reduced_fit = fit_ols(reduced, y)
    fitted0 = reduced_fit["fitted"]
    residual0 = reduced_fit["residual"]
    unique_cells = sorted(set(cells))
    cell_index = {cell: index for index, cell in enumerate(unique_cells)}
    row_cluster = np.array([cell_index[cell] for cell in cells], dtype=int)
    cluster_count = len(unique_cells)

    if cluster_count <= 15:
        sign_vectors = itertools.product((-1.0, 1.0), repeat=cluster_count)
        method = f"exact_wild_cluster_sign_flip_{2 ** cluster_count}_patterns"
    else:
        rng = np.random.default_rng(seed)
        sign_vectors = (rng.choice((-1.0, 1.0), size=cluster_count) for _ in range(iterations))
        method = f"monte_carlo_wild_cluster_sign_flip_{iterations}_draws"

    extreme = 0
    total = 0
    for signs in sign_vectors:
        sign_array = np.asarray(signs, dtype=float)
        y_star = fitted0 + residual0 * sign_array[row_cluster]
        try:
            beta_star = float(fit_ols(design, y_star)["beta"][1])
        except (RuntimeError, np.linalg.LinAlgError):
            continue
        total += 1
        if abs(beta_star) >= abs(observed) - 1e-12:
            extreme += 1
    pvalue = (extreme + 1.0) / (total + 1.0)
    return pvalue, total, method


def residualize(values: np.ndarray, controls: np.ndarray) -> np.ndarray:
    return values - controls @ np.linalg.lstsq(controls, values, rcond=None)[0]


def classify(beta: float, ci_low: float, ci_high: float, pvalue: float) -> str:
    if beta > 0 and ci_low > 0 and pvalue < 0.05:
        return "EXPLORATORY_DIRECTIONAL_SUPPORT_FOR_H3"
    if beta > 0:
        return "POSITIVE_BUT_UNCERTAIN"
    if beta < 0 and ci_high < 0 and pvalue < 0.05:
        return "EXPLORATORY_OPPOSITE_DIRECTION"
    if beta < 0:
        return "NO_DIRECTIONAL_SUPPORT_NEGATIVE_ESTIMATE"
    return "NO_DIRECTIONAL_SUPPORT_ZERO_ESTIMATE"


def run_model(rows, predictor, response, periods, bootstrap_iterations, permutation_iterations, seed):
    built = build_design(rows, predictor, response, periods)
    fit = fit_ols(built["design"], built["y"])
    cells = [item["cell"] for item in built["retained"]]
    bootstrap = cluster_bootstrap(built["design"], built["y"], cells, bootstrap_iterations, seed)
    if len(bootstrap) < max(500, bootstrap_iterations // 2):
        raise RuntimeError(f"Too few valid cluster-bootstrap fits: {len(bootstrap)}")
    ci_low, ci_high = np.quantile(bootstrap, [0.025, 0.975])
    observed = float(fit["beta"][1])
    pvalue, permutation_n, permutation_method = wild_cluster_pvalue(
        built["design"], built["y"], cells, observed, permutation_iterations, seed + 101
    )
    return {
        "built": built,
        "fit": fit,
        "bootstrap": bootstrap,
        "beta": observed,
        "ci_low": float(ci_low),
        "ci_high": float(ci_high),
        "pvalue": pvalue,
        "permutation_n": permutation_n,
        "permutation_method": permutation_method,
        "conclusion": classify(observed, float(ci_low), float(ci_high), pvalue),
    }


def write_primary_plot(outdir: Path, model: dict[str, Any]) -> list[str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    design = model["built"]["design"]
    y = model["built"]["y"]
    x = design[:, 1]
    controls = np.delete(design, 1, axis=1)
    x_resid = residualize(x, controls)
    y_resid = residualize(y, controls)
    order = np.argsort(x_resid)
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    ax.scatter(x_resid, y_resid, s=48)
    ax.plot(x_resid[order], model["beta"] * x_resid[order], linewidth=1.8)
    for xr, yr, item in zip(x_resid, y_resid, model["built"]["retained"]):
        ax.annotate(item["row"].get("pair_id", ""), (xr, yr), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.axhline(0, linewidth=0.8)
    ax.axvline(0, linewidth=0.8)
    ax.set_xlabel("Worsening composite stress (partial residual, SD units)")
    ax.set_ylabel("C3 − N0 temporal Simpson replacement (partial residual)")
    ax.set_title("Phase 14D exploratory H3 temporal stress test")
    note = f"β={model['beta']:.3f}; 95% cell-bootstrap CI [{model['ci_low']:.3f}, {model['ci_high']:.3f}]; p={model['pvalue']:.3f}"
    ax.text(0.02, 0.98, note, transform=ax.transAxes, va="top", fontsize=9)
    fig.tight_layout()
    outputs = []
    for suffix in ("png", "svg"):
        path = outdir / f"14D_primary_partial_effect.{suffix}"
        fig.savefig(path, dpi=300 if suffix == "png" else None)
        outputs.append(str(path))
    plt.close(fig)
    return outputs


def main() -> int:
    args = parse_args()
    if args.bootstrap_iterations < 500:
        raise ValueError("Use at least 500 cluster-bootstrap iterations.")
    root = args.project_root.expanduser().resolve()
    base = default_analysis_output_root(root)
    input_path = args.model_input.expanduser().resolve() if args.model_input else (
        base / "14C_temporal_stress_join" / "14C_H3_model_input_complete_cases.csv"
    )
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else base / "14D_recent_stress_tracking_test"
    outdir.mkdir(parents=True, exist_ok=True)
    if not input_path.exists():
        raise FileNotFoundError(f"Phase 14C complete-case model input not found: {input_path}")

    _, rows = read_delimited(input_path)
    periods_list = load_periods(args.period_config.expanduser().resolve())
    periods = {str(period["period_id"]): period for period in periods_list}
    available_primary = [
        row for row in rows
        if as_float(row.get(PRIMARY_RESPONSE)) is not None and as_float(row.get(PRIMARY_PREDICTOR)) is not None
    ]
    cells = {row.get("grid_cell_id", "") for row in available_primary}
    if len(available_primary) < args.min_pairs or len(cells) < args.min_cells:
        status = {
            "phase": "14D",
            "script_version": SCRIPT_VERSION,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "status": "NOT_RUN_INSUFFICIENT_COMPLETE_CASES",
            "complete_pairs": len(available_primary),
            "complete_cells": len(cells),
            "minimum_pairs": args.min_pairs,
            "minimum_cells": args.min_cells,
        }
        write_json(outdir / "14D_run_status.json", status)
        print(status)
        return 0

    primary = run_model(
        rows,
        PRIMARY_PREDICTOR,
        PRIMARY_RESPONSE,
        periods,
        args.bootstrap_iterations,
        args.permutation_iterations,
        args.seed,
    )

    coefficient_rows = []
    for name, beta, se, tvalue in zip(
        primary["built"]["names"], primary["fit"]["beta"], primary["fit"]["se"], primary["fit"]["t"]
    ):
        coefficient_rows.append(
            {
                "model": "primary_simpson_signed_composite_stress",
                "term": name,
                "estimate": float(beta),
                "classical_se_descriptive_only": float(se),
                "classical_t_descriptive_only": float(tvalue) if math.isfinite(float(tvalue)) else "",
            }
        )
    write_csv(outdir / "14D_primary_model_coefficients.csv", coefficient_rows)

    sensitivity_specs = [(PRIMARY_RESPONSE, predictor) for predictor in SECONDARY_PREDICTORS]
    sensitivity_specs.append((SECONDARY_RESPONSE, PRIMARY_PREDICTOR))
    sensitivity_rows = []
    for index, (response, predictor) in enumerate(sensitivity_specs, start=1):
        try:
            model = run_model(
                rows,
                predictor,
                response,
                periods,
                max(1000, args.bootstrap_iterations // 2),
                args.permutation_iterations,
                args.seed + index * 1000,
            )
            sensitivity_rows.append(
                {
                    "response": response,
                    "predictor": predictor,
                    "n_pairs": len(model["built"]["retained"]),
                    "n_cells": len({item["cell"] for item in model["built"]["retained"]}),
                    "coefficient_per_predictor_sd": model["beta"],
                    "cluster_bootstrap_q025": model["ci_low"],
                    "cluster_bootstrap_q975": model["ci_high"],
                    "wild_cluster_p_two_sided": model["pvalue"],
                    "conclusion": model["conclusion"],
                }
            )
        except Exception as exc:
            sensitivity_rows.append(
                {
                    "response": response,
                    "predictor": predictor,
                    "n_pairs": "",
                    "n_cells": "",
                    "coefficient_per_predictor_sd": "",
                    "cluster_bootstrap_q025": "",
                    "cluster_bootstrap_q975": "",
                    "wild_cluster_p_two_sided": "",
                    "conclusion": f"MODEL_FAILED: {exc}",
                }
            )
    write_csv(outdir / "14D_sensitivity_models.csv", sensitivity_rows)

    model_summary = {
        "model": "primary_simpson_signed_composite_stress",
        "response": PRIMARY_RESPONSE,
        "predictor": PRIMARY_PREDICTOR,
        "predictor_orientation": "positive means period 2 was more stressful than period 1",
        "coefficient_scale": "change in C3-minus-N0 Simpson replacement per 1 SD increase in stress worsening",
        "n_pairs": len(primary["built"]["retained"]),
        "n_cells": len({item["cell"] for item in primary["built"]["retained"]}),
        "latitude_band_reference": primary["built"]["reference_band"],
        "controls": "; ".join(primary["built"]["names"][2:]) or "none",
        "omitted_controls": "; ".join(primary["built"]["omitted_controls"]),
        "coefficient": primary["beta"],
        "cluster_bootstrap_q025": primary["ci_low"],
        "cluster_bootstrap_q975": primary["ci_high"],
        "wild_cluster_p_two_sided": primary["pvalue"],
        "wild_cluster_permutations": primary["permutation_n"],
        "wild_cluster_method": primary["permutation_method"],
        "r_squared_descriptive": primary["fit"]["r2"],
        "conclusion": primary["conclusion"],
    }
    write_csv(outdir / "14D_primary_model_summary.csv", [model_summary])
    plot_outputs = write_primary_plot(outdir, primary)

    if primary["conclusion"] == "EXPLORATORY_DIRECTIONAL_SUPPORT_FOR_H3":
        wording = (
            "Within the limited repeatedly sampled cells, worsening recent environmental stress was associated "
            "with a larger C3-minus-N0 temporal replacement response. This is exploratory observational support, "
            "not a peninsula-wide or causal demonstration."
        )
    elif primary["conclusion"] == "POSITIVE_BUT_UNCERTAIN":
        wording = (
            "The estimated relationship was in the H3-predicted direction, but uncertainty included no effect. "
            "The limited opportunistic data therefore do not demonstrate recent stress tracking."
        )
    elif primary["beta"] < 0:
        wording = (
            "The estimated relationship was opposite the H3 prediction. Under the limited exploratory design, "
            "there is no evidence that C3 temporal replacement tracked worsening recent stress more strongly than N0."
        )
    else:
        wording = "The limited exploratory analysis provided no directional evidence for H3."

    conclusion_rows = [
        {"item": "phase14A_scope", "value": "CONDITIONAL_LIMITED_H3_TEST"},
        {"item": "primary_result", "value": primary["conclusion"]},
        {"item": "manuscript_interpretation", "value": wording},
        {"item": "causal_guardrail", "value": "Association cannot establish that recent stress caused redistribution."},
        {"item": "geographic_guardrail", "value": "Inference is restricted to the eligible repeated cells and represented latitude bands."},
    ]
    write_csv(outdir / "14D_conclusion_summary.csv", conclusion_rows)
    readme = f"""PHASE 14D — RECENT ENVIRONMENTAL-STRESS TRACKING (H3)
=========================================================
Status: COMPLETED_LIMITED_EXPLORATORY_TEST
Pairs: {model_summary['n_pairs']}
Cells: {model_summary['n_cells']}

Primary coefficient per 1 SD worsening composite stress: {primary['beta']:.6f}
Cell-cluster bootstrap 95% interval: [{primary['ci_low']:.6f}, {primary['ci_high']:.6f}]
Wild-cluster two-sided p: {primary['pvalue']:.6f}
Conclusion: {primary['conclusion']}

Interpretation:
{wording}

This result must not be described as a peninsula-wide confirmatory test or as
causal evidence. The Phase 14A audit found only conditional temporal coverage.
"""
    (outdir / "14D_README.txt").write_text(readme, encoding="utf-8")
    status = {
        "phase": "14D",
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETED_LIMITED_EXPLORATORY_TEST",
        "primary": model_summary,
        "plot_outputs": plot_outputs,
        "bootstrap_iterations_requested": args.bootstrap_iterations,
        "bootstrap_iterations_valid": len(primary["bootstrap"]),
        "seed": args.seed,
    }
    write_json(outdir / "14D_run_status.json", status)
    print(readme)
    print(f"OUTPUT_DIR={outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
