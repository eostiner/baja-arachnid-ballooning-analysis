#!/usr/bin/env python3
"""Guarded small-cluster regression utilities for Phase 14 exploratory models."""
from __future__ import annotations

import itertools
import math
from typing import Any

import numpy as np


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
    with np.errstate(divide="ignore", invalid="ignore"):
        t_values = beta / se
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


def build_design(
    rows: list[dict[str, str]],
    predictor: str,
    response: str,
    periods: dict[str, dict[str, Any]],
    include_band: bool = True,
    include_time: bool = True,
    include_effort: bool = True,
) -> dict[str, Any]:
    retained: list[dict[str, Any]] = []
    for row in rows:
        y = as_float(row.get(response))
        x = as_float(row.get(predictor))
        events = as_float(row.get("common_event_resample_n"))
        if y is None or x is None or events is None:
            continue
        p1 = periods.get(str(row.get("period_1", "")))
        p2 = periods.get(str(row.get("period_2", "")))
        if not p1 or not p2:
            continue
        midpoint = (float(p1["start_year"]) + float(p2["end_year"])) / 2.0
        retained.append({
            "row": row,
            "y": y,
            "x": x,
            "log_events": math.log1p(events),
            "midpoint": midpoint,
            "band": str(row.get("latitude_band", "unknown")),
            "cell": str(row.get("grid_cell_id", "")),
        })
    if not retained:
        raise RuntimeError(f"No complete rows for response={response}, predictor={predictor}")

    y = np.array([item["y"] for item in retained], dtype=float)
    x_raw = np.array([item["x"] for item in retained], dtype=float)
    x_z, x_mean, x_sd = zscore(x_raw)
    names = ["intercept", "stress_predictor_1SD"]
    design = np.column_stack([np.ones(len(retained)), x_z])
    omitted_controls: list[str] = []
    scaling: dict[str, float | str] = {"predictor_mean": x_mean, "predictor_sd": x_sd}

    def try_add_numeric(name: str, raw_values: np.ndarray) -> None:
        nonlocal design
        if len(raw_values) < 2:
            omitted_controls.append(f"{name}:too_few_values")
            return
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
        if candidate.shape[0] - candidate.shape[1] < 3:
            omitted_controls.append(f"{name}:insufficient_residual_df")
            return
        design = candidate
        names.append(name)
        scaling[f"{name}_mean"] = mean_value
        scaling[f"{name}_sd"] = sd

    if include_effort:
        try_add_numeric("log_common_events_1SD", np.array([item["log_events"] for item in retained], dtype=float))
    if include_time:
        try_add_numeric("transition_midpoint_1SD", np.array([item["midpoint"] for item in retained], dtype=float))

    bands = sorted({item["band"] for item in retained})
    reference_band = bands[0] if bands else ""
    if include_band:
        for band in bands[1:]:
            name = f"latitude_band_{band}_vs_{reference_band}"
            column = np.array([1.0 if item["band"] == band else 0.0 for item in retained])
            candidate = np.column_stack([design, column])
            if np.linalg.matrix_rank(candidate) <= np.linalg.matrix_rank(design):
                omitted_controls.append(f"{name}:collinear")
                continue
            if candidate.shape[0] - candidate.shape[1] < 3:
                omitted_controls.append(f"{name}:insufficient_residual_df")
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
    indices_by_cell = {
        cell: np.array([i for i, value in enumerate(cells) if value == cell], dtype=int)
        for cell in unique_cells
    }
    draws: list[float] = []
    for _ in range(iterations):
        sampled = rng.choice(unique_cells, size=len(unique_cells), replace=True)
        indices = np.concatenate([indices_by_cell[cell] for cell in sampled])
        try:
            fit = fit_ols(design[indices], y[indices])
        except (RuntimeError, np.linalg.LinAlgError):
            continue
        coefficient = float(fit["beta"][1])
        if math.isfinite(coefficient):
            draws.append(coefficient)
    return np.asarray(draws, dtype=float)


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
    return (extreme + 1.0) / (total + 1.0), total, method


def run_model(
    rows: list[dict[str, str]],
    predictor: str,
    response: str,
    periods: dict[str, dict[str, Any]],
    bootstrap_iterations: int,
    permutation_iterations: int,
    seed: int,
    include_band: bool = True,
    include_time: bool = True,
    include_effort: bool = True,
) -> dict[str, Any]:
    built = build_design(
        rows, predictor, response, periods,
        include_band=include_band,
        include_time=include_time,
        include_effort=include_effort,
    )
    fit = fit_ols(built["design"], built["y"])
    cells = [item["cell"] for item in built["retained"]]
    bootstrap = cluster_bootstrap(built["design"], built["y"], cells, bootstrap_iterations, seed)
    if len(bootstrap) < max(250, bootstrap_iterations // 3):
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
    }


def leave_one_cell_out(
    rows: list[dict[str, str]],
    predictor: str,
    response: str,
    periods: dict[str, dict[str, Any]],
    include_band: bool = True,
    include_time: bool = True,
    include_effort: bool = True,
) -> list[dict[str, Any]]:
    complete_cells = sorted({
        str(row.get("grid_cell_id", ""))
        for row in rows
        if as_float(row.get(predictor)) is not None and as_float(row.get(response)) is not None
    })
    output: list[dict[str, Any]] = []
    for omitted in complete_cells:
        subset = [row for row in rows if str(row.get("grid_cell_id", "")) != omitted]
        try:
            built = build_design(
                subset, predictor, response, periods,
                include_band=include_band,
                include_time=include_time,
                include_effort=include_effort,
            )
            fit = fit_ols(built["design"], built["y"])
            output.append({
                "omitted_cell": omitted,
                "n_pairs": len(built["retained"]),
                "n_cells": len({item["cell"] for item in built["retained"]}),
                "coefficient": float(fit["beta"][1]),
                "status": "OK",
            })
        except Exception as exc:
            output.append({
                "omitted_cell": omitted,
                "n_pairs": "",
                "n_cells": "",
                "coefficient": "",
                "status": f"FAILED: {exc}",
            })
    return output


def bh_fdr(pvalues: list[float | None]) -> list[float | None]:
    indexed = [(i, p) for i, p in enumerate(pvalues) if p is not None and math.isfinite(p)]
    if not indexed:
        return [None] * len(pvalues)
    indexed.sort(key=lambda item: item[1])
    m = len(indexed)
    adjusted = [None] * len(pvalues)
    running = 1.0
    for rank_rev, (index, pvalue) in enumerate(reversed(indexed), start=1):
        rank = m - rank_rev + 1
        qvalue = min(running, pvalue * m / rank)
        running = qvalue
        adjusted[index] = qvalue
    return adjusted


def classify(beta: float, ci_low: float, ci_high: float, pvalue: float, qvalue: float | None = None) -> str:
    if beta > 0 and ci_low > 0 and pvalue < 0.05 and (qvalue is None or qvalue < 0.10):
        return "EXPLORATORY_DIRECTIONAL_SUPPORT"
    if beta > 0:
        return "POSITIVE_BUT_UNCERTAIN"
    if beta < 0 and ci_high < 0 and pvalue < 0.05 and (qvalue is None or qvalue < 0.10):
        return "EXPLORATORY_OPPOSITE_DIRECTION"
    if beta < 0:
        return "NEGATIVE_BUT_UNCERTAIN"
    return "NO_DIRECTIONAL_SIGNAL"
