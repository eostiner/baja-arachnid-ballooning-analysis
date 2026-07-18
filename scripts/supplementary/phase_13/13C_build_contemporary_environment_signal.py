#!/usr/bin/env python3
"""
13C_build_contemporary_environment_signal.py

Purpose
-------
Build a frozen, outcome-independent contemporary-environment signal from the
retained Phase 12C 205-cell environmental table.

The predictor set is frozen by ecological domain BEFORE Phase 13 compares
ballooning (C3) and non-ballooning (N0) assemblages.

This script deliberately does NOT:
- inspect C3/N0 coefficients or P-values,
- select predictors based on trait outcomes,
- use Step 12D/12E model ranking to choose variables,
- use ECOSTRESS because the retained Step 12L audit failed coverage thresholds,
- test historical boundaries against assemblage outcomes.

Primary contemporary signal
---------------------------
Four balanced domains:
  thermal    : tmean_c, tseason_monthly_sd_c
  moisture   : precip_annual_mean_mm, precip_interannual_cv_pct,
               vpd_mean_kpa, soil_water_mean_frac
  wind       : wind_speed_mean_ms, wind_monthly_sd_ms
  vegetation : evi_mean, evi_interannual_cv_pct

Each predictor is z-standardized across eligible cells.
For every cell pair, Euclidean distance is calculated within each domain and
divided by sqrt(number of predictors in that domain), so domains with more
variables do not receive extra weight. The primary environmental distance is
the unweighted mean of the four domain distances.

Outputs
-------
13C_predictor_config_frozen.csv
13C_cell_environment_standardized.csv
13C_pairwise_environment_distance.csv
13C_environment_domain_summary.csv
13C_predictor_correlation.csv
13C_run_manifest.txt
13C_README.txt
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args():
    p = argparse.ArgumentParser(
        description="Build frozen Phase 13 contemporary-environment signal."
    )
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument(
        "--config",
        type=Path,
        default=Path.home()
        / "Downloads/PHASE_13_STARTER/configs/phase_13_contemporary_predictors_frozen.csv",
    )
    p.add_argument("--environment-table", type=Path, default=None)
    p.add_argument("--output-dir", type=Path, default=None)
    return p.parse_args()


def validate_config(cfg: pd.DataFrame) -> None:
    required = {"domain", "predictor", "direction", "role", "notes"}
    missing = required - set(cfg.columns)
    if missing:
        raise ValueError(f"Predictor config missing columns: {sorted(missing)}")

    if cfg["predictor"].duplicated().any():
        dups = cfg.loc[cfg["predictor"].duplicated(), "predictor"].tolist()
        raise ValueError(f"Predictors must be unique: {dups}")

    if set(cfg["role"]) != {"primary"}:
        raise ValueError("13C v1 config must contain only frozen primary predictors.")

    expected_domains = {"thermal", "moisture", "wind", "vegetation"}
    if set(cfg["domain"]) != expected_domains:
        raise ValueError(
            f"Domains must be exactly {sorted(expected_domains)}; "
            f"found {sorted(set(cfg['domain']))}"
        )


def zscore(series: pd.Series):
    mean = float(series.mean())
    sd = float(series.std(ddof=1))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError(f"Cannot standardize {series.name}: SD={sd}")
    return (series - mean) / sd, mean, sd


def main():
    args = parse_args()
    root = args.project_root.expanduser().resolve()

    env_path = args.environment_table or (
        root
        / "04_analysis_USE _THIS"
        / "12C_cell_environment_model_table"
        / "12C_cell_environment_model_table.csv"
    )
    config_path = args.config.expanduser().resolve()

    outdir = args.output_dir or (
        root
        / "04_analysis_USE _THIS"
        / "13_historical_vs_contemporary"
        / "13C_contemporary_environment_signal"
    )
    outdir = outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print("PHASE 13C — CONTEMPORARY ENVIRONMENT SIGNAL")
    print("=" * 76)
    print(f"PROJECT ROOT : {root}")
    print(f"ENV TABLE    : {env_path}")
    print(f"CONFIG       : {config_path}")
    print(f"OUTPUT DIR   : {outdir}")

    if not env_path.exists():
        raise FileNotFoundError(f"Environment table not found: {env_path}")
    if not config_path.exists():
        raise FileNotFoundError(f"Frozen predictor config not found: {config_path}")

    env = pd.read_csv(env_path)
    cfg = pd.read_csv(config_path)
    validate_config(cfg)

    if len(env) != 205:
        raise ValueError(
            f"Expected retained 205-cell Phase 12C table; found {len(env)}."
        )
    if env["grid_cell_id"].duplicated().any():
        raise ValueError("Environment table has duplicated grid_cell_id values.")

    predictors = cfg["predictor"].tolist()
    missing_cols = [c for c in predictors if c not in env.columns]
    if missing_cols:
        raise ValueError(f"Environment table missing predictors: {missing_cols}")

    # Freeze eligibility without referring to C3/N0 outcomes.
    # Prefer the retained Phase 12 recommendation flag; otherwise require complete predictors.
    if "recommended_primary_model_cell" in env.columns:
        eligible_flag = env["recommended_primary_model_cell"].fillna(False).astype(bool)
        eligibility_rule = "recommended_primary_model_cell == True"
    else:
        eligible_flag = env[predictors].notna().all(axis=1)
        eligibility_rule = "complete selected predictors"

    eligible = env.loc[eligible_flag].copy()
    eligible = eligible.loc[eligible[predictors].notna().all(axis=1)].copy()

    if len(eligible) < 150:
        raise ValueError(
            f"Only {len(eligible)} eligible cells remain; expected broad retained coverage."
        )

    print("\nFROZEN CONTEMPORARY PREDICTOR DESIGN")
    print(cfg[["domain", "predictor", "role", "notes"]].to_string(index=False))
    print(f"\nELIGIBILITY RULE: {eligibility_rule}")
    print(f"ELIGIBLE CELLS  : {len(eligible)} / {len(env)}")

    # Standardize predictors using eligible cells only.
    standardized = eligible[
        [
            "grid_cell_id", "grid_cell_order", "centroid_latitude",
            "centroid_longitude", "latitude_band"
        ]
    ].copy()

    stats_rows = []
    for pred in predictors:
        z, mean, sd = zscore(eligible[pred].astype(float))
        standardized[f"z_{pred}"] = z.values
        stats_rows.append({
            "predictor": pred,
            "domain": cfg.loc[cfg["predictor"] == pred, "domain"].iloc[0],
            "n_cells": len(eligible),
            "mean_raw": mean,
            "sd_raw": sd,
            "min_raw": float(eligible[pred].min()),
            "max_raw": float(eligible[pred].max()),
        })

    stats = pd.DataFrame(stats_rows)

    # Correlation is diagnostic only; no predictor is dropped based on it.
    corr = eligible[predictors].corr(method="pearson")
    corr.index.name = "predictor"

    # Pairwise, domain-balanced environmental distance.
    domains = cfg["domain"].drop_duplicates().tolist()
    pred_by_domain = {
        d: cfg.loc[cfg["domain"] == d, "predictor"].tolist() for d in domains
    }

    pair_rows = []
    std_index = standardized.set_index("grid_cell_id")

    meta = eligible.set_index("grid_cell_id")[
        ["centroid_latitude", "centroid_longitude", "latitude_band"]
    ]

    ids = standardized["grid_cell_id"].tolist()
    for a, b in combinations(ids, 2):
        row = {
            "cell_i": a,
            "cell_j": b,
            "lat_i": float(meta.loc[a, "centroid_latitude"]),
            "lat_j": float(meta.loc[b, "centroid_latitude"]),
            "lon_i": float(meta.loc[a, "centroid_longitude"]),
            "lon_j": float(meta.loc[b, "centroid_longitude"]),
            "band_i": meta.loc[a, "latitude_band"],
            "band_j": meta.loc[b, "latitude_band"],
        }

        domain_distances = []
        for domain in domains:
            cols = [f"z_{p}" for p in pred_by_domain[domain]]
            va = std_index.loc[a, cols].to_numpy(dtype=float)
            vb = std_index.loc[b, cols].to_numpy(dtype=float)
            dist = float(np.sqrt(np.sum((va - vb) ** 2)) / np.sqrt(len(cols)))
            row[f"envdist_{domain}"] = dist
            domain_distances.append(dist)

        row["envdist_primary_balanced"] = float(np.mean(domain_distances))
        pair_rows.append(row)

    pairwise = pd.DataFrame(pair_rows)

    # Add an all-predictor standardized Euclidean distance as sensitivity descriptor,
    # not the primary signal.
    all_z_cols = [f"z_{p}" for p in predictors]
    for idx, r in pairwise.iterrows():
        va = std_index.loc[r["cell_i"], all_z_cols].to_numpy(dtype=float)
        vb = std_index.loc[r["cell_j"], all_z_cols].to_numpy(dtype=float)
        pairwise.at[idx, "envdist_all_predictors_std_euclidean"] = float(
            np.sqrt(np.sum((va - vb) ** 2)) / np.sqrt(len(all_z_cols))
        )

    summary_rows = []
    for col in [f"envdist_{d}" for d in domains] + [
        "envdist_primary_balanced",
        "envdist_all_predictors_std_euclidean",
    ]:
        summary_rows.append({
            "distance_metric": col,
            "n_pairs": len(pairwise),
            "mean": float(pairwise[col].mean()),
            "sd": float(pairwise[col].std(ddof=1)),
            "median": float(pairwise[col].median()),
            "q025": float(pairwise[col].quantile(0.025)),
            "q975": float(pairwise[col].quantile(0.975)),
            "min": float(pairwise[col].min()),
            "max": float(pairwise[col].max()),
        })
    distance_summary = pd.DataFrame(summary_rows)

    frozen_cfg_path = outdir / "13C_predictor_config_frozen.csv"
    std_path = outdir / "13C_cell_environment_standardized.csv"
    pair_path = outdir / "13C_pairwise_environment_distance.csv"
    stats_path = outdir / "13C_environment_domain_summary.csv"
    corr_path = outdir / "13C_predictor_correlation.csv"
    dist_summary_path = outdir / "13C_environment_distance_summary.csv"
    manifest_path = outdir / "13C_run_manifest.txt"
    readme_path = outdir / "13C_README.txt"

    shutil.copy2(config_path, frozen_cfg_path)
    standardized.to_csv(std_path, index=False)
    pairwise.to_csv(pair_path, index=False)
    stats.to_csv(stats_path, index=False)
    corr.to_csv(corr_path)
    distance_summary.to_csv(dist_summary_path, index=False)

    readme = """PHASE 13C — FROZEN CONTEMPORARY ENVIRONMENT SIGNAL

Purpose:
Create an outcome-independent contemporary environmental dissimilarity signal
for the retained Baja 25-km cells.

Primary environmental distance:
1. Z-standardize frozen predictors across eligible cells.
2. Compute standardized Euclidean distance separately within thermal, moisture,
   wind, and vegetation domains.
3. Divide each domain distance by sqrt(number of variables in that domain).
4. Average the four domain distances with equal domain weight.

This prevents the moisture domain (which has more predictors) from dominating
the environmental signal simply because it contains more variables.

Predictor selection is not based on C3/N0 Phase 13 outcomes.
ECOSTRESS is not included because the retained Step 12L feasibility analysis
failed coverage thresholds. Existing ERA5-Land/MODIS variables provide broader
and already-audited peninsula-wide coverage.

Do not alter predictor membership after inspecting Phase 13 trait results.
Any alternative predictor sets must be labeled sensitivity analyses.
"""
    readme_path.write_text(readme, encoding="utf-8")

    manifest = [
        "PHASE 13C RUN MANIFEST",
        f"python={sys.version.replace(chr(10), ' ')}",
        f"project_root={root}",
        f"environment_table={env_path}",
        f"environment_table_sha256={sha256(env_path)}",
        f"config_source={config_path}",
        f"config_source_sha256={sha256(config_path)}",
        f"n_cells_input={len(env)}",
        f"n_cells_eligible={len(eligible)}",
        f"n_cell_pairs={len(pairwise)}",
        f"eligibility_rule={eligibility_rule}",
        "primary_environment_signal=equal-weight mean of thermal/moisture/wind/vegetation domain distances",
        "outcome_based_predictor_selection=NO",
        "C3_N0_outcome_testing_performed=NO",
        "ECOSTRESS_included=NO (retained Step 12L coverage failure)",
    ]
    manifest_path.write_text("\n".join(manifest) + "\n", encoding="utf-8")

    print("\nENVIRONMENT DISTANCE SUMMARY")
    print(distance_summary.to_string(index=False))

    print("\nFILES WRITTEN")
    for p in [
        frozen_cfg_path, std_path, pair_path, stats_path, corr_path,
        dist_summary_path, manifest_path, readme_path
    ]:
        print(f"  {p}")

    print("\nSTATUS: PASS")
    print("13C built and locked the outcome-independent contemporary-environment signal.")
    print("No C3/N0 outcome testing was performed.")
    print("Next step: 13D construct paired C3/N0 community dissimilarities on the same cells.")


if __name__ == "__main__":
    main()
