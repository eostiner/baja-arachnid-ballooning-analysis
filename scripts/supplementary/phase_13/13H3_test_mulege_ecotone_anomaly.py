#!/usr/bin/env python3
"""
13H3_test_mulege_ecotone_anomaly.py

Exploratory/post-hoc local-scale test of the apparent Mulegé / central-Gulf
increase in C3 ballooning representation.

The focal location was noticed visually. Therefore:
- all location-specific results are exploratory;
- the focal center/radii and junction predictor were frozen in 13H2;
- no result should replace the confirmatory Phase 13F/13G conclusion.

Primary response
----------------
C3_richness successes out of C3_richness + N0_richness classified genera,
modeled with a binomial-logit GLM.

Broad spatial/sampling baseline
-------------------------------
z(latitude), z(longitude), quadratic latitude/longitude, latitude×longitude,
and z(log1p classified richness).

Primary exploratory additions
-----------------------------
1. Fixed 75-km Mulegé focal indicator.
2. Frozen ecoregion-junction score.
3. Costa Central del Golfo–Sierra de la Giganta overlap/interface flag.
4. Mulegé indicator + junction score jointly.

Uncertainty
-----------
Geographic block bootstrap using 1° latitude × 1° longitude blocks. Bootstrap
percentile intervals and two-sided sign probabilities are reported. Ordinary
iid row-wise P-values are diagnostic only and are not used for conclusions.

Sensitivity checks
------------------
- 50, 100 and 125-km focal radii.
- Minimum classified richness of 5 and 10 genera per cell.
- Spatial + environmental-PC baseline.
- Descriptive spatial scan of 75-km residual windows.
- Exploratory quadratic Gulf-side curve and Mulegé-vs-flanks contrast.

Outputs
-------
13H3_model_results.csv
13H3_block_bootstrap_results.csv
13H3_spatial_scan_75km.csv
13H3_gulf_curve.csv
13H3_gulf_local_contrast.csv
13H3_conclusion_summary.csv
13H3_mulege_residual_scan.png/.svg
13H3_gulf_curve.png/.svg
13H3_README.txt
13H3_run_manifest.txt
"""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd


EPS = 1e-9


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def expit(x):
    x = np.clip(np.asarray(x, dtype=float), -35, 35)
    return 1.0 / (1.0 + np.exp(-x))


def logit(p):
    p = np.clip(np.asarray(p, dtype=float), 1e-7, 1 - 1e-7)
    return np.log(p / (1 - p))


def safe_z(series):
    """
    Z-standardize either a pandas Series or a NumPy-like array.
    """
    numeric = pd.to_numeric(series, errors="coerce")
    x = np.asarray(numeric, dtype=float)
    mu = np.nanmean(x)
    sd = np.nanstd(x, ddof=1)
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("Cannot standardize a constant/nonfinite variable.")
    return (x - mu) / sd, mu, sd


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0088
    p1 = np.radians(np.asarray(lat1, dtype=float))
    p2 = np.radians(np.asarray(lat2, dtype=float))
    dp = p2 - p1
    dl = np.radians(np.asarray(lon2, dtype=float) - np.asarray(lon1, dtype=float))
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * r * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def fit_binomial_irls(X, y, n, max_iter=100, tol=1e-9):
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    n = np.asarray(n, dtype=float)

    if np.any(n <= 0):
        raise ValueError("Binomial denominators must be positive.")
    if np.any(y < 0) or np.any(y > n):
        raise ValueError("Binomial successes must satisfy 0 <= y <= n.")

    p0 = (y + 0.5) / (n + 1.0)
    eta = logit(p0)
    beta = np.linalg.pinv(X) @ eta

    converged = False
    for _ in range(max_iter):
        eta = X @ beta
        p = expit(eta)
        w = np.maximum(n * p * (1 - p), 1e-8)
        z = eta + (y - n * p) / w

        xtwx = X.T @ (w[:, None] * X)
        xtwz = X.T @ (w * z)
        ridge = np.eye(X.shape[1]) * 1e-9
        ridge[0, 0] = 0.0
        beta_new = np.linalg.pinv(xtwx + ridge) @ xtwz

        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            converged = True
            break
        beta = beta_new

    eta = X @ beta
    p = expit(eta)
    w = np.maximum(n * p * (1 - p), 1e-8)
    xtwx = X.T @ (w[:, None] * X)
    cov_base = np.linalg.pinv(xtwx)

    resid_pearson = (y - n * p) / np.sqrt(np.maximum(n * p * (1 - p), 1e-8))
    rank = np.linalg.matrix_rank(X)
    df_resid = max(len(y) - rank, 1)
    phi = max(float(np.sum(resid_pearson ** 2) / df_resid), 1.0)
    cov_quasi = cov_base * phi
    se = np.sqrt(np.maximum(np.diag(cov_quasi), 0))

    loglik = float(np.sum(y * np.log(p + EPS) + (n - y) * np.log(1 - p + EPS)))
    aic = float(-2 * loglik + 2 * X.shape[1])

    return {
        "beta": beta,
        "se_quasi": se,
        "p_fitted": p,
        "eta": eta,
        "resid_pearson": resid_pearson / math.sqrt(phi),
        "phi": phi,
        "loglik": loglik,
        "aic": aic,
        "converged": converged,
        "rank": rank,
        "df_resid": df_resid,
    }


def build_baseline(df, include_environment=False):
    lat_z, lat_mu, lat_sd = safe_z(df["centroid_latitude"])
    lon_z, lon_mu, lon_sd = safe_z(df["centroid_longitude"])
    logn_z, logn_mu, logn_sd = safe_z(np.log1p(df["classified_richness_C3_plus_N0"]))

    names = [
        "Intercept",
        "z_latitude",
        "z_longitude",
        "z_latitude_sq",
        "z_longitude_sq",
        "z_latitude_x_longitude",
        "z_log_classified_richness",
    ]
    cols = [
        np.ones(len(df)),
        lat_z,
        lon_z,
        lat_z ** 2,
        lon_z ** 2,
        lat_z * lon_z,
        logn_z,
    ]

    if include_environment:
        for pc in ["environment_PC1", "environment_PC2"]:
            if pc not in df.columns:
                raise ValueError(f"Missing environmental score: {pc}")
            zpc, _, _ = safe_z(df[pc])
            names.append(pc)
            cols.append(zpc)

    return np.column_stack(cols), names, {
        "lat_mu": lat_mu, "lat_sd": lat_sd,
        "lon_mu": lon_mu, "lon_sd": lon_sd,
        "logn_mu": logn_mu, "logn_sd": logn_sd,
    }


def add_terms(X_base, names_base, df, terms):
    X = X_base.copy()
    names = list(names_base)
    for name, col in terms:
        vals = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(vals).all():
            raise ValueError(f"Nonfinite values in predictor {col}")
        if np.nanstd(vals) <= 0:
            raise ValueError(f"Predictor {col} has no variation.")
        if col == "ecoregion_junction_score_primary":
            vals, _, _ = safe_z(vals)
        X = np.column_stack([X, vals])
        names.append(name)
    return X, names


def block_ids(df):
    lat_block = np.floor(pd.to_numeric(df["centroid_latitude"], errors="coerce")).astype(int)
    lon_block = np.floor(pd.to_numeric(df["centroid_longitude"], errors="coerce")).astype(int)
    return lat_block.astype(str) + "_" + lon_block.astype(str)


def bootstrap_coefficients(df, terms, interest_terms, include_environment,
                           n_boot, seed):
    df = df.reset_index(drop=True).copy()
    blocks = block_ids(df)
    unique_blocks = sorted(blocks.unique())
    block_to_idx = {
        b: np.where(blocks.to_numpy() == b)[0] for b in unique_blocks
    }

    rng = np.random.default_rng(seed)
    collected = {term: [] for term in interest_terms}
    valid = 0

    for _ in range(n_boot):
        sampled_blocks = rng.choice(unique_blocks, size=len(unique_blocks), replace=True)
        idx = np.concatenate([block_to_idx[b] for b in sampled_blocks])
        boot = df.iloc[idx].reset_index(drop=True)

        try:
            Xb, nb, _ = build_baseline(boot, include_environment=include_environment)
            Xb, names = add_terms(Xb, nb, boot, terms)
            fit = fit_binomial_irls(
                Xb,
                boot["C3_richness"].to_numpy(dtype=float),
                boot["classified_richness_C3_plus_N0"].to_numpy(dtype=float),
            )
            bmap = dict(zip(names, fit["beta"]))
            if all(np.isfinite(bmap[t]) for t in interest_terms):
                for t in interest_terms:
                    collected[t].append(float(bmap[t]))
                valid += 1
        except Exception:
            continue

    rows = []
    for term in interest_terms:
        arr = np.asarray(collected[term], dtype=float)
        if len(arr) < max(200, int(0.5 * n_boot)):
            rows.append({
                "term": term,
                "bootstrap_valid": len(arr),
                "bootstrap_ci_low_2p5": np.nan,
                "bootstrap_ci_high_97p5": np.nan,
                "bootstrap_p_two_sided_sign": np.nan,
                "bootstrap_probability_positive": np.nan,
            })
            continue
        prob_pos = (np.sum(arr > 0) + 1) / (len(arr) + 1)
        prob_neg = (np.sum(arr < 0) + 1) / (len(arr) + 1)
        rows.append({
            "term": term,
            "bootstrap_valid": len(arr),
            "bootstrap_ci_low_2p5": float(np.quantile(arr, 0.025)),
            "bootstrap_ci_high_97p5": float(np.quantile(arr, 0.975)),
            "bootstrap_p_two_sided_sign": float(min(1.0, 2 * min(prob_pos, prob_neg))),
            "bootstrap_probability_positive": float(prob_pos),
        })
    return pd.DataFrame(rows), len(unique_blocks)


def run_model_spec(df, spec_name, family, status, min_classified,
                   terms, interest_terms, include_environment,
                   n_boot, seed):
    d = df.loc[
        df["classified_richness_C3_plus_N0"] >= min_classified
    ].dropna(subset=[
        "C3_richness",
        "classified_richness_C3_plus_N0",
        "centroid_latitude",
        "centroid_longitude",
    ] + [c for _, c in terms]).copy()

    if len(d) < 40:
        return [], [], f"SKIP: only {len(d)} eligible cells"

    X0, names0, _ = build_baseline(d, include_environment=include_environment)
    fit0 = fit_binomial_irls(
        X0,
        d["C3_richness"].to_numpy(float),
        d["classified_richness_C3_plus_N0"].to_numpy(float),
    )

    X1, names1 = add_terms(X0, names0, d, terms)
    fit1 = fit_binomial_irls(
        X1,
        d["C3_richness"].to_numpy(float),
        d["classified_richness_C3_plus_N0"].to_numpy(float),
    )

    bmap = dict(zip(names1, fit1["beta"]))
    semap = dict(zip(names1, fit1["se_quasi"]))

    model_rows = []
    for term in interest_terms:
        beta = float(bmap[term])
        se = float(semap[term])
        z = beta / se if se > 0 else np.nan
        p_diag = math.erfc(abs(z) / math.sqrt(2)) if np.isfinite(z) else np.nan
        model_rows.append({
            "spec_name": spec_name,
            "spec_family": family,
            "analysis_status": status,
            "min_classified_richness": min_classified,
            "include_environment_PCs": include_environment,
            "n_cells": len(d),
            "n_focal_50km": int(d["Mulege_focal_sensitivity_50km"].sum()),
            "n_focal_75km": int(d["Mulege_focal_primary_75km"].sum()),
            "n_focal_100km": int(d["Mulege_focal_sensitivity_100km"].sum()),
            "n_focal_125km": int(d["Mulege_focal_sensitivity_125km"].sum()),
            "term": term,
            "beta_log_odds": beta,
            "odds_ratio": float(np.exp(np.clip(beta, -20, 20))),
            "quasibinomial_se_diagnostic": se,
            "quasibinomial_z_diagnostic": z,
            "quasibinomial_p_diagnostic": p_diag,
            "overdispersion_phi": fit1["phi"],
            "AIC_reduced_descriptive": fit0["aic"],
            "AIC_full_descriptive": fit1["aic"],
            "delta_AIC_full_minus_reduced": fit1["aic"] - fit0["aic"],
            "converged": fit1["converged"],
        })

    boot, n_blocks = bootstrap_coefficients(
        d, terms, interest_terms, include_environment,
        n_boot=n_boot, seed=seed
    )
    boot["spec_name"] = spec_name
    boot["spec_family"] = family
    boot["analysis_status"] = status
    boot["min_classified_richness"] = min_classified
    boot["include_environment_PCs"] = include_environment
    boot["n_cells"] = len(d)
    boot["n_geographic_blocks"] = n_blocks
    boot["n_bootstrap_requested"] = n_boot

    return model_rows, boot.to_dict("records"), "PASS"


def compute_environment_pcs(env, predictor_df):
    env = env.copy()
    zcols = [c for c in env.columns if c.startswith("z_")]
    if len(zcols) < 4:
        raise ValueError("13C standardized table lacks sufficient z predictors.")

    env["grid_cell_id"] = env["grid_cell_id"].astype(str)
    x = env[zcols].apply(pd.to_numeric, errors="coerce")
    complete = x.notna().all(axis=1)
    env2 = env.loc[complete, ["grid_cell_id"]].copy()
    X = x.loc[complete].to_numpy(dtype=float)
    X = X - X.mean(axis=0)
    _, _, vt = np.linalg.svd(X, full_matrices=False)
    scores = X @ vt[:2].T
    env2["environment_PC1"] = scores[:, 0]
    env2["environment_PC2"] = scores[:, 1]

    out = predictor_df.merge(env2, on="grid_cell_id", how="left", validate="one_to_one")
    return out, zcols


def spatial_scan(df, radius_km=75.0):
    d = df.loc[
        df["classified_richness_C3_plus_N0"] >= 1
    ].dropna(subset=[
        "C3_richness", "classified_richness_C3_plus_N0",
        "centroid_latitude", "centroid_longitude"
    ]).copy().reset_index(drop=True)

    X0, names0, _ = build_baseline(d, include_environment=False)
    fit0 = fit_binomial_irls(
        X0,
        d["C3_richness"].to_numpy(float),
        d["classified_richness_C3_plus_N0"].to_numpy(float),
    )
    d["baseline_pearson_residual"] = fit0["resid_pearson"]

    lat = d["centroid_latitude"].to_numpy(float)
    lon = d["centroid_longitude"].to_numpy(float)
    rows = []

    for _, center in d.iterrows():
        dist = haversine_km(
            lat, lon,
            float(center["centroid_latitude"]),
            float(center["centroid_longitude"])
        )
        mask = dist <= radius_km
        rows.append({
            "center_grid_cell_id": center["grid_cell_id"],
            "center_latitude": center["centroid_latitude"],
            "center_longitude": center["centroid_longitude"],
            "radius_km": radius_km,
            "n_cells_in_window": int(mask.sum()),
            "mean_baseline_residual": float(d.loc[mask, "baseline_pearson_residual"].mean()),
            "median_baseline_residual": float(d.loc[mask, "baseline_pearson_residual"].median()),
        })

    scan = pd.DataFrame(rows)

    mulege_dist = haversine_km(lat, lon, 26.89, -111.98)
    mulege_mask = mulege_dist <= radius_km
    mulege_n = int(mulege_mask.sum())
    mulege_stat = float(d.loc[mulege_mask, "baseline_pearson_residual"].mean())

    valid_all = scan.loc[scan["n_cells_in_window"] >= 5].copy()
    matched = valid_all.loc[
        (valid_all["n_cells_in_window"] >= max(5, mulege_n - 3))
        & (valid_all["n_cells_in_window"] <= mulege_n + 3)
    ].copy()

    percentile_all = float(
        (valid_all["mean_baseline_residual"] <= mulege_stat).mean()
    ) if len(valid_all) else np.nan
    percentile_matched = float(
        (matched["mean_baseline_residual"] <= mulege_stat).mean()
    ) if len(matched) else np.nan

    summary = {
        "mulege_radius_km": radius_km,
        "mulege_n_cells": mulege_n,
        "mulege_mean_baseline_residual": mulege_stat,
        "candidate_windows_n_ge5": len(valid_all),
        "mulege_percentile_among_all_windows": percentile_all,
        "matched_window_count": len(matched),
        "mulege_percentile_among_similar_n_windows": percentile_matched,
        "note": "Descriptive spatial rarity scan, not a confirmatory P-value.",
    }
    return d, scan, summary


def gulf_curve_and_contrast(df, n_boot=1999, seed=20260718):
    d = df.loc[
        df["overlaps_Costa_Central_del_Golfo"].astype(bool)
        & (df["classified_richness_C3_plus_N0"] >= 1)
    ].copy()

    if len(d) < 15:
        return pd.DataFrame(), pd.DataFrame([{
            "status": "SKIP", "reason": f"Only {len(d)} Costa Central del Golfo cells"
        }])

    lat_z, lat_mu, lat_sd = safe_z(d["centroid_latitude"])
    logn_z, _, _ = safe_z(np.log1p(d["classified_richness_C3_plus_N0"]))
    X = np.column_stack([np.ones(len(d)), lat_z, lat_z ** 2, logn_z])
    fit = fit_binomial_irls(
        X, d["C3_richness"].to_numpy(float),
        d["classified_richness_C3_plus_N0"].to_numpy(float)
    )

    grid_lat = np.linspace(d["centroid_latitude"].min(),
                           d["centroid_latitude"].max(), 150)
    gz = (grid_lat - lat_mu) / lat_sd
    Xg = np.column_stack([np.ones(len(grid_lat)), gz, gz ** 2,
                          np.zeros(len(grid_lat))])
    pred = expit(Xg @ fit["beta"])
    curve = pd.DataFrame({
        "latitude": grid_lat,
        "predicted_C3_fraction_at_mean_log_richness": pred
    })

    def pred_at(beta, lat):
        z = (lat - lat_mu) / lat_sd
        return float(expit(np.array([1.0, z, z*z, 0.0]) @ beta))

    center = 26.89
    south = 26.20
    north = 27.58
    observed_contrast = (
        pred_at(fit["beta"], center)
        - 0.5 * (pred_at(fit["beta"], south) + pred_at(fit["beta"], north))
    )

    blocks = (
        np.floor(d["centroid_latitude"] * 2) / 2
    ).astype(str)
    unique_blocks = sorted(blocks.unique())
    block_to_idx = {b: np.where(blocks.to_numpy() == b)[0] for b in unique_blocks}
    rng = np.random.default_rng(seed)
    vals = []

    for _ in range(n_boot):
        sampled = rng.choice(unique_blocks, size=len(unique_blocks), replace=True)
        idx = np.concatenate([block_to_idx[b] for b in sampled])
        bdf = d.iloc[idx].reset_index(drop=True)
        try:
            bz, bmu, bsd = safe_z(bdf["centroid_latitude"])
            bln, _, _ = safe_z(np.log1p(bdf["classified_richness_C3_plus_N0"]))
            BX = np.column_stack([np.ones(len(bdf)), bz, bz**2, bln])
            bf = fit_binomial_irls(
                BX, bdf["C3_richness"].to_numpy(float),
                bdf["classified_richness_C3_plus_N0"].to_numpy(float)
            )

            def bp(lat):
                z = (lat - bmu) / bsd
                return float(expit(np.array([1.0, z, z*z, 0.0]) @ bf["beta"]))

            vals.append(bp(center) - 0.5 * (bp(south) + bp(north)))
        except Exception:
            continue

    arr = np.asarray(vals, dtype=float)
    if len(arr) >= 200:
        prob_pos = (np.sum(arr > 0) + 1) / (len(arr) + 1)
        p_two = min(1.0, 2 * min(prob_pos, 1 - prob_pos + 1/(len(arr)+1)))
        lo, hi = np.quantile(arr, [0.025, 0.975])
    else:
        prob_pos = p_two = lo = hi = np.nan

    contrast = pd.DataFrame([{
        "status": "PASS",
        "n_Costa_Central_del_Golfo_cells": len(d),
        "center_latitude": center,
        "south_flank_latitude": south,
        "north_flank_latitude": north,
        "predicted_center_minus_mean_flanks": observed_contrast,
        "bootstrap_valid": len(arr),
        "bootstrap_ci_low_2p5": lo,
        "bootstrap_ci_high_97p5": hi,
        "bootstrap_probability_positive": prob_pos,
        "bootstrap_p_two_sided_sign": p_two,
        "analysis_status": "EXPLORATORY_POST_HOC",
    }])
    return curve, contrast


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--predictors", type=Path, default=None)
    p.add_argument("--environment", type=Path, default=None)
    p.add_argument("--output-dir", type=Path, default=None)
    p.add_argument("--bootstrap", type=int, default=1999)
    p.add_argument("--seed", type=int, default=20260718)
    return p.parse_args()


def main():
    args = parse_args()
    root = args.project_root.expanduser().resolve()

    predictor_path = args.predictors or (
        root / "04_analysis_USE _THIS" / "13_historical_vs_contemporary"
        / "13H_mulege_ecotone" / "13H2_ecotone_predictors"
        / "13H2_ecotone_predictors_analysis_cells.csv"
    )
    env_path = args.environment or (
        root / "04_analysis_USE _THIS" / "13_historical_vs_contemporary"
        / "13C_contemporary_environment_signal"
        / "13C_cell_environment_standardized.csv"
    )
    outdir = args.output_dir or (
        root / "04_analysis_USE _THIS" / "13_historical_vs_contemporary"
        / "13H_mulege_ecotone" / "13H3_local_inference"
    )

    predictor_path = predictor_path.expanduser().resolve()
    env_path = env_path.expanduser().resolve()
    outdir = outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print("PHASE 13H3 — MULEGÉ / ECOREGION-JUNCTION LOCAL INFERENCE")
    print("=" * 80)
    print(f"PREDICTORS : {predictor_path}")
    print(f"ENVIRONMENT: {env_path}")
    print(f"OUTPUT DIR : {outdir}")
    print(f"BOOTSTRAP  : {args.bootstrap}")
    print(f"SEED       : {args.seed}")

    for p in [predictor_path, env_path]:
        if not p.exists():
            raise FileNotFoundError(p)

    pred = pd.read_csv(predictor_path, low_memory=False)
    env = pd.read_csv(env_path, low_memory=False)
    pred["grid_cell_id"] = pred["grid_cell_id"].astype(str)

    pred, env_cols = compute_environment_pcs(env, pred)

    specs = [
        {
            "name": "PRIMARY_Mulege_75km",
            "family": "fixed_focal_anomaly",
            "status": "PRIMARY_EXPLORATORY",
            "min_n": 1,
            "terms": [("Mulege_75km", "Mulege_focal_primary_75km")],
            "interest": ["Mulege_75km"],
            "env": False,
        },
        {
            "name": "PRIMARY_junction_score",
            "family": "ecoregion_junction",
            "status": "PRIMARY_EXPLORATORY",
            "min_n": 1,
            "terms": [("junction_score", "ecoregion_junction_score_primary")],
            "interest": ["junction_score"],
            "env": False,
        },
        {
            "name": "PRIMARY_Costa_Golfo_x_Giganta_overlap",
            "family": "ecoregion_interface",
            "status": "PRIMARY_EXPLORATORY",
            "min_n": 1,
            "terms": [("Costa_Golfo_Giganta_overlap",
                       "overlaps_Costa_Golfo_and_Sierra_Giganta")],
            "interest": ["Costa_Golfo_Giganta_overlap"],
            "env": False,
        },
        {
            "name": "PRIMARY_Mulege75_plus_junction",
            "family": "joint_local_and_junction",
            "status": "PRIMARY_EXPLORATORY",
            "min_n": 1,
            "terms": [
                ("Mulege_75km", "Mulege_focal_primary_75km"),
                ("junction_score", "ecoregion_junction_score_primary"),
            ],
            "interest": ["Mulege_75km", "junction_score"],
            "env": False,
        },
        {
            "name": "SENS_radius_50km",
            "family": "focal_radius",
            "status": "SENSITIVITY",
            "min_n": 1,
            "terms": [("Mulege_50km", "Mulege_focal_sensitivity_50km")],
            "interest": ["Mulege_50km"],
            "env": False,
        },
        {
            "name": "SENS_radius_100km",
            "family": "focal_radius",
            "status": "SENSITIVITY",
            "min_n": 1,
            "terms": [("Mulege_100km", "Mulege_focal_sensitivity_100km")],
            "interest": ["Mulege_100km"],
            "env": False,
        },
        {
            "name": "SENS_radius_125km",
            "family": "focal_radius",
            "status": "SENSITIVITY",
            "min_n": 1,
            "terms": [("Mulege_125km", "Mulege_focal_sensitivity_125km")],
            "interest": ["Mulege_125km"],
            "env": False,
        },
        {
            "name": "SENS_Mulege75_min5",
            "family": "denominator_threshold",
            "status": "SENSITIVITY",
            "min_n": 5,
            "terms": [("Mulege_75km", "Mulege_focal_primary_75km")],
            "interest": ["Mulege_75km"],
            "env": False,
        },
        {
            "name": "SENS_junction_min5",
            "family": "denominator_threshold",
            "status": "SENSITIVITY",
            "min_n": 5,
            "terms": [("junction_score", "ecoregion_junction_score_primary")],
            "interest": ["junction_score"],
            "env": False,
        },
        {
            "name": "SENS_Mulege75_min10",
            "family": "denominator_threshold",
            "status": "SENSITIVITY",
            "min_n": 10,
            "terms": [("Mulege_75km", "Mulege_focal_primary_75km")],
            "interest": ["Mulege_75km"],
            "env": False,
        },
        {
            "name": "SENS_junction_min10",
            "family": "denominator_threshold",
            "status": "SENSITIVITY",
            "min_n": 10,
            "terms": [("junction_score", "ecoregion_junction_score_primary")],
            "interest": ["junction_score"],
            "env": False,
        },
        {
            "name": "SENS_Mulege75_environment_PCs",
            "family": "environment_adjusted",
            "status": "SENSITIVITY",
            "min_n": 1,
            "terms": [("Mulege_75km", "Mulege_focal_primary_75km")],
            "interest": ["Mulege_75km"],
            "env": True,
        },
        {
            "name": "SENS_junction_environment_PCs",
            "family": "environment_adjusted",
            "status": "SENSITIVITY",
            "min_n": 1,
            "terms": [("junction_score", "ecoregion_junction_score_primary")],
            "interest": ["junction_score"],
            "env": True,
        },
    ]

    model_rows = []
    boot_rows = []
    status_rows = []

    for i, spec in enumerate(specs):
        mr, br, status = run_model_spec(
            pred,
            spec_name=spec["name"],
            family=spec["family"],
            status=spec["status"],
            min_classified=spec["min_n"],
            terms=spec["terms"],
            interest_terms=spec["interest"],
            include_environment=spec["env"],
            n_boot=args.bootstrap,
            seed=args.seed + i * 1000,
        )
        model_rows.extend(mr)
        boot_rows.extend(br)
        status_rows.append({"spec_name": spec["name"], "run_status": status})

    models = pd.DataFrame(model_rows)
    boots = pd.DataFrame(boot_rows)
    statuses = pd.DataFrame(status_rows)

    residual_df, scan, scan_summary = spatial_scan(pred, radius_km=75.0)
    curve, contrast = gulf_curve_and_contrast(
        pred, n_boot=args.bootstrap, seed=args.seed + 50000
    )

    merged = models.merge(
        boots,
        on=[
            "spec_name", "spec_family", "analysis_status",
            "min_classified_richness", "include_environment_PCs",
            "n_cells", "term"
        ],
        how="left",
        validate="one_to_one"
    )

    def get_row(spec, term):
        x = merged.loc[(merged["spec_name"] == spec) & (merged["term"] == term)]
        return x.iloc[0] if len(x) else None

    focal = get_row("PRIMARY_Mulege_75km", "Mulege_75km")
    junction = get_row("PRIMARY_junction_score", "junction_score")
    interface = get_row(
        "PRIMARY_Costa_Golfo_x_Giganta_overlap",
        "Costa_Golfo_Giganta_overlap"
    )

    def supported(row):
        if row is None:
            return False
        return (
            row["beta_log_odds"] > 0
            and np.isfinite(row["bootstrap_ci_low_2p5"])
            and row["bootstrap_ci_low_2p5"] > 0
        )

    focal_support = supported(focal)
    junction_support = supported(junction)
    interface_support = supported(interface)

    if focal_support and junction_support:
        conclusion = "localized_Mulege_increase_and_junction_association_supported_exploratorily"
    elif focal_support:
        conclusion = "localized_Mulege_increase_supported_but_not_general_junction_effect"
    elif junction_support:
        conclusion = "junction_association_supported_but_fixed_Mulege_anomaly_not_supported"
    elif interface_support:
        conclusion = "specific_ecoregion_interface_supported_without_broader_Mulege_or_junction_signal"
    else:
        conclusion = "no_block_bootstrap_supported_local_bump_or_junction_effect"

    conclusion_df = pd.DataFrame([{
        "analysis_status": "EXPLORATORY_POST_HOC",
        "primary_focal_75km_positive_CI_excludes_zero": focal_support,
        "primary_junction_positive_CI_excludes_zero": junction_support,
        "primary_interface_positive_CI_excludes_zero": interface_support,
        "spatial_scan_Mulege_percentile_all_windows":
            scan_summary["mulege_percentile_among_all_windows"],
        "spatial_scan_Mulege_percentile_similar_n_windows":
            scan_summary["mulege_percentile_among_similar_n_windows"],
        "gulf_curve_local_contrast_positive_CI_excludes_zero": (
            len(contrast)
            and contrast.iloc[0].get("status") == "PASS"
            and np.isfinite(contrast.iloc[0].get("bootstrap_ci_low_2p5", np.nan))
            and contrast.iloc[0]["bootstrap_ci_low_2p5"] > 0
        ),
        "overall_exploratory_conclusion": conclusion,
        "reporting_rule": (
            "These tests are post hoc/exploratory because the Mulege feature "
            "was identified visually. They cannot overturn Phase 13F/13G."
        )
    }])

    models_path = outdir / "13H3_model_results.csv"
    boots_path = outdir / "13H3_block_bootstrap_results.csv"
    merged_path = outdir / "13H3_model_and_bootstrap_summary.csv"
    status_path = outdir / "13H3_spec_run_status.csv"
    scan_path = outdir / "13H3_spatial_scan_75km.csv"
    scan_summary_path = outdir / "13H3_spatial_scan_summary.csv"
    curve_path = outdir / "13H3_gulf_curve.csv"
    contrast_path = outdir / "13H3_gulf_local_contrast.csv"
    conclusion_path = outdir / "13H3_conclusion_summary.csv"
    readme_path = outdir / "13H3_README.txt"
    manifest_path = outdir / "13H3_run_manifest.txt"

    models.to_csv(models_path, index=False)
    boots.to_csv(boots_path, index=False)
    merged.to_csv(merged_path, index=False)
    statuses.to_csv(status_path, index=False)
    scan.to_csv(scan_path, index=False)
    pd.DataFrame([scan_summary]).to_csv(scan_summary_path, index=False)
    curve.to_csv(curve_path, index=False)
    contrast.to_csv(contrast_path, index=False)
    conclusion_df.to_csv(conclusion_path, index=False)

    scan_png = outdir / "13H3_mulege_residual_scan.png"
    scan_svg = outdir / "13H3_mulege_residual_scan.svg"
    curve_png = outdir / "13H3_gulf_curve.png"
    curve_svg = outdir / "13H3_gulf_curve.svg"

    try:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 9))
        sc = ax.scatter(
            residual_df["centroid_longitude"],
            residual_df["centroid_latitude"],
            c=residual_df["baseline_pearson_residual"],
            s=40,
        )
        ax.scatter([-111.98], [26.89], marker="*", s=180, label="Heroica Mulegé")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title("C3 representation residuals after broad spatial baseline")
        ax.legend()
        fig.colorbar(sc, ax=ax, label="Standardized Pearson residual")
        fig.tight_layout()
        fig.savefig(scan_png, dpi=300, bbox_inches="tight")
        fig.savefig(scan_svg, bbox_inches="tight")
        plt.close(fig)

        if len(curve):
            gulf = pred.loc[
                pred["overlaps_Costa_Central_del_Golfo"].astype(bool)
                & (pred["classified_richness_C3_plus_N0"] >= 1)
            ].copy()
            fig, ax = plt.subplots(figsize=(9, 6))
            sizes = np.maximum(
                20,
                np.sqrt(gulf["classified_richness_C3_plus_N0"]) * 10
            )
            ax.scatter(
                gulf["centroid_latitude"], gulf["C3_fraction"], s=sizes,
                alpha=0.65, label="Occupied cells"
            )
            ax.plot(
                curve["latitude"],
                curve["predicted_C3_fraction_at_mean_log_richness"],
                linewidth=2, label="Quadratic binomial fit"
            )
            ax.axvline(26.89, linestyle="--", label="Mulegé latitude")
            ax.set_xlabel("Latitude")
            ax.set_ylabel("C3 representation")
            ax.set_title("Exploratory Costa Central del Golfo C3 curve")
            ax.legend()
            fig.tight_layout()
            fig.savefig(curve_png, dpi=300, bbox_inches="tight")
            fig.savefig(curve_svg, bbox_inches="tight")
            plt.close(fig)
    except Exception as e:
        print(f"WARNING: figure creation failed: {e}")

    readme_path.write_text(
        """PHASE 13H3 — EXPLORATORY MULEGÉ / ECOREGION-JUNCTION TEST

This is a post-hoc exploratory analysis because the apparent Mulegé feature
was noticed visually. The location, radii and junction score were frozen in
13H2 before this test.

Primary response:
C3 genera out of C3 + N0 classified genera per cell.

Primary baseline:
quadratic latitude/longitude spatial surface plus classified richness.

Uncertainty:
1-degree geographic block bootstrap. Ordinary iid quasibinomial P-values are
diagnostic only. Interpret bootstrap intervals first.

Primary exploratory evidence requires a positive coefficient whose geographic
block-bootstrap 95% interval excludes zero.

The residual-window scan is descriptive and reports the percentile of the
Mulegé 75-km window among equivalent windows across the occupied peninsula.
It is not a confirmatory P-value.

The Gulf curve and local center-vs-flanks contrast are also exploratory.
None of Phase 13H can overturn the confirmatory Phase 13F/13G null mechanism
test.
""",
        encoding="utf-8"
    )

    manifest_path.write_text(
        "\n".join([
            "PHASE 13H3 RUN MANIFEST",
            f"python={sys.version.replace(chr(10), ' ')}",
            f"project_root={root}",
            f"predictors={predictor_path}",
            f"predictors_sha256={sha256(predictor_path)}",
            f"environment={env_path}",
            f"environment_sha256={sha256(env_path)}",
            f"bootstrap_requested={args.bootstrap}",
            f"seed={args.seed}",
            f"environment_z_columns={';'.join(env_cols)}",
            "analysis_status=EXPLORATORY_POST_HOC",
            "primary_focal_radius_km=75",
            "spatial_block_size=1_degree_latitude_x_1_degree_longitude",
            "ordinary_iid_p_values_used_for_conclusion=NO",
        ]) + "\n",
        encoding="utf-8"
    )

    print("\nPRIMARY MODEL + BLOCK BOOTSTRAP RESULTS")
    show = merged.loc[
        merged["analysis_status"] == "PRIMARY_EXPLORATORY",
        [
            "spec_name", "term", "n_cells", "beta_log_odds", "odds_ratio",
            "bootstrap_ci_low_2p5", "bootstrap_ci_high_97p5",
            "bootstrap_p_two_sided_sign", "bootstrap_probability_positive",
            "n_geographic_blocks", "bootstrap_valid"
        ]
    ]
    print(show.to_string(index=False))

    print("\nSPATIAL SCAN SUMMARY")
    print(pd.DataFrame([scan_summary]).to_string(index=False))

    print("\nGULF LOCAL CONTRAST")
    print(contrast.to_string(index=False))

    print("\nCONCLUSION SUMMARY")
    print(conclusion_df.to_string(index=False))

    print("\nFILES WRITTEN")
    for p in [
        models_path, boots_path, merged_path, status_path,
        scan_path, scan_summary_path, curve_path, contrast_path,
        conclusion_path, scan_png, scan_svg, curve_png, curve_svg,
        readme_path, manifest_path
    ]:
        if p.exists():
            print(f"  {p}")

    print("\nSTATUS: PASS")
    print("13H3 completed the exploratory local Mulegé/ecotone tests.")
    print("Review 13H3_conclusion_summary.csv first.")
    print("Next: decide whether the local pattern is unsupported, suggestive, or supported exploratorily.")


if __name__ == "__main__":
    main()
