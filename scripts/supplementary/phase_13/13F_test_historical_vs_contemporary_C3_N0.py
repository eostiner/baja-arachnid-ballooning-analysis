#!/usr/bin/env python3
"""
13F_test_historical_vs_contemporary_C3_N0.py

Purpose
-------
Primary Phase 13 inference.

Question:
Do ballooning-capable assemblages (C3 = D1+D2+D3) show a weaker relative
association with frozen historical biogeographic boundaries and/or a stronger
relative association with contemporary environmental dissimilarity than
non-ballooning assemblages (N0)?

Primary response:
  delta_simpson = C3 Simpson replacement - N0 Simpson replacement

Secondary response:
  delta_jaccard = C3 Jaccard total dissimilarity - N0 Jaccard total dissimilarity

Primary cell set:
Cells where BOTH C3 richness > 0 and N0 richness > 0, intersected with the
189 Phase 13C eligible cells. This yields complete C3 and N0 dissimilarity
matrices for node-label permutation and avoids interpreting empty assemblages
as replacement signal.

Primary model:
  delta ~ z(log1p(geographic_distance_km))
          + z(envdist_primary_balanced)
          + B01_strict_cross
          + B03_strict_cross

Interpretation of paired-difference coefficients:
  environment beta > 0:
      C3 dissimilarity increases more strongly with contemporary environmental
      distance than N0 (direction expected under stronger modern tracking).

  B01/B03 beta < 0:
      C3 shows a weaker across-boundary discontinuity than N0 (direction
      expected under weaker retention of historical-boundary structure).

Inference:
Freedman-Lane-style residual permutation with simultaneous row/column
(cell-label) permutation of the residual matrix. This preserves dyadic/matrix
dependence. P-values are permutation-based, not ordinary iid regression
P-values.

Tests:
1. Contemporary environment incremental effect (1 df), controlling geography
   and B01/B03.
2. Historical boundaries joint incremental effect (2 df), controlling geography
   and environment.
3. B01 individual coefficient, controlling all other primary predictors.
4. B03 individual coefficient, controlling all other primary predictors.

Primary metric = Simpson replacement.
Jaccard = secondary.
All tests are two-sided for coefficients; joint historical F test is upper-tail.

Outputs
-------
13F_analysis_cells_primary.csv
13F_primary_pair_table.csv
13F_model_coefficients.csv
13F_permutation_tests.csv
13F_trait_specific_descriptive_models.csv
13F_collinearity_diagnostics.csv
13F_hypothesis_summary.csv
13F_run_manifest.txt
13F_README.txt
"""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PRIMARY_BOUNDARIES = ["B01_strict_cross", "B03_strict_cross"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args():
    p = argparse.ArgumentParser(
        description="Phase 13F paired historical-vs-contemporary permutation test."
    )
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--master", type=Path, default=None)
    p.add_argument("--cell-richness", type=Path, default=None)
    p.add_argument("--output-dir", type=Path, default=None)
    p.add_argument("--permutations", type=int, default=4999)
    p.add_argument("--seed", type=int, default=20260718)
    return p.parse_args()


def zscore(x: np.ndarray):
    x = np.asarray(x, dtype=float)
    mu = float(np.mean(x))
    sd = float(np.std(x, ddof=1))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("Cannot z-standardize predictor with zero/nonfinite SD.")
    return (x - mu) / sd, mu, sd


def prepare_ols(X: np.ndarray):
    X = np.asarray(X, dtype=float)
    xtx = X.T @ X
    inv = np.linalg.pinv(xtx)
    beta_map = inv @ X.T
    rank = np.linalg.matrix_rank(X)
    df_resid = X.shape[0] - rank
    if df_resid <= 0:
        raise ValueError("Nonpositive regression residual degrees of freedom.")
    return {
        "X": X,
        "inv": inv,
        "beta_map": beta_map,
        "rank": rank,
        "df_resid": df_resid,
    }


def fit_prepared(prep, y: np.ndarray):
    X = prep["X"]
    y = np.asarray(y, dtype=float)
    beta = prep["beta_map"] @ y
    fitted = X @ beta
    resid = y - fitted
    rss = float(resid @ resid)
    tss = float(((y - np.mean(y)) ** 2).sum())
    r2 = np.nan if tss <= 0 else 1.0 - rss / tss
    s2 = rss / prep["df_resid"]
    se = np.sqrt(np.maximum(np.diag(prep["inv"]) * s2, 0.0))
    t = np.divide(beta, se, out=np.full_like(beta, np.nan), where=se > 0)
    return {
        "beta": beta,
        "fitted": fitted,
        "resid": resid,
        "rss": rss,
        "r2": r2,
        "se": se,
        "t": t,
        "df_resid": prep["df_resid"],
    }


def f_compare(reduced_fit, full_fit, q: int):
    rss_r = reduced_fit["rss"]
    rss_f = full_fit["rss"]
    df_f = full_fit["df_resid"]
    if q <= 0 or rss_f <= 0:
        return np.nan
    num = max(rss_r - rss_f, 0.0) / q
    den = rss_f / df_f
    return np.nan if den <= 0 else num / den


def partial_r2(reduced_fit, full_fit):
    rss_r = reduced_fit["rss"]
    rss_f = full_fit["rss"]
    if rss_r <= 0:
        return np.nan
    return max(0.0, (rss_r - rss_f) / rss_r)


def make_symmetric_matrix(values, a_idx, b_idx, n):
    m = np.zeros((n, n), dtype=float)
    m[a_idx, b_idx] = values
    m[b_idx, a_idx] = values
    np.fill_diagonal(m, 0.0)
    return m


def vif_table(X_no_intercept: pd.DataFrame):
    rows = []
    cols = list(X_no_intercept.columns)
    arr = X_no_intercept.to_numpy(dtype=float)

    for j, name in enumerate(cols):
        y = arr[:, j]
        other_idx = [k for k in range(arr.shape[1]) if k != j]
        if not other_idx:
            vif = 1.0
            r2 = 0.0
        else:
            Xo = np.column_stack([np.ones(len(y)), arr[:, other_idx]])
            prep = prepare_ols(Xo)
            fit = fit_prepared(prep, y)
            r2 = fit["r2"]
            vif = np.inf if r2 >= 1.0 else 1.0 / (1.0 - r2)
        rows.append({"predictor": name, "vif": vif, "r2_against_other_predictors": r2})

    return pd.DataFrame(rows)


def freedman_lane_term_test(
    y,
    X_full,
    full_names,
    reduced_keep_names,
    tested_name,
    cell_n,
    tri_i,
    tri_j,
    permutations,
    rng,
):
    """
    Single-coefficient Freedman-Lane-style node-label permutation test.
    """
    full_names = list(full_names)
    reduced_keep_names = list(reduced_keep_names)

    keep_idx = [full_names.index(n) for n in reduced_keep_names]
    test_idx = full_names.index(tested_name)

    prep_f = prepare_ols(X_full)
    fit_f = fit_prepared(prep_f, y)

    X_r = X_full[:, keep_idx]
    prep_r = prepare_ols(X_r)
    fit_r = fit_prepared(prep_r, y)

    obs_t = float(fit_f["t"][test_idx])
    obs_beta = float(fit_f["beta"][test_idx])

    fitted_r = fit_r["fitted"]
    resid_mat = make_symmetric_matrix(fit_r["resid"], tri_i, tri_j, cell_n)

    extreme = 0
    perm_t_sum = 0.0
    perm_t_sq = 0.0

    for _ in range(permutations):
        perm = rng.permutation(cell_n)
        rp = resid_mat[np.ix_(perm, perm)][tri_i, tri_j]
        y_star = fitted_r + rp
        fit_star = fit_prepared(prep_f, y_star)
        t_star = float(fit_star["t"][test_idx])

        if np.isfinite(t_star):
            perm_t_sum += t_star
            perm_t_sq += t_star * t_star
            if abs(t_star) >= abs(obs_t) - 1e-12:
                extreme += 1

    p_two = (extreme + 1.0) / (permutations + 1.0)
    mean_t = perm_t_sum / permutations
    var_t = max(perm_t_sq / permutations - mean_t * mean_t, 0.0)

    return {
        "test": tested_name,
        "observed_beta": obs_beta,
        "observed_t": obs_t,
        "permutation_p_two_sided": p_two,
        "n_permutations": permutations,
        "null_t_mean": mean_t,
        "null_t_sd": math.sqrt(var_t),
        "reduced_model": " + ".join(reduced_keep_names),
    }


def freedman_lane_joint_test(
    y,
    X_full,
    full_names,
    reduced_keep_names,
    tested_group_name,
    q,
    cell_n,
    tri_i,
    tri_j,
    permutations,
    rng,
):
    """
    Joint incremental F test via Freedman-Lane-style node-label permutation.
    """
    full_names = list(full_names)
    reduced_keep_names = list(reduced_keep_names)
    keep_idx = [full_names.index(n) for n in reduced_keep_names]

    prep_f = prepare_ols(X_full)
    fit_f = fit_prepared(prep_f, y)

    X_r = X_full[:, keep_idx]
    prep_r = prepare_ols(X_r)
    fit_r = fit_prepared(prep_r, y)

    obs_f = f_compare(fit_r, fit_f, q=q)
    pr2 = partial_r2(fit_r, fit_f)

    fitted_r = fit_r["fitted"]
    resid_mat = make_symmetric_matrix(fit_r["resid"], tri_i, tri_j, cell_n)

    extreme = 0
    perm_f_sum = 0.0
    perm_f_sq = 0.0

    for _ in range(permutations):
        perm = rng.permutation(cell_n)
        rp = resid_mat[np.ix_(perm, perm)][tri_i, tri_j]
        y_star = fitted_r + rp

        fit_f_star = fit_prepared(prep_f, y_star)
        fit_r_star = fit_prepared(prep_r, y_star)
        f_star = f_compare(fit_r_star, fit_f_star, q=q)

        if np.isfinite(f_star):
            perm_f_sum += f_star
            perm_f_sq += f_star * f_star
            if f_star >= obs_f - 1e-12:
                extreme += 1

    p_upper = (extreme + 1.0) / (permutations + 1.0)
    mean_f = perm_f_sum / permutations
    var_f = max(perm_f_sq / permutations - mean_f * mean_f, 0.0)

    return {
        "test": tested_group_name,
        "observed_F": obs_f,
        "df_test": q,
        "partial_R2": pr2,
        "permutation_p_upper": p_upper,
        "n_permutations": permutations,
        "null_F_mean": mean_f,
        "null_F_sd": math.sqrt(var_f),
        "reduced_model": " + ".join(reduced_keep_names),
    }


def main():
    args = parse_args()
    if args.permutations < 99:
        raise ValueError("Use at least 99 permutations; 4999 is the planned primary run.")

    root = args.project_root.expanduser().resolve()

    master_path = args.master or (
        root
        / "04_analysis_USE _THIS"
        / "13_historical_vs_contemporary"
        / "13E_pairwise_master"
        / "13E_pairwise_master_all_pairs.csv"
    )
    richness_path = args.cell_richness or (
        root
        / "04_analysis_USE _THIS"
        / "13_historical_vs_contemporary"
        / "13D_paired_C3_N0_community_dissimilarity"
        / "13D_cell_trait_richness.csv"
    )
    outdir = args.output_dir or (
        root
        / "04_analysis_USE _THIS"
        / "13_historical_vs_contemporary"
        / "13F_historical_vs_contemporary_test"
    )

    master_path = master_path.expanduser().resolve()
    richness_path = richness_path.expanduser().resolve()
    outdir = outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print("PHASE 13F — TEST HISTORICAL VS CONTEMPORARY SIGNAL")
    print("=" * 80)
    print(f"PROJECT ROOT : {root}")
    print(f"MASTER PAIRS : {master_path}")
    print(f"CELL RICHNESS: {richness_path}")
    print(f"OUTPUT DIR   : {outdir}")
    print(f"PERMUTATIONS : {args.permutations}")
    print(f"SEED         : {args.seed}")

    for p in [master_path, richness_path]:
        if not p.exists():
            raise FileNotFoundError(f"Required input not found: {p}")

    master = pd.read_csv(master_path)
    richness = pd.read_csv(richness_path)

    required_master = {
        "cell_a", "cell_b",
        "geographic_distance_km",
        "envdist_primary_balanced",
        "B01_strict_cross", "B03_strict_cross",
        "C3_jaccard", "N0_jaccard",
        "C3_simpson_replacement", "N0_simpson_replacement",
    }
    missing = required_master - set(master.columns)
    if missing:
        raise ValueError(f"13E master missing columns: {sorted(missing)}")

    required_rich = {
        "grid_cell_id", "C3_richness", "N0_richness",
        "centroid_latitude", "centroid_longitude"
    }
    missing_r = required_rich - set(richness.columns)
    if missing_r:
        raise ValueError(f"13D cell richness missing columns: {sorted(missing_r)}")

    richness["grid_cell_id"] = richness["grid_cell_id"].astype(str)
    master["cell_a"] = master["cell_a"].astype(str)
    master["cell_b"] = master["cell_b"].astype(str)

    # Primary complete-matrix cell set: both trait assemblages observed.
    primary_cells = richness.loc[
        (richness["C3_richness"] > 0) & (richness["N0_richness"] > 0)
    ].copy()
    primary_ids = sorted(primary_cells["grid_cell_id"].tolist())
    primary_set = set(primary_ids)

    if len(primary_ids) < 50:
        raise ValueError(
            f"Only {len(primary_ids)} cells have both C3 and N0 richness > 0; "
            "unexpectedly small primary cell set."
        )

    pairs = master.loc[
        master["cell_a"].isin(primary_set) & master["cell_b"].isin(primary_set)
    ].copy()

    expected_pairs = len(primary_ids) * (len(primary_ids) - 1) // 2
    if len(pairs) != expected_pairs:
        raise ValueError(
            f"Primary pair table has {len(pairs)} rows; expected complete matrix "
            f"of {expected_pairs} pairs from {len(primary_ids)} cells."
        )

    outcome_cols = [
        "C3_jaccard", "N0_jaccard",
        "C3_simpson_replacement", "N0_simpson_replacement"
    ]
    if pairs[outcome_cols].isna().any().any():
        raise ValueError(
            "Primary both-traits-present cell set still contains undefined "
            "community dissimilarities; stop and audit."
        )

    print("\nPRIMARY COMPLETE-MATRIX CELL SET")
    print(f"  Phase 13C eligible cells              : {richness['grid_cell_id'].nunique()}")
    print(f"  Cells with C3 > 0 AND N0 > 0         : {len(primary_ids)}")
    print(f"  Complete unordered pairs             : {len(pairs)}")

    # Canonical cell order and pair indices.
    cell_to_idx = {c: i for i, c in enumerate(primary_ids)}
    a_idx = pairs["cell_a"].map(cell_to_idx).to_numpy(dtype=int)
    b_idx = pairs["cell_b"].map(cell_to_idx).to_numpy(dtype=int)

    # Make sure every unordered pair appears exactly once.
    pair_keys = {
        tuple(sorted((a, b)))
        for a, b in zip(pairs["cell_a"], pairs["cell_b"])
    }
    if len(pair_keys) != expected_pairs:
        raise ValueError("Duplicate or missing unordered pairs in primary pair table.")

    # Continuous predictors: freeze standardization on primary pair set.
    log_geo = np.log1p(pairs["geographic_distance_km"].to_numpy(dtype=float))
    z_geo, geo_mu, geo_sd = zscore(log_geo)
    env = pairs["envdist_primary_balanced"].to_numpy(dtype=float)
    z_env, env_mu, env_sd = zscore(env)

    b01 = pairs["B01_strict_cross"].to_numpy(dtype=float)
    b03 = pairs["B03_strict_cross"].to_numpy(dtype=float)

    full_names = ["Intercept", "z_log_geographic_distance", "z_environment", "B01_LaPaz", "B03_Vizcaino"]
    X_full = np.column_stack([
        np.ones(len(pairs)),
        z_geo,
        z_env,
        b01,
        b03,
    ])

    # Collinearity diagnostics.
    Xdiag = pd.DataFrame({
        "z_log_geographic_distance": z_geo,
        "z_environment": z_env,
        "B01_LaPaz": b01,
        "B03_Vizcaino": b03,
    })
    vifs = vif_table(Xdiag)
    corr = Xdiag.corr()
    condition_number = float(np.linalg.cond(X_full))

    print("\nPRIMARY MODEL")
    print("  delta ~ z(log1p(geographic distance)) + z(environment)")
    print("          + B01 La Paz strict crossing + B03 Vizcaino strict crossing")
    print("\nCOLLINEARITY DIAGNOSTICS")
    print(vifs.to_string(index=False))
    print(f"  design condition number: {condition_number:.3f}")

    # Outcomes.
    pairs["delta_simpson_C3_minus_N0"] = (
        pairs["C3_simpson_replacement"] - pairs["N0_simpson_replacement"]
    )
    pairs["delta_jaccard_C3_minus_N0"] = (
        pairs["C3_jaccard"] - pairs["N0_jaccard"]
    )

    metrics = [
        {
            "metric": "simpson_replacement",
            "status": "PRIMARY",
            "delta_col": "delta_simpson_C3_minus_N0",
            "c3_col": "C3_simpson_replacement",
            "n0_col": "N0_simpson_replacement",
        },
        {
            "metric": "jaccard_total",
            "status": "SECONDARY",
            "delta_col": "delta_jaccard_C3_minus_N0",
            "c3_col": "C3_jaccard",
            "n0_col": "N0_jaccard",
        },
    ]

    coeff_rows = []
    perm_rows = []
    trait_rows = []
    hypothesis_rows = []

    base_rng = np.random.default_rng(args.seed)

    for m_idx, spec in enumerate(metrics):
        metric = spec["metric"]
        status = spec["status"]
        y = pairs[spec["delta_col"]].to_numpy(dtype=float)

        prep_f = prepare_ols(X_full)
        fit_f = fit_prepared(prep_f, y)

        # Observed coefficient table.
        for name, beta, se, t in zip(
            full_names, fit_f["beta"], fit_f["se"], fit_f["t"]
        ):
            coeff_rows.append({
                "metric": metric,
                "analysis_status": status,
                "response": "C3_minus_N0",
                "term": name,
                "beta": float(beta),
                "ordinary_ols_se_for_studentization": float(se),
                "ordinary_ols_t_for_studentization": float(t),
                "model_R2": float(fit_f["r2"]),
                "n_cells": len(primary_ids),
                "n_pairs": len(pairs),
            })

        # Trait-specific descriptive fits on exactly the same pairs.
        for trait, col in [("C3", spec["c3_col"]), ("N0", spec["n0_col"])]:
            yt = pairs[col].to_numpy(dtype=float)
            ft = fit_prepared(prep_f, yt)

            # Reduced models for descriptive incremental partial R2.
            idx_no_env = [0, 1, 3, 4]
            idx_no_hist = [0, 1, 2]
            fr_env = fit_prepared(prepare_ols(X_full[:, idx_no_env]), yt)
            fr_hist = fit_prepared(prepare_ols(X_full[:, idx_no_hist]), yt)

            env_pr2 = partial_r2(fr_env, ft)
            hist_pr2 = partial_r2(fr_hist, ft)

            for name, beta in zip(full_names, ft["beta"]):
                trait_rows.append({
                    "metric": metric,
                    "analysis_status": status,
                    "trait_class": trait,
                    "term": name,
                    "beta_descriptive": float(beta),
                    "model_R2": float(ft["r2"]),
                    "environment_partial_R2_descriptive": env_pr2,
                    "historical_B01_B03_partial_R2_descriptive": hist_pr2,
                    "n_cells": len(primary_ids),
                    "n_pairs": len(pairs),
                    "note": "Descriptive coefficients; primary inference is on paired C3-N0 delta.",
                })

        # Independent reproducible RNG stream per metric/test.
        def rng_for(offset):
            return np.random.default_rng(args.seed + m_idx * 100 + offset)

        # Environment coefficient test controlling geography + both boundaries.
        env_test = freedman_lane_term_test(
            y=y,
            X_full=X_full,
            full_names=full_names,
            reduced_keep_names=["Intercept", "z_log_geographic_distance", "B01_LaPaz", "B03_Vizcaino"],
            tested_name="z_environment",
            cell_n=len(primary_ids),
            tri_i=a_idx,
            tri_j=b_idx,
            permutations=args.permutations,
            rng=rng_for(1),
        )
        env_test.update({
            "metric": metric,
            "analysis_status": status,
            "effect_family": "contemporary_environment",
            "expected_direction_for_hypothesis": "positive",
        })
        perm_rows.append(env_test)

        # Joint historical boundary test controlling geography + environment.
        hist_test = freedman_lane_joint_test(
            y=y,
            X_full=X_full,
            full_names=full_names,
            reduced_keep_names=["Intercept", "z_log_geographic_distance", "z_environment"],
            tested_group_name="B01_B03_joint_historical",
            q=2,
            cell_n=len(primary_ids),
            tri_i=a_idx,
            tri_j=b_idx,
            permutations=args.permutations,
            rng=rng_for(2),
        )
        hist_test.update({
            "metric": metric,
            "analysis_status": status,
            "effect_family": "historical_boundaries",
            "expected_direction_for_hypothesis": "negative C3-N0 boundary coefficients",
        })
        perm_rows.append(hist_test)

        # Individual boundary coefficient tests.
        b01_test = freedman_lane_term_test(
            y=y,
            X_full=X_full,
            full_names=full_names,
            reduced_keep_names=["Intercept", "z_log_geographic_distance", "z_environment", "B03_Vizcaino"],
            tested_name="B01_LaPaz",
            cell_n=len(primary_ids),
            tri_i=a_idx,
            tri_j=b_idx,
            permutations=args.permutations,
            rng=rng_for(3),
        )
        b01_test.update({
            "metric": metric,
            "analysis_status": status,
            "effect_family": "historical_boundary_individual",
            "expected_direction_for_hypothesis": "negative",
        })
        perm_rows.append(b01_test)

        b03_test = freedman_lane_term_test(
            y=y,
            X_full=X_full,
            full_names=full_names,
            reduced_keep_names=["Intercept", "z_log_geographic_distance", "z_environment", "B01_LaPaz"],
            tested_name="B03_Vizcaino",
            cell_n=len(primary_ids),
            tri_i=a_idx,
            tri_j=b_idx,
            permutations=args.permutations,
            rng=rng_for(4),
        )
        b03_test.update({
            "metric": metric,
            "analysis_status": status,
            "effect_family": "historical_boundary_individual",
            "expected_direction_for_hypothesis": "negative",
        })
        perm_rows.append(b03_test)

        # Directional hypothesis summary, without converting directional
        # expectations into one-sided p-values.
        beta_map = dict(zip(full_names, fit_f["beta"]))
        env_p = env_test["permutation_p_two_sided"]
        hist_p = hist_test["permutation_p_upper"]

        env_direction = beta_map["z_environment"] > 0
        b01_direction = beta_map["B01_LaPaz"] < 0
        b03_direction = beta_map["B03_Vizcaino"] < 0

        if env_direction and env_p < 0.05 and hist_p < 0.05 and (b01_direction or b03_direction):
            support = "strong_joint_support"
        elif env_direction and env_p < 0.05:
            support = "supports_stronger_contemporary_environment_association_only"
        elif hist_p < 0.05 and (b01_direction or b03_direction):
            support = "supports_weaker_historical_boundary_association_only"
        elif env_direction and (b01_direction or b03_direction):
            support = "directionally_consistent_but_not_jointly_significant"
        else:
            support = "not_directionally_consistent_with_full_hypothesis"

        hypothesis_rows.append({
            "metric": metric,
            "analysis_status": status,
            "n_cells": len(primary_ids),
            "n_pairs": len(pairs),
            "beta_environment_C3_minus_N0": float(beta_map["z_environment"]),
            "environment_expected_positive": bool(env_direction),
            "environment_permutation_p_two_sided": env_p,
            "beta_B01_LaPaz_C3_minus_N0": float(beta_map["B01_LaPaz"]),
            "B01_expected_negative": bool(b01_direction),
            "B01_permutation_p_two_sided": b01_test["permutation_p_two_sided"],
            "beta_B03_Vizcaino_C3_minus_N0": float(beta_map["B03_Vizcaino"]),
            "B03_expected_negative": bool(b03_direction),
            "B03_permutation_p_two_sided": b03_test["permutation_p_two_sided"],
            "historical_joint_partial_R2": hist_test["partial_R2"],
            "historical_joint_permutation_p": hist_p,
            "interpretive_class": support,
        })

    coeff_df = pd.DataFrame(coeff_rows)
    perm_df = pd.DataFrame(perm_rows)
    trait_df = pd.DataFrame(trait_rows)
    hypothesis_df = pd.DataFrame(hypothesis_rows)

    # Save primary pair table with standardized predictors.
    pairs_out = pairs.copy()
    pairs_out["z_log_geographic_distance"] = z_geo
    pairs_out["z_environment"] = z_env

    cells_path = outdir / "13F_analysis_cells_primary.csv"
    pair_path = outdir / "13F_primary_pair_table.csv"
    coeff_path = outdir / "13F_model_coefficients.csv"
    perm_path = outdir / "13F_permutation_tests.csv"
    trait_path = outdir / "13F_trait_specific_descriptive_models.csv"
    vif_path = outdir / "13F_collinearity_diagnostics.csv"
    corr_path = outdir / "13F_predictor_correlation.csv"
    hyp_path = outdir / "13F_hypothesis_summary.csv"
    manifest_path = outdir / "13F_run_manifest.txt"
    readme_path = outdir / "13F_README.txt"

    primary_cells.sort_values("grid_cell_id").to_csv(cells_path, index=False)
    pairs_out.to_csv(pair_path, index=False)
    coeff_df.to_csv(coeff_path, index=False)
    perm_df.to_csv(perm_path, index=False)
    trait_df.to_csv(trait_path, index=False)
    vifs.to_csv(vif_path, index=False)
    corr.to_csv(corr_path)
    hypothesis_df.to_csv(hyp_path, index=False)

    readme = f"""PHASE 13F — HISTORICAL VS CONTEMPORARY TEST

PRIMARY QUESTION
Do C3 ballooning assemblages show stronger relative contemporary-environment
tracking and weaker relative historical-boundary structure than N0
non-ballooning assemblages?

PRIMARY RESPONSE
C3 minus N0 Simpson replacement dissimilarity.

SECONDARY RESPONSE
C3 minus N0 Jaccard total dissimilarity.

PRIMARY CELL SET
Only cells with at least one C3 genus and at least one N0 genus are retained.
This creates complete matched dissimilarity matrices and avoids interpreting
empty assemblages as Simpson replacement.

MODEL
delta ~ z(log1p geographic distance) + z(contemporary environmental distance)
        + B01 La Paz strict crossing + B03 Vizcaino strict crossing

INFERENCE
Freedman-Lane-style residual permutation with simultaneous row/column
(cell-label) permutation. Pairwise rows are NOT treated as iid observations.
Permutation p-values are the inferential p-values.

PERMUTATIONS
{args.permutations}
Seed: {args.seed}

DIRECTIONAL EXPECTATIONS
- Environment coefficient > 0: C3 is relatively more associated with modern
  environmental dissimilarity than N0.
- Boundary coefficient < 0: C3 is relatively less discontinuous across that
  historical boundary than N0.

All reported coefficient permutation tests are two-sided. Direction is
interpreted from the sign of the observed coefficient.

IMPORTANT
The primary historical hypothesis uses only B01 (La Paz) and B03
(Vizcaino/mid-peninsular), frozen before Phase 13 outcome testing.
B02 Loreto and B04 northern transition are not included in the primary model.
"""
    readme_path.write_text(readme, encoding="utf-8")

    manifest = [
        "PHASE 13F RUN MANIFEST",
        f"python={sys.version.replace(chr(10), ' ')}",
        f"project_root={root}",
        f"master_input={master_path}",
        f"master_input_sha256={sha256(master_path)}",
        f"cell_richness_input={richness_path}",
        f"cell_richness_input_sha256={sha256(richness_path)}",
        f"n_phase13C_eligible_cells={richness['grid_cell_id'].nunique()}",
        f"n_primary_both_traits_present_cells={len(primary_ids)}",
        f"n_primary_complete_pairs={len(pairs)}",
        f"permutations={args.permutations}",
        f"seed={args.seed}",
        f"log_geo_mean_before_z={geo_mu}",
        f"log_geo_sd_before_z={geo_sd}",
        f"environment_mean_before_z={env_mu}",
        f"environment_sd_before_z={env_sd}",
        f"condition_number={condition_number}",
        "primary_metric=Simpson replacement",
        "secondary_metric=Jaccard total dissimilarity",
        "primary_historical_boundaries=B01_LaPaz,B03_Vizcaino",
        "primary_environment=13C domain-balanced environmental distance",
        "geographic_control=z(log1p great-circle distance km)",
        "inference=Freedman-Lane-style cell-label residual permutation",
        "ordinary_iid_regression_p_values_used=NO",
    ]
    manifest_path.write_text("\n".join(manifest) + "\n", encoding="utf-8")

    print("\nPERMUTATION TEST RESULTS")
    display_cols = [
        c for c in [
            "metric", "analysis_status", "test", "effect_family",
            "observed_beta", "observed_t", "permutation_p_two_sided",
            "observed_F", "partial_R2", "permutation_p_upper"
        ] if c in perm_df.columns
    ]
    print(perm_df[display_cols].to_string(index=False))

    print("\nHYPOTHESIS SUMMARY")
    print(hypothesis_df.to_string(index=False))

    print("\nFILES WRITTEN")
    for p in [
        cells_path, pair_path, coeff_path, perm_path, trait_path,
        vif_path, corr_path, hyp_path, manifest_path, readme_path
    ]:
        print(f"  {p}")

    print("\nSTATUS: PASS")
    print("13F completed the primary historical-vs-contemporary paired permutation test.")
    print("Interpret Simpson replacement first; Jaccard is secondary.")
    print("Next step: 13G robustness/sensitivity and publication figure synthesis.")


if __name__ == "__main__":
    main()
