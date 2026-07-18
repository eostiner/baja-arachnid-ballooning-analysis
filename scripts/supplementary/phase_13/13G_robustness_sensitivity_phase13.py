#!/usr/bin/env python3
"""
13G_robustness_sensitivity_phase13.py

Purpose
-------
Pre-specified robustness and sensitivity analysis for Phase 13F.

This step does NOT search for a significant formulation. It asks whether the
13F conclusion changes materially under reasonable, clearly labeled
alternatives.

Primary 13F result being stress-tested:
- Response: C3 - N0 Simpson replacement (primary)
- Secondary response: C3 - N0 Jaccard total dissimilarity
- Predictors: geography + contemporary environment + frozen B01/B03 boundaries
- Inference: cell-label residual permutation on complete matched cell matrices

Sensitivity families
--------------------
A. Cell-set thresholds (complete matrices only)
   1. C3>=1 and N0>=1 genus per cell  [13F primary]
   2. C3>=2 and N0>=2
   3. C3>=3 and N0>=3

B. Historical boundary coding
   4. B01/B03 strict crossing          [13F primary]
   5. B01/B03 touches-or-crosses       [alternative zone coding]

C. Contemporary environment representation
   6. balanced composite               [13F primary]
   7. thermal domain only
   8. moisture domain only
   9. wind domain only
  10. vegetation domain only

D. Secondary/contextual transitions (exploratory only)
  11. B02 Loreto + B04 northern transition

E. All paired-valid pairs
  12. All 16,400 paired-valid pairs are fit descriptively only.
      No permutation p-value is attached because the pair set is incomplete
      and ordinary row-wise inference would be anti-conservative.

Primary sensitivity metric:
- Simpson replacement

Jaccard is repeated as a secondary metric.

Inference
---------
For complete-matrix sensitivities, the script uses the same Freedman-Lane-style
cell-label residual permutation logic as 13F.

Outputs
-------
13G_sensitivity_models.csv
13G_sensitivity_permutation_tests.csv
13G_all_paired_valid_descriptive.csv
13G_robustness_summary.csv
13G_primary_conclusion_check.csv
13G_sensitivity_effects.png
13G_sensitivity_effects.svg
13G_run_manifest.txt
13G_README.txt
"""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
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
    p = argparse.ArgumentParser(description="Phase 13G robustness/sensitivity.")
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--master", type=Path, default=None)
    p.add_argument("--cell-richness", type=Path, default=None)
    p.add_argument("--output-dir", type=Path, default=None)
    p.add_argument("--permutations", type=int, default=1999)
    p.add_argument("--seed", type=int, default=20260718)
    return p.parse_args()


def zscore(x):
    x = np.asarray(x, dtype=float)
    mu = float(np.mean(x))
    sd = float(np.std(x, ddof=1))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("Cannot z-standardize predictor.")
    return (x - mu) / sd


def prepare_ols(X):
    X = np.asarray(X, dtype=float)
    inv = np.linalg.pinv(X.T @ X)
    beta_map = inv @ X.T
    rank = np.linalg.matrix_rank(X)
    df_resid = X.shape[0] - rank
    if df_resid <= 0:
        raise ValueError("Nonpositive residual df.")
    return {"X": X, "inv": inv, "beta_map": beta_map, "df_resid": df_resid}


def fit(prep, y):
    y = np.asarray(y, dtype=float)
    X = prep["X"]
    beta = prep["beta_map"] @ y
    fitted = X @ beta
    resid = y - fitted
    rss = float(resid @ resid)
    tss = float(((y - y.mean()) ** 2).sum())
    r2 = np.nan if tss <= 0 else 1 - rss / tss
    s2 = rss / prep["df_resid"]
    se = np.sqrt(np.maximum(np.diag(prep["inv"]) * s2, 0.0))
    t = np.divide(beta, se, out=np.full_like(beta, np.nan), where=se > 0)
    return {"beta": beta, "fitted": fitted, "resid": resid, "rss": rss,
            "r2": r2, "se": se, "t": t, "df_resid": prep["df_resid"]}


def make_matrix(values, ai, bi, n):
    m = np.zeros((n, n), dtype=float)
    m[ai, bi] = values
    m[bi, ai] = values
    np.fill_diagonal(m, 0.0)
    return m


def f_compare(fr, ff, q):
    if ff["rss"] <= 0:
        return np.nan
    num = max(fr["rss"] - ff["rss"], 0.0) / q
    den = ff["rss"] / ff["df_resid"]
    return np.nan if den <= 0 else num / den


def partial_r2(fr, ff):
    return np.nan if fr["rss"] <= 0 else max(0.0, (fr["rss"] - ff["rss"]) / fr["rss"])


def perm_term(y, Xfull, names, reduced_names, term, ai, bi, ncell, nperm, seed):
    names = list(names)
    keep = [names.index(x) for x in reduced_names]
    tidx = names.index(term)

    pf = prepare_ols(Xfull)
    ff = fit(pf, y)
    pr = prepare_ols(Xfull[:, keep])
    fr = fit(pr, y)

    obs_t = float(ff["t"][tidx])
    obs_beta = float(ff["beta"][tidx])

    rmat = make_matrix(fr["resid"], ai, bi, ncell)
    rng = np.random.default_rng(seed)
    extreme = 0

    for _ in range(nperm):
        perm = rng.permutation(ncell)
        rp = rmat[np.ix_(perm, perm)][ai, bi]
        ys = fr["fitted"] + rp
        fs = fit(pf, ys)
        ts = float(fs["t"][tidx])
        if np.isfinite(ts) and abs(ts) >= abs(obs_t) - 1e-12:
            extreme += 1

    return {
        "term": term,
        "beta": obs_beta,
        "t": obs_t,
        "p_perm_two_sided": (extreme + 1) / (nperm + 1),
    }


def perm_joint(y, Xfull, names, reduced_names, group_name, q, ai, bi, ncell, nperm, seed):
    names = list(names)
    keep = [names.index(x) for x in reduced_names]

    pf = prepare_ols(Xfull)
    ff = fit(pf, y)
    pr = prepare_ols(Xfull[:, keep])
    fr = fit(pr, y)

    obs_f = f_compare(fr, ff, q)
    pr2 = partial_r2(fr, ff)
    rmat = make_matrix(fr["resid"], ai, bi, ncell)
    rng = np.random.default_rng(seed)
    extreme = 0

    for _ in range(nperm):
        perm = rng.permutation(ncell)
        rp = rmat[np.ix_(perm, perm)][ai, bi]
        ys = fr["fitted"] + rp
        ffs = fit(pf, ys)
        frs = fit(pr, ys)
        fs = f_compare(frs, ffs, q)
        if np.isfinite(fs) and fs >= obs_f - 1e-12:
            extreme += 1

    return {
        "term": group_name,
        "F": obs_f,
        "partial_R2": pr2,
        "p_perm_upper": (extreme + 1) / (nperm + 1),
    }


def complete_pair_subset(master, richness, threshold):
    ids = richness.loc[
        (richness["C3_richness"] >= threshold) &
        (richness["N0_richness"] >= threshold),
        "grid_cell_id"
    ].astype(str).tolist()
    ids = sorted(ids)
    s = set(ids)

    p = master.loc[
        master["cell_a"].astype(str).isin(s) &
        master["cell_b"].astype(str).isin(s)
    ].copy()

    expected = len(ids) * (len(ids) - 1) // 2
    if len(p) != expected:
        raise ValueError(
            f"Threshold {threshold}: found {len(p)} pairs, expected {expected} "
            f"for complete matrix of {len(ids)} cells."
        )
    return ids, p


def run_complete_spec(
    spec_name, spec_family, analysis_status, metric_name, response_col,
    cell_ids, pairs, env_col, b1_col, b2_col, b1_label, b2_label,
    nperm, seed
):
    if len(cell_ids) < 30:
        return None, None

    cell_to_idx = {c: i for i, c in enumerate(cell_ids)}
    ai = pairs["cell_a"].astype(str).map(cell_to_idx).to_numpy(dtype=int)
    bi = pairs["cell_b"].astype(str).map(cell_to_idx).to_numpy(dtype=int)

    y = pairs[response_col].to_numpy(dtype=float)
    if not np.isfinite(y).all():
        raise ValueError(f"{spec_name}: nonfinite response in complete-matrix spec.")

    zgeo = zscore(np.log1p(pairs["geographic_distance_km"].to_numpy(dtype=float)))
    zenv = zscore(pairs[env_col].to_numpy(dtype=float))
    b1 = pairs[b1_col].to_numpy(dtype=float)
    b2 = pairs[b2_col].to_numpy(dtype=float)

    names = ["Intercept", "z_geo", "z_env", b1_label, b2_label]
    X = np.column_stack([np.ones(len(pairs)), zgeo, zenv, b1, b2])

    pf = prepare_ols(X)
    ff = fit(pf, y)
    beta_map = dict(zip(names, ff["beta"]))

    model_row = {
        "spec_name": spec_name,
        "spec_family": spec_family,
        "analysis_status": analysis_status,
        "metric": metric_name,
        "n_cells": len(cell_ids),
        "n_pairs": len(pairs),
        "environment_variable": env_col,
        "boundary_1_variable": b1_col,
        "boundary_2_variable": b2_col,
        "model_R2": float(ff["r2"]),
        "beta_environment": float(beta_map["z_env"]),
        "beta_boundary_1": float(beta_map[b1_label]),
        "beta_boundary_2": float(beta_map[b2_label]),
    }

    env_test = perm_term(
        y, X, names,
        ["Intercept", "z_geo", b1_label, b2_label],
        "z_env", ai, bi, len(cell_ids), nperm, seed + 1
    )
    hist_test = perm_joint(
        y, X, names,
        ["Intercept", "z_geo", "z_env"],
        "historical_joint", 2, ai, bi, len(cell_ids), nperm, seed + 2
    )
    b1_test = perm_term(
        y, X, names,
        ["Intercept", "z_geo", "z_env", b2_label],
        b1_label, ai, bi, len(cell_ids), nperm, seed + 3
    )
    b2_test = perm_term(
        y, X, names,
        ["Intercept", "z_geo", "z_env", b1_label],
        b2_label, ai, bi, len(cell_ids), nperm, seed + 4
    )

    tests = []
    for kind, obj in [
        ("environment", env_test),
        ("historical_joint", hist_test),
        ("boundary_1", b1_test),
        ("boundary_2", b2_test),
    ]:
        r = {
            "spec_name": spec_name,
            "spec_family": spec_family,
            "analysis_status": analysis_status,
            "metric": metric_name,
            "test_family": kind,
            "n_cells": len(cell_ids),
            "n_pairs": len(pairs),
            "n_permutations": nperm,
        }
        r.update(obj)
        tests.append(r)

    return model_row, tests


def main():
    args = parse_args()
    if args.permutations < 199:
        raise ValueError("Use at least 199 permutations.")

    root = args.project_root.expanduser().resolve()
    master_path = args.master or (
        root / "04_analysis_USE _THIS" / "13_historical_vs_contemporary"
        / "13E_pairwise_master" / "13E_pairwise_master_all_pairs.csv"
    )
    richness_path = args.cell_richness or (
        root / "04_analysis_USE _THIS" / "13_historical_vs_contemporary"
        / "13D_paired_C3_N0_community_dissimilarity" / "13D_cell_trait_richness.csv"
    )
    outdir = args.output_dir or (
        root / "04_analysis_USE _THIS" / "13_historical_vs_contemporary"
        / "13G_robustness_sensitivity"
    )

    master_path = master_path.expanduser().resolve()
    richness_path = richness_path.expanduser().resolve()
    outdir = outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print("PHASE 13G — ROBUSTNESS / SENSITIVITY")
    print("=" * 76)
    print(f"PROJECT ROOT : {root}")
    print(f"MASTER PAIRS : {master_path}")
    print(f"CELL RICHNESS: {richness_path}")
    print(f"OUTPUT DIR   : {outdir}")
    print(f"PERMUTATIONS : {args.permutations}")
    print(f"SEED         : {args.seed}")

    master = pd.read_csv(master_path)
    richness = pd.read_csv(richness_path)

    master["cell_a"] = master["cell_a"].astype(str)
    master["cell_b"] = master["cell_b"].astype(str)
    richness["grid_cell_id"] = richness["grid_cell_id"].astype(str)

    # Recreate paired deltas if needed.
    master["delta_simpson_C3_minus_N0"] = (
        master["C3_simpson_replacement"] - master["N0_simpson_replacement"]
    )
    master["delta_jaccard_C3_minus_N0"] = (
        master["C3_jaccard"] - master["N0_jaccard"]
    )

    # All-paired-valid descriptive comparison only.
    all_valid = master.loc[
        master["paired_valid_jaccard"].astype(bool) &
        master["paired_valid_simpson"].astype(bool)
    ].copy()

    descriptive_rows = []
    for metric, response in [
        ("simpson_replacement", "delta_simpson_C3_minus_N0"),
        ("jaccard_total", "delta_jaccard_C3_minus_N0"),
    ]:
        d = all_valid.dropna(subset=[
            response, "geographic_distance_km", "envdist_primary_balanced",
            "B01_strict_cross", "B03_strict_cross"
        ]).copy()

        X = np.column_stack([
            np.ones(len(d)),
            zscore(np.log1p(d["geographic_distance_km"].to_numpy(dtype=float))),
            zscore(d["envdist_primary_balanced"].to_numpy(dtype=float)),
            d["B01_strict_cross"].to_numpy(dtype=float),
            d["B03_strict_cross"].to_numpy(dtype=float),
        ])
        names = ["Intercept", "z_geo", "z_env", "B01_LaPaz", "B03_Vizcaino"]
        f = fit(prepare_ols(X), d[response].to_numpy(dtype=float))

        for name, beta in zip(names, f["beta"]):
            descriptive_rows.append({
                "dataset": "all_paired_valid_incomplete_pair_set",
                "metric": metric,
                "n_pairs": len(d),
                "term": name,
                "beta_descriptive_only": float(beta),
                "model_R2_descriptive_only": float(f["r2"]),
                "permutation_inference": "NOT_PERFORMED",
                "reason": "Incomplete pair graph; ordinary row-wise p-values would be anti-conservative.",
            })

    # Complete-matrix specs.
    model_rows = []
    test_rows = []

    metric_specs = [
        ("simpson_replacement", "PRIMARY", "delta_simpson_C3_minus_N0"),
        ("jaccard_total", "SECONDARY", "delta_jaccard_C3_minus_N0"),
    ]

    # Spec definitions generated for threshold 1 unless threshold-specific.
    env_specs = [
        ("balanced_environment", "envdist_primary_balanced", "PRIMARY"),
        ("thermal_only", "envdist_thermal", "SENSITIVITY"),
        ("moisture_only", "envdist_moisture", "SENSITIVITY"),
        ("wind_only", "envdist_wind", "SENSITIVITY"),
        ("vegetation_only", "envdist_vegetation", "SENSITIVITY"),
    ]

    spec_counter = 0

    for metric_name, metric_status, response_col in metric_specs:
        # A. thresholds 1,2,3 with primary model.
        for threshold in [1, 2, 3]:
            ids, pairs = complete_pair_subset(master, richness, threshold)
            spec_counter += 1
            row, tests = run_complete_spec(
                spec_name=f"threshold_ge_{threshold}_primary",
                spec_family="cell_threshold",
                analysis_status=("PRIMARY_REPLICATION" if threshold == 1 else "SENSITIVITY"),
                metric_name=metric_name,
                response_col=response_col,
                cell_ids=ids,
                pairs=pairs,
                env_col="envdist_primary_balanced",
                b1_col="B01_strict_cross",
                b2_col="B03_strict_cross",
                b1_label="B01_LaPaz",
                b2_label="B03_Vizcaino",
                nperm=args.permutations,
                seed=args.seed + spec_counter * 1000,
            )
            if row is not None:
                model_rows.append(row)
                test_rows.extend(tests)

        # Threshold 1 base for coding/domain/context sensitivities.
        ids1, pairs1 = complete_pair_subset(master, richness, 1)

        # B. alternate boundary coding.
        spec_counter += 1
        row, tests = run_complete_spec(
            spec_name="touches_or_crosses_B01_B03",
            spec_family="boundary_coding",
            analysis_status="SENSITIVITY",
            metric_name=metric_name,
            response_col=response_col,
            cell_ids=ids1,
            pairs=pairs1,
            env_col="envdist_primary_balanced",
            b1_col="B01_touches_or_crosses",
            b2_col="B03_touches_or_crosses",
            b1_label="B01_LaPaz",
            b2_label="B03_Vizcaino",
            nperm=args.permutations,
            seed=args.seed + spec_counter * 1000,
        )
        model_rows.append(row)
        test_rows.extend(tests)

        # C. environmental domain sensitivities.
        for env_name, env_col, env_status in env_specs[1:]:
            spec_counter += 1
            row, tests = run_complete_spec(
                spec_name=env_name,
                spec_family="environment_domain",
                analysis_status="SENSITIVITY",
                metric_name=metric_name,
                response_col=response_col,
                cell_ids=ids1,
                pairs=pairs1,
                env_col=env_col,
                b1_col="B01_strict_cross",
                b2_col="B03_strict_cross",
                b1_label="B01_LaPaz",
                b2_label="B03_Vizcaino",
                nperm=args.permutations,
                seed=args.seed + spec_counter * 1000,
            )
            model_rows.append(row)
            test_rows.extend(tests)

        # D. secondary/contextual transitions.
        spec_counter += 1
        row, tests = run_complete_spec(
            spec_name="secondary_B02_Loreto_B04_North",
            spec_family="secondary_contextual_boundaries",
            analysis_status="EXPLORATORY",
            metric_name=metric_name,
            response_col=response_col,
            cell_ids=ids1,
            pairs=pairs1,
            env_col="envdist_primary_balanced",
            b1_col="B02_strict_cross",
            b2_col="B04_strict_cross",
            b1_label="B02_Loreto",
            b2_label="B04_North",
            nperm=args.permutations,
            seed=args.seed + spec_counter * 1000,
        )
        model_rows.append(row)
        test_rows.extend(tests)

    models = pd.DataFrame(model_rows)
    tests = pd.DataFrame(test_rows)
    descriptive = pd.DataFrame(descriptive_rows)

    # Robustness summary for primary Simpson metric.
    sim_models = models.loc[models["metric"] == "simpson_replacement"].copy()
    sim_env_tests = tests.loc[
        (tests["metric"] == "simpson_replacement") &
        (tests["test_family"] == "environment")
    ].copy()
    sim_hist_tests = tests.loc[
        (tests["metric"] == "simpson_replacement") &
        (tests["test_family"] == "historical_joint")
    ].copy()

    summary_rows = []
    for _, m in sim_models.iterrows():
        e = sim_env_tests.loc[sim_env_tests["spec_name"] == m["spec_name"]].iloc[0]
        h = sim_hist_tests.loc[sim_hist_tests["spec_name"] == m["spec_name"]].iloc[0]
        summary_rows.append({
            "spec_name": m["spec_name"],
            "spec_family": m["spec_family"],
            "analysis_status": m["analysis_status"],
            "n_cells": m["n_cells"],
            "n_pairs": m["n_pairs"],
            "beta_environment": m["beta_environment"],
            "environment_direction_expected_positive": bool(m["beta_environment"] > 0),
            "environment_perm_p": e.get("p_perm_two_sided", np.nan),
            "beta_boundary_1": m["beta_boundary_1"],
            "beta_boundary_2": m["beta_boundary_2"],
            "any_boundary_direction_expected_negative": bool(
                (m["beta_boundary_1"] < 0) or (m["beta_boundary_2"] < 0)
            ),
            "historical_joint_partial_R2": h.get("partial_R2", np.nan),
            "historical_joint_perm_p": h.get("p_perm_upper", np.nan),
        })

    robustness = pd.DataFrame(summary_rows)

    # Conclusion check excludes exploratory secondary-boundary spec.
    confirmatory = robustness.loc[
        robustness["analysis_status"].isin(["PRIMARY_REPLICATION", "SENSITIVITY"])
    ].copy()

    n_env_sig_expected = int(
        ((confirmatory["beta_environment"] > 0) &
         (confirmatory["environment_perm_p"] < 0.05)).sum()
    )
    n_hist_sig = int((confirmatory["historical_joint_perm_p"] < 0.05).sum())

    primary_conclusion = pd.DataFrame([{
        "primary_metric": "simpson_replacement",
        "n_confirmatory_sensitivity_specs": len(confirmatory),
        "n_specs_supporting_significant_positive_environment_difference": n_env_sig_expected,
        "n_specs_with_significant_joint_historical_difference": n_hist_sig,
        "robustness_conclusion": (
            "13F_null_conclusion_robust"
            if n_env_sig_expected == 0 and n_hist_sig == 0
            else "some_sensitivity_specs_change_inference_review_individually"
        ),
        "interpretive_rule": (
            "Do not overturn 13F based on isolated exploratory results; "
            "review pre-specified confirmatory sensitivities and effect directions."
        )
    }])

    # Figure synthesis.
    png_path = outdir / "13G_sensitivity_effects.png"
    svg_path = outdir / "13G_sensitivity_effects.svg"
    try:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(11, 6))
        plot_df = robustness.loc[
            robustness["analysis_status"] != "EXPLORATORY"
        ].copy().reset_index(drop=True)

        x = np.arange(len(plot_df))
        ax.axhline(0, linewidth=1)
        ax.plot(x, plot_df["beta_environment"], marker="o", label="Environment β (C3−N0)")
        ax.plot(x, plot_df["beta_boundary_1"], marker="s", label="Boundary 1 β (C3−N0)")
        ax.plot(x, plot_df["beta_boundary_2"], marker="^", label="Boundary 2 β (C3−N0)")
        ax.set_xticks(x)
        ax.set_xticklabels(plot_df["spec_name"], rotation=45, ha="right")
        ax.set_ylabel("Paired-difference coefficient")
        ax.set_title("Phase 13G robustness: Simpson replacement effect directions")
        ax.legend()
        fig.tight_layout()
        fig.savefig(png_path, dpi=300, bbox_inches="tight")
        fig.savefig(svg_path, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        print(f"\nWARNING: Figure synthesis skipped: {e}")

    # Write outputs.
    models_path = outdir / "13G_sensitivity_models.csv"
    tests_path = outdir / "13G_sensitivity_permutation_tests.csv"
    desc_path = outdir / "13G_all_paired_valid_descriptive.csv"
    robust_path = outdir / "13G_robustness_summary.csv"
    conclusion_path = outdir / "13G_primary_conclusion_check.csv"
    manifest_path = outdir / "13G_run_manifest.txt"
    readme_path = outdir / "13G_README.txt"

    models.to_csv(models_path, index=False)
    tests.to_csv(tests_path, index=False)
    descriptive.to_csv(desc_path, index=False)
    robustness.to_csv(robust_path, index=False)
    primary_conclusion.to_csv(conclusion_path, index=False)

    readme = f"""PHASE 13G — ROBUSTNESS / SENSITIVITY

This step was specified to test whether the Phase 13F conclusion is stable,
not to search for a significant alternative.

Permutation count per complete-matrix test: {args.permutations}

CONFIRMATORY SENSITIVITIES
- Cell thresholds >=1, >=2, >=3 genera in both C3 and N0
- Strict boundary crossing vs touches-or-crosses
- Thermal, moisture, wind, vegetation domain distances separately

EXPLORATORY ONLY
- B02 Loreto + B04 northern transition

ALL 16,400 PAIRED-VALID PAIRS
These are fit descriptively only. No ordinary row-wise p-values are reported
because the incomplete dyadic graph does not satisfy iid assumptions and does
not support the same simple complete-matrix cell-label permutation scheme.

PRIMARY INTERPRETATION RULE
The Phase 13F conclusion should only be considered materially changed if a
pre-specified confirmatory sensitivity shows a coherent directional and
permutation-supported reversal, not because one exploratory model happens to
cross P < 0.05.
"""
    readme_path.write_text(readme, encoding="utf-8")

    manifest = [
        "PHASE 13G RUN MANIFEST",
        f"python={sys.version.replace(chr(10), ' ')}",
        f"project_root={root}",
        f"master_input={master_path}",
        f"master_input_sha256={sha256(master_path)}",
        f"richness_input={richness_path}",
        f"richness_input_sha256={sha256(richness_path)}",
        f"permutations_per_test={args.permutations}",
        f"seed={args.seed}",
        "purpose=robustness not significance search",
        "primary_metric=Simpson replacement",
        "secondary_metric=Jaccard total dissimilarity",
        "all_paired_valid_inference=descriptive_only",
    ]
    manifest_path.write_text("\n".join(manifest) + "\n", encoding="utf-8")

    print("\nPRIMARY SIMPSON ROBUSTNESS SUMMARY")
    print(robustness.to_string(index=False))

    print("\nPRIMARY CONCLUSION CHECK")
    print(primary_conclusion.to_string(index=False))

    print("\nFILES WRITTEN")
    for p in [
        models_path, tests_path, desc_path, robust_path, conclusion_path,
        png_path, svg_path, manifest_path, readme_path
    ]:
        if p.exists():
            print(f"  {p}")

    print("\nSTATUS: PASS")
    print("13G robustness/sensitivity analysis completed.")
    print("Review 13G_primary_conclusion_check.csv first.")
    print("Next step: lock final Phase 13 interpretation and decide main-text vs supplement.")


if __name__ == "__main__":
    main()
