#!/usr/bin/env python3
"""
13E_join_pairwise_historical_environment_geography_community.py

Purpose
-------
Construct the single paired analysis table for Phase 13 by joining, on exactly
the same unordered cell pairs:

1. Historical/contextual boundary separation signals from 13B
2. Contemporary environmental distance from 13C
3. Geographic distance between cell centroids
4. Paired C3 and N0 community dissimilarities from 13D

This is a DATA-CONSTRUCTION step only. It performs no inferential model fitting.

Historical pairwise coding
--------------------------
For each frozen boundary:
- strict_cross = 1 only when one cell is south of the boundary zone/line and
  the other is north. Cells inside a latitude zone do not count as a strict
  crossing.
- touches_or_crosses = 1 when the two cells occupy different boundary-side
  categories, including cases where one lies inside a frozen zone.
- same_side = 1 when both cells share the same boundary-side category.

Primary historical separation signal:
- primary_strict_cross_count = B01 + B03 strict crossings (0, 1, or 2)
- primary_any_strict_cross = 1 if either B01 or B03 is strictly crossed

Secondary/contextual signals B02 and B04 are retained separately and are not
combined into the primary historical-vicariance signal.

Geographic distance
-------------------
Great-circle (haversine) distance between 25-km cell centroids in kilometers.

Expected retained pair counts
-----------------------------
189 eligible cells -> 17,766 unordered pairs
Paired-valid C3/N0 community comparisons expected from 13D: 16,400

Outputs
-------
13E_pairwise_master_all_pairs.csv
13E_pairwise_master_paired_valid.csv
13E_pairwise_join_audit.csv
13E_pairwise_signal_summary.csv
13E_run_manifest.txt
13E_README.txt
"""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd


EXPECTED_CELLS = 189
EXPECTED_ALL_PAIRS = EXPECTED_CELLS * (EXPECTED_CELLS - 1) // 2  # 17,766
EXPECTED_PAIRED_VALID = 16400
PRIMARY_BOUNDARIES = ["B01", "B03"]
SECONDARY_BOUNDARIES = ["B02", "B04"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args():
    p = argparse.ArgumentParser(
        description="Join Phase 13 pairwise historical, environmental, geographic, and community signals."
    )
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--boundary-long", type=Path, default=None)
    p.add_argument("--environment-pairs", type=Path, default=None)
    p.add_argument("--community-pairs", type=Path, default=None)
    p.add_argument("--output-dir", type=Path, default=None)
    return p.parse_args()


def canonicalize_pairs(df: pd.DataFrame, i_col: str, j_col: str) -> pd.DataFrame:
    x = df.copy()
    a = x[i_col].astype(str)
    b = x[j_col].astype(str)
    x["_pair_a"] = np.where(a <= b, a, b)
    x["_pair_b"] = np.where(a <= b, b, a)
    x["pair_id"] = x["_pair_a"] + "||" + x["_pair_b"]

    if x["pair_id"].duplicated().any():
        dups = x.loc[x["pair_id"].duplicated(), "pair_id"].head(20).tolist()
        raise ValueError(f"Duplicate unordered pair IDs found: {dups}")

    return x


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0088
    phi1 = math.radians(float(lat1))
    phi2 = math.radians(float(lat2))
    dphi = math.radians(float(lat2) - float(lat1))
    dlambda = math.radians(float(lon2) - float(lon1))
    a = (
        math.sin(dphi / 2.0) ** 2
        + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2.0) ** 2
    )
    return 2.0 * r * math.asin(math.sqrt(a))


def build_boundary_pairs(boundary_long: pd.DataFrame, eligible_ids: set[str]) -> pd.DataFrame:
    required = {
        "grid_cell_id", "boundary_id", "boundary_name", "boundary_role",
        "boundary_side", "distance_to_boundary_zone_km"
    }
    missing = required - set(boundary_long.columns)
    if missing:
        raise ValueError(f"13B boundary long table missing columns: {sorted(missing)}")

    x = boundary_long.copy()
    x["grid_cell_id"] = x["grid_cell_id"].astype(str)
    x = x.loc[x["grid_cell_id"].isin(eligible_ids)].copy()

    observed_boundaries = set(x["boundary_id"].unique())
    expected_boundaries = set(PRIMARY_BOUNDARIES + SECONDARY_BOUNDARIES)
    if observed_boundaries != expected_boundaries:
        raise ValueError(
            f"Expected boundaries {sorted(expected_boundaries)}; "
            f"found {sorted(observed_boundaries)}"
        )

    counts = x.groupby("boundary_id")["grid_cell_id"].nunique().to_dict()
    for bid in sorted(expected_boundaries):
        if counts.get(bid, 0) != EXPECTED_CELLS:
            raise ValueError(
                f"{bid} has {counts.get(bid, 0)} eligible cells; expected {EXPECTED_CELLS}."
            )

    pair_base = None

    for bid in PRIMARY_BOUNDARIES + SECONDARY_BOUNDARIES:
        b = (
            x.loc[x["boundary_id"] == bid, [
                "grid_cell_id", "boundary_side", "distance_to_boundary_zone_km"
            ]]
            .drop_duplicates("grid_cell_id")
            .set_index("grid_cell_id")
        )

        ids = sorted(b.index.tolist())
        rows = []
        for ii in range(len(ids)):
            ci = ids[ii]
            si = str(b.loc[ci, "boundary_side"])
            di = float(b.loc[ci, "distance_to_boundary_zone_km"])

            for jj in range(ii + 1, len(ids)):
                cj = ids[jj]
                sj = str(b.loc[cj, "boundary_side"])
                dj = float(b.loc[cj, "distance_to_boundary_zone_km"])

                strict = int({si, sj} == {"south", "north"})
                different = int(si != sj)
                same = int(si == sj)

                rows.append({
                    "cell_i": ci,
                    "cell_j": cj,
                    "pair_id": ci + "||" + cj,
                    f"{bid}_strict_cross": strict,
                    f"{bid}_touches_or_crosses": different,
                    f"{bid}_same_side": same,
                    f"{bid}_min_cell_distance_to_boundary_km": min(di, dj),
                    f"{bid}_mean_cell_distance_to_boundary_km": (di + dj) / 2.0,
                })

        p = pd.DataFrame(rows)
        if len(p) != EXPECTED_ALL_PAIRS:
            raise ValueError(
                f"{bid} pair count is {len(p)}; expected {EXPECTED_ALL_PAIRS}."
            )

        if pair_base is None:
            pair_base = p
        else:
            pair_base = pair_base.merge(
                p.drop(columns=["cell_i", "cell_j"]),
                on="pair_id",
                how="inner",
                validate="one_to_one",
            )

    pair_base["primary_strict_cross_count"] = (
        pair_base["B01_strict_cross"] + pair_base["B03_strict_cross"]
    )
    pair_base["primary_any_strict_cross"] = (
        pair_base["primary_strict_cross_count"] > 0
    ).astype(int)

    pair_base["primary_touches_or_crosses_count"] = (
        pair_base["B01_touches_or_crosses"]
        + pair_base["B03_touches_or_crosses"]
    )
    pair_base["primary_any_touches_or_crosses"] = (
        pair_base["primary_touches_or_crosses_count"] > 0
    ).astype(int)

    pair_base["secondary_strict_cross_count"] = (
        pair_base["B02_strict_cross"] + pair_base["B04_strict_cross"]
    )

    return pair_base


def main():
    args = parse_args()
    root = args.project_root.expanduser().resolve()

    boundary_path = args.boundary_long or (
        root
        / "04_analysis_USE _THIS"
        / "13_historical_vs_contemporary"
        / "13B_historical_boundary_signals"
        / "13B_cell_boundary_signals_long.csv"
    )
    environment_path = args.environment_pairs or (
        root
        / "04_analysis_USE _THIS"
        / "13_historical_vs_contemporary"
        / "13C_contemporary_environment_signal"
        / "13C_pairwise_environment_distance.csv"
    )
    community_path = args.community_pairs or (
        root
        / "04_analysis_USE _THIS"
        / "13_historical_vs_contemporary"
        / "13D_paired_C3_N0_community_dissimilarity"
        / "13D_pairwise_C3_N0_dissimilarity_wide.csv"
    )
    outdir = args.output_dir or (
        root
        / "04_analysis_USE _THIS"
        / "13_historical_vs_contemporary"
        / "13E_pairwise_master"
    )

    boundary_path = boundary_path.expanduser().resolve()
    environment_path = environment_path.expanduser().resolve()
    community_path = community_path.expanduser().resolve()
    outdir = outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print("PHASE 13E — JOIN PAIRWISE HISTORICAL / ENVIRONMENT / GEOGRAPHY / COMMUNITY")
    print("=" * 86)
    print(f"PROJECT ROOT    : {root}")
    print(f"BOUNDARY SIGNAL : {boundary_path}")
    print(f"ENV PAIRS       : {environment_path}")
    print(f"COMMUNITY PAIRS : {community_path}")
    print(f"OUTPUT DIR      : {outdir}")

    for p in [boundary_path, environment_path, community_path]:
        if not p.exists():
            raise FileNotFoundError(f"Required input not found: {p}")

    boundaries = pd.read_csv(boundary_path)
    env = pd.read_csv(environment_path)
    comm = pd.read_csv(community_path)

    env = canonicalize_pairs(env, "cell_i", "cell_j")
    comm = canonicalize_pairs(comm, "cell_i", "cell_j")

    if len(env) != EXPECTED_ALL_PAIRS:
        raise ValueError(
            f"13C environment pair count = {len(env)}; expected {EXPECTED_ALL_PAIRS}."
        )
    if len(comm) != EXPECTED_ALL_PAIRS:
        raise ValueError(
            f"13D community pair count = {len(comm)}; expected {EXPECTED_ALL_PAIRS}."
        )

    env_pairs = set(env["pair_id"])
    comm_pairs = set(comm["pair_id"])
    if env_pairs != comm_pairs:
        raise ValueError(
            f"13C/13D pair sets differ: env-only={len(env_pairs - comm_pairs)}, "
            f"community-only={len(comm_pairs - env_pairs)}"
        )

    eligible_ids = set(env["_pair_a"]) | set(env["_pair_b"])
    if len(eligible_ids) != EXPECTED_CELLS:
        raise ValueError(
            f"Expected {EXPECTED_CELLS} unique eligible cells; found {len(eligible_ids)}."
        )

    hist = build_boundary_pairs(boundaries, eligible_ids)

    # Select environment columns, preserving coordinates for geography.
    required_env = {
        "pair_id", "_pair_a", "_pair_b",
        "lat_i", "lat_j", "lon_i", "lon_j",
        "envdist_thermal", "envdist_moisture", "envdist_wind",
        "envdist_vegetation", "envdist_primary_balanced",
        "envdist_all_predictors_std_euclidean"
    }
    missing_env = required_env - set(env.columns)
    if missing_env:
        raise ValueError(f"13C pair table missing columns: {sorted(missing_env)}")

    env_keep = env[[
        "pair_id", "_pair_a", "_pair_b",
        "lat_i", "lat_j", "lon_i", "lon_j",
        "envdist_thermal", "envdist_moisture", "envdist_wind",
        "envdist_vegetation", "envdist_primary_balanced",
        "envdist_all_predictors_std_euclidean"
    ]].copy()
    env_keep = env_keep.rename(columns={"_pair_a": "cell_a", "_pair_b": "cell_b"})

    # Coordinates in 13C were stored relative to original cell_i/cell_j.
    # Rebuild geography robustly from canonicalized pair endpoints.
    # Create a cell coordinate lookup from both sides of the original env table.
    coord_rows = []
    for _, r in env.iterrows():
        coord_rows.append((str(r["cell_i"]), float(r["lat_i"]), float(r["lon_i"])))
        coord_rows.append((str(r["cell_j"]), float(r["lat_j"]), float(r["lon_j"])))
    coords = pd.DataFrame(coord_rows, columns=["grid_cell_id", "lat", "lon"]).drop_duplicates()
    conflicts = coords.groupby("grid_cell_id").size()
    if (conflicts > 1).any():
        # Check whether duplicates are exact coordinate duplicates after drop_duplicates.
        bad = conflicts[conflicts > 1]
        raise ValueError(f"Conflicting coordinates found for cells: {bad.index.tolist()[:20]}")
    coords = coords.set_index("grid_cell_id")

    env_keep["geographic_distance_km"] = [
        haversine_km(
            coords.loc[a, "lat"], coords.loc[a, "lon"],
            coords.loc[b, "lat"], coords.loc[b, "lon"]
        )
        for a, b in zip(env_keep["cell_a"], env_keep["cell_b"])
    ]

    # Community columns.
    comm_keep_cols = [
        "pair_id",
        "C3_shared_a", "C3_unique_i_b", "C3_unique_j_c",
        "C3_richness_i", "C3_richness_j",
        "C3_jaccard", "C3_simpson_replacement",
        "N0_shared_a", "N0_unique_i_b", "N0_unique_j_c",
        "N0_richness_i", "N0_richness_j",
        "N0_jaccard", "N0_simpson_replacement",
        "delta_jaccard_C3_minus_N0",
        "delta_simpson_C3_minus_N0",
        "paired_valid_jaccard", "paired_valid_simpson",
    ]
    missing_comm = [c for c in comm_keep_cols if c not in comm.columns]
    if missing_comm:
        raise ValueError(f"13D pair table missing columns: {missing_comm}")

    comm_keep = comm[comm_keep_cols].copy()

    master = (
        env_keep
        .merge(hist.drop(columns=["cell_i", "cell_j"]), on="pair_id", how="inner", validate="one_to_one")
        .merge(comm_keep, on="pair_id", how="inner", validate="one_to_one")
    )

    if len(master) != EXPECTED_ALL_PAIRS:
        raise ValueError(
            f"Joined master has {len(master)} pairs; expected {EXPECTED_ALL_PAIRS}."
        )

    paired_valid = master.loc[
        master["paired_valid_jaccard"].astype(bool)
        & master["paired_valid_simpson"].astype(bool)
    ].copy()

    if len(paired_valid) != EXPECTED_PAIRED_VALID:
        raise ValueError(
            f"Paired-valid master has {len(paired_valid)} pairs; "
            f"expected {EXPECTED_PAIRED_VALID} from 13D."
        )

    # Join audit.
    audit_rows = [
        {"check": "unique_eligible_cells", "observed": len(eligible_ids), "expected": EXPECTED_CELLS, "status": "PASS"},
        {"check": "environment_pairs", "observed": len(env), "expected": EXPECTED_ALL_PAIRS, "status": "PASS"},
        {"check": "community_pairs", "observed": len(comm), "expected": EXPECTED_ALL_PAIRS, "status": "PASS"},
        {"check": "historical_pairs", "observed": len(hist), "expected": EXPECTED_ALL_PAIRS, "status": "PASS"},
        {"check": "joined_master_pairs", "observed": len(master), "expected": EXPECTED_ALL_PAIRS, "status": "PASS"},
        {"check": "paired_valid_pairs", "observed": len(paired_valid), "expected": EXPECTED_PAIRED_VALID, "status": "PASS"},
    ]
    audit = pd.DataFrame(audit_rows)

    # Signal summaries.
    summary_vars = [
        "geographic_distance_km",
        "envdist_primary_balanced",
        "envdist_thermal",
        "envdist_moisture",
        "envdist_wind",
        "envdist_vegetation",
        "primary_strict_cross_count",
        "primary_any_strict_cross",
        "B01_strict_cross",
        "B03_strict_cross",
        "B02_strict_cross",
        "B04_strict_cross",
        "C3_jaccard",
        "N0_jaccard",
        "C3_simpson_replacement",
        "N0_simpson_replacement",
        "delta_jaccard_C3_minus_N0",
        "delta_simpson_C3_minus_N0",
    ]

    summary_rows = []
    for dataset_name, d in [("all_pairs", master), ("paired_valid", paired_valid)]:
        for col in summary_vars:
            s = pd.to_numeric(d[col], errors="coerce").dropna()
            summary_rows.append({
                "dataset": dataset_name,
                "variable": col,
                "n": len(s),
                "mean": float(s.mean()),
                "sd": float(s.std(ddof=1)) if len(s) > 1 else np.nan,
                "median": float(s.median()),
                "q025": float(s.quantile(0.025)),
                "q975": float(s.quantile(0.975)),
                "min": float(s.min()),
                "max": float(s.max()),
            })
    summary = pd.DataFrame(summary_rows)

    all_path = outdir / "13E_pairwise_master_all_pairs.csv"
    valid_path = outdir / "13E_pairwise_master_paired_valid.csv"
    audit_path = outdir / "13E_pairwise_join_audit.csv"
    summary_path = outdir / "13E_pairwise_signal_summary.csv"
    manifest_path = outdir / "13E_run_manifest.txt"
    readme_path = outdir / "13E_README.txt"

    master.to_csv(all_path, index=False)
    paired_valid.to_csv(valid_path, index=False)
    audit.to_csv(audit_path, index=False)
    summary.to_csv(summary_path, index=False)

    readme = """PHASE 13E — PAIRED MASTER TABLE

This step joins four pre-locked information layers on exactly the same unordered
25-km cell pairs:

1. Historical boundary separation:
   Primary = B01 Isthmus of La Paz + B03 Vizcaino/mid-peninsular.
   Secondary = B02 Loreto + B04 northern climatic transition.
2. Contemporary environmental distance from 13C.
3. Great-circle geographic distance between cell centroids.
4. Paired C3 and N0 community dissimilarities from 13D.

Primary historical coding:
- strict_cross requires one cell south and one north of the frozen boundary.
- For broad primary zones, a cell inside the zone does NOT create a strict
  crossing. A separate touches_or_crosses variable retains those cases.

No statistical inference is performed in 13E.
No boundary location, environmental predictor, or pair is selected based on
C3/N0 outcome magnitude.

The paired-valid table is the primary input for 13F and contains only pairs for
which both C3 and N0 Jaccard/Simpson values are defined.
"""

    readme_path.write_text(readme, encoding="utf-8")

    manifest = [
        "PHASE 13E RUN MANIFEST",
        f"python={sys.version.replace(chr(10), ' ')}",
        f"project_root={root}",
        f"boundary_input={boundary_path}",
        f"boundary_input_sha256={sha256(boundary_path)}",
        f"environment_input={environment_path}",
        f"environment_input_sha256={sha256(environment_path)}",
        f"community_input={community_path}",
        f"community_input_sha256={sha256(community_path)}",
        f"n_cells={EXPECTED_CELLS}",
        f"n_all_pairs={len(master)}",
        f"n_paired_valid_pairs={len(paired_valid)}",
        "primary_historical_boundaries=B01,B03",
        "secondary_contextual_boundaries=B02,B04",
        "primary_environment_signal=envdist_primary_balanced",
        "geographic_control=great-circle centroid distance km",
        "inference_performed=NO",
    ]
    manifest_path.write_text("\n".join(manifest) + "\n", encoding="utf-8")

    print("\nJOIN AUDIT")
    print(audit.to_string(index=False))

    print("\nPRIMARY SIGNAL COUNTS — PAIRED VALID")
    print(
        paired_valid[
            ["primary_strict_cross_count", "primary_any_strict_cross",
             "B01_strict_cross", "B03_strict_cross"]
        ]
        .value_counts()
        .head(20)
        .to_string()
    )

    print("\nFILES WRITTEN")
    for p in [all_path, valid_path, audit_path, summary_path, manifest_path, readme_path]:
        print(f"  {p}")

    print("\nSTATUS: PASS")
    print("13E joined historical, contemporary, geographic, and paired community signals.")
    print(f"All pairs: {len(master)}")
    print(f"Paired-valid C3/N0 pairs: {len(paired_valid)}")
    print("No inferential testing was performed.")
    print("Next step: 13F test whether C3 vs N0 differ in historical-boundary versus contemporary-environment association.")


if __name__ == "__main__":
    main()
