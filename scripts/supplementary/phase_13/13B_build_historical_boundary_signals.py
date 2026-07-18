#!/usr/bin/env python3
"""
13B_build_historical_boundary_signals.py

Purpose
-------
Construct a locked, cell-level historical/contextual boundary signal table for
Phase 13 WITHOUT testing any arachnid C3/N0 outcome.

Primary historical hypotheses are frozen before outcome testing:
  B01 Isthmus of La Paz: 24.0-25.0 N zone
  B03 Vizcaino / mid-peninsular: 27.0-28.0 N zone

Secondary contextual transitions:
  B02 Loreto: 26.0 N line
  B04 Northern climatic transition: 30.0 N line

This script deliberately does NOT:
- calculate turnover by trait,
- compare ballooning vs non-ballooning,
- optimize boundary locations,
- move boundaries toward observed arachnid peaks.

Outputs
-------
13B_boundary_config_frozen.csv
13B_cell_boundary_signals_long.csv
13B_cell_boundary_signals_wide.csv
13B_boundary_cell_summary.csv
13B_run_manifest.txt
13B_README.txt
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

KM_PER_DEG_LAT = 111.195


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args():
    p = argparse.ArgumentParser(
        description="Build frozen Phase 13 historical/contextual boundary signals."
    )
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument(
        "--config",
        type=Path,
        default=Path.home()
        / "Downloads/PHASE_13_STARTER/configs/phase_13_boundaries_frozen.csv",
        help="Frozen boundary config CSV.",
    )
    p.add_argument(
        "--cell-lookup",
        type=Path,
        default=None,
        help="Optional override for 10_common_grid25km_cell_lookup.csv",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Optional output directory override.",
    )
    return p.parse_args()


def validate_config(cfg: pd.DataFrame) -> None:
    required = {
        "boundary_id", "name", "role", "representation", "latitude_or_path",
        "source_short", "citation", "status", "notes"
    }
    missing = required - set(cfg.columns)
    if missing:
        raise ValueError(f"Boundary config missing columns: {sorted(missing)}")

    if cfg["boundary_id"].duplicated().any():
        raise ValueError("boundary_id values must be unique.")

    allowed_roles = {"primary", "secondary"}
    bad_roles = sorted(set(cfg["role"]) - allowed_roles)
    if bad_roles:
        raise ValueError(f"Unsupported roles: {bad_roles}")

    allowed_rep = {"latitude_zone", "latitude_line"}
    bad_rep = sorted(set(cfg["representation"]) - allowed_rep)
    if bad_rep:
        raise ValueError(f"Unsupported representations: {bad_rep}")

    unresolved = cfg["status"].astype(str).str.contains("UNRESOLVED", case=False, na=False)
    if unresolved.any():
        ids = cfg.loc[unresolved, "boundary_id"].tolist()
        raise ValueError(f"Unresolved boundaries remain in config: {ids}")

    primary = cfg.loc[cfg["role"] == "primary", "boundary_id"].tolist()
    if primary != ["B01", "B03"]:
        raise ValueError(
            "Primary historical boundaries must be exactly B01 and B03 in that order. "
            f"Found: {primary}"
        )


def parse_boundary(row: pd.Series):
    rep = row["representation"]
    raw = str(row["latitude_or_path"]).strip()

    if rep == "latitude_line":
        center = float(raw)
        return center, center, center

    if rep == "latitude_zone":
        parts = raw.split(":")
        if len(parts) != 2:
            raise ValueError(
                f"{row['boundary_id']}: latitude_zone must be LOW:HIGH, got {raw}"
            )
        low, high = map(float, parts)
        if low >= high:
            raise ValueError(f"{row['boundary_id']}: zone low must be < high.")
        return low, high, (low + high) / 2.0

    raise ValueError(f"Unsupported representation {rep}")


def build_signals(cells: pd.DataFrame, cfg: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for _, b in cfg.iterrows():
        low, high, center = parse_boundary(b)

        for _, c in cells.iterrows():
            lat = float(c["centroid_latitude"])
            signed_center_deg = lat - center
            abs_center_deg = abs(signed_center_deg)

            if b["representation"] == "latitude_zone":
                inside = low <= lat <= high
                if lat < low:
                    distance_zone_deg = low - lat
                    side = "south"
                elif lat > high:
                    distance_zone_deg = lat - high
                    side = "north"
                else:
                    distance_zone_deg = 0.0
                    side = "inside_zone"
            else:
                inside = False
                distance_zone_deg = abs_center_deg
                side = "south" if lat < center else ("north" if lat > center else "on_line")

            distance_km = distance_zone_deg * KM_PER_DEG_LAT

            rows.append({
                "grid_cell_id": c["grid_cell_id"],
                "grid_cell_order": c.get("grid_cell_order", np.nan),
                "centroid_latitude": lat,
                "centroid_longitude": float(c["centroid_longitude"]),
                "centroid_latitude_band": c.get(
                    "centroid_latitude_band", c.get("latitude_band", np.nan)
                ),
                "boundary_id": b["boundary_id"],
                "boundary_name": b["name"],
                "boundary_role": b["role"],
                "boundary_representation": b["representation"],
                "boundary_low_lat": low,
                "boundary_high_lat": high,
                "boundary_center_lat": center,
                "signed_distance_to_center_deg": signed_center_deg,
                "abs_distance_to_center_deg": abs_center_deg,
                "distance_to_boundary_zone_deg": distance_zone_deg,
                "distance_to_boundary_zone_km": distance_km,
                "boundary_side": side,
                "inside_boundary_zone": bool(inside),
                "within_25km": bool(distance_km <= 25.0),
                "within_50km": bool(distance_km <= 50.0),
                "within_100km": bool(distance_km <= 100.0),
                "source_short": b["source_short"],
                "status": b["status"],
            })

    return pd.DataFrame(rows)


def make_wide(long: pd.DataFrame) -> pd.DataFrame:
    base_cols = [
        "grid_cell_id", "grid_cell_order", "centroid_latitude",
        "centroid_longitude", "centroid_latitude_band"
    ]
    base = long[base_cols].drop_duplicates("grid_cell_id").copy()

    metrics = [
        "signed_distance_to_center_deg",
        "abs_distance_to_center_deg",
        "distance_to_boundary_zone_km",
        "inside_boundary_zone",
        "within_25km",
        "within_50km",
        "within_100km",
        "boundary_side",
    ]

    for bid in long["boundary_id"].drop_duplicates():
        sub = long.loc[long["boundary_id"] == bid, ["grid_cell_id"] + metrics].copy()
        sub = sub.rename(columns={m: f"{bid}_{m}" for m in metrics})
        base = base.merge(sub, on="grid_cell_id", how="left", validate="one_to_one")

    return base.sort_values("grid_cell_order").reset_index(drop=True)


def main():
    args = parse_args()
    root = args.project_root.expanduser().resolve()

    cell_lookup = args.cell_lookup or (
        root / "02_data_clean/08_grid25km_incidence/10_common_grid25km_cell_lookup.csv"
    )
    config = args.config.expanduser().resolve()

    outdir = args.output_dir or (
        root
        / "04_analysis_USE _THIS"
        / "13_historical_vs_contemporary"
        / "13B_historical_boundary_signals"
    )
    outdir = outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print("PHASE 13B — HISTORICAL BOUNDARY SIGNALS")
    print("=" * 72)
    print(f"PROJECT ROOT : {root}")
    print(f"CELL LOOKUP  : {cell_lookup}")
    print(f"CONFIG       : {config}")
    print(f"OUTPUT DIR   : {outdir}")

    if not cell_lookup.exists():
        raise FileNotFoundError(f"Cell lookup not found: {cell_lookup}")
    if not config.exists():
        raise FileNotFoundError(f"Frozen config not found: {config}")

    cells = pd.read_csv(cell_lookup)
    cfg = pd.read_csv(config)

    required_cells = {"grid_cell_id", "centroid_latitude", "centroid_longitude"}
    missing_cells = required_cells - set(cells.columns)
    if missing_cells:
        raise ValueError(f"Cell lookup missing columns: {sorted(missing_cells)}")

    if cells["grid_cell_id"].duplicated().any():
        raise ValueError("Cell lookup contains duplicated grid_cell_id values.")

    validate_config(cfg)

    if len(cells) != 205:
        raise ValueError(
            f"Expected retained 205 occupied 25-km cells; found {len(cells)}. "
            "Stop and audit before continuing."
        )

    print("\nLOCKED BOUNDARY DESIGN")
    print(cfg[
        ["boundary_id", "name", "role", "representation", "latitude_or_path", "status"]
    ].to_string(index=False))

    long = build_signals(cells, cfg)
    wide = make_wide(long)

    summary = (
        long.groupby(
            [
                "boundary_id", "boundary_name", "boundary_role",
                "boundary_representation", "boundary_low_lat",
                "boundary_high_lat", "boundary_center_lat"
            ],
            dropna=False,
        )
        .agg(
            n_cells=("grid_cell_id", "nunique"),
            n_inside_zone=("inside_boundary_zone", "sum"),
            n_within_25km=("within_25km", "sum"),
            n_within_50km=("within_50km", "sum"),
            n_within_100km=("within_100km", "sum"),
            min_distance_km=("distance_to_boundary_zone_km", "min"),
            median_distance_km=("distance_to_boundary_zone_km", "median"),
            max_distance_km=("distance_to_boundary_zone_km", "max"),
        )
        .reset_index()
    )

    frozen_cfg_path = outdir / "13B_boundary_config_frozen.csv"
    long_path = outdir / "13B_cell_boundary_signals_long.csv"
    wide_path = outdir / "13B_cell_boundary_signals_wide.csv"
    summary_path = outdir / "13B_boundary_cell_summary.csv"
    manifest_path = outdir / "13B_run_manifest.txt"
    readme_path = outdir / "13B_README.txt"

    shutil.copy2(config, frozen_cfg_path)
    long.to_csv(long_path, index=False)
    wide.to_csv(wide_path, index=False)
    summary.to_csv(summary_path, index=False)

    readme = """PHASE 13B — FROZEN HISTORICAL BOUNDARY SIGNALS

This step creates predictor/signal geometry only. It does not inspect or test
ballooning/non-ballooning outcomes.

PRIMARY historical/vicariant hypotheses:
- B01 Isthmus of La Paz: 24-25 N a-priori zone
- B03 Vizcaino/mid-peninsular: 27-28 N a-priori zone

SECONDARY contextual transitions:
- B02 Loreto: 26 N operational line; climate-associated contact context
- B04 Northern transition: 30 N operational climatic/vegetational line

Important:
Boundary locations must not be moved after seeing arachnid outcome results.
Any alternative boundary definitions must be labeled sensitivity analyses.
"""
    readme_path.write_text(readme, encoding="utf-8")

    manifest = [
        "PHASE 13B RUN MANIFEST",
        f"python={sys.version.replace(chr(10), ' ')}",
        f"project_root={root}",
        f"cell_lookup={cell_lookup}",
        f"cell_lookup_sha256={sha256(cell_lookup)}",
        f"config_source={config}",
        f"config_source_sha256={sha256(config)}",
        f"n_cells={len(cells)}",
        f"n_boundaries={len(cfg)}",
        "primary_boundaries=B01,B03",
        "secondary_boundaries=B02,B04",
        "outcome_testing_performed=NO",
    ]
    manifest_path.write_text("\n".join(manifest) + "\n", encoding="utf-8")

    print("\nBOUNDARY CELL SUMMARY")
    print(summary.to_string(index=False))

    print("\nFILES WRITTEN")
    for p in [frozen_cfg_path, long_path, wide_path, summary_path, manifest_path, readme_path]:
        print(f"  {p}")

    print("\nSTATUS: PASS")
    print("13B built frozen boundary geometry/signals only.")
    print("No C3/N0 outcome testing was performed.")
    print("Next step: 13C contemporary-environment signal construction / locking.")


if __name__ == "__main__":
    main()
