#!/usr/bin/env python3
"""Phase 14C0 — freeze exact 25-km polygons for cells eligible for temporal H3 analysis."""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from phase14_common import (
    cell_polygon_lonlat,
    default_analysis_output_root,
    read_delimited,
    write_csv,
    write_json,
)

SCRIPT_VERSION = "14C0_v0.2.0_2026-07-18"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare exact eligible 25-km cell polygons for Earth Engine extraction.")
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--audit-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    return parser.parse_args()


def mean_or_blank(values: list[float]) -> float | str:
    return sum(values) / len(values) if values else ""


def main() -> int:
    args = parse_args()
    root = args.project_root.expanduser().resolve()
    base = default_analysis_output_root(root)
    audit_dir = args.audit_dir.expanduser().resolve() if args.audit_dir else base / "14A_temporal_feasibility_audit"
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else base / "14C_real_stress_anomalies"
    outdir.mkdir(parents=True, exist_ok=True)

    pairs_path = audit_dir / "14A_adjacent_period_pair_eligibility.csv"
    coverage_path = audit_dir / "14A_cell_period_coverage.csv"
    if not pairs_path.exists() or not coverage_path.exists():
        raise FileNotFoundError("Phase 14A eligibility outputs not found. Run 14A before 14C0.")

    _, pair_rows = read_delimited(pairs_path)
    _, coverage_rows = read_delimited(coverage_path)
    eligible_pairs = [row for row in pair_rows if int(row.get("temporal_pair_eligible", 0)) == 1]
    eligible_cells = sorted({row["grid_cell_id"] for row in eligible_pairs})
    if not eligible_cells:
        raise RuntimeError("No eligible temporal cells were found in Phase 14A outputs.")

    coords: dict[str, dict[str, list[float]]] = defaultdict(lambda: {"lat": [], "lon": []})
    bands: dict[str, Counter[str]] = defaultdict(Counter)
    for row in coverage_rows:
        cell = row.get("grid_cell_id", "")
        if cell not in eligible_cells:
            continue
        try:
            coords[cell]["lat"].append(float(row["centroid_latitude_mean"]))
            coords[cell]["lon"].append(float(row["centroid_longitude_mean"]))
        except (TypeError, ValueError, KeyError):
            pass
        if row.get("latitude_band"):
            bands[cell][row["latitude_band"]] += 1

    pair_counts = Counter(row["grid_cell_id"] for row in eligible_pairs)
    transitions: dict[str, list[str]] = defaultdict(list)
    for row in eligible_pairs:
        transitions[row["grid_cell_id"]].append(f"{row['period_1']}->{row['period_2']}")

    csv_rows: list[dict[str, Any]] = []
    features: list[dict[str, Any]] = []
    for cell in eligible_cells:
        centroid_lat = mean_or_blank(coords[cell]["lat"])
        centroid_lon = mean_or_blank(coords[cell]["lon"])
        band = bands[cell].most_common(1)[0][0] if bands[cell] else ""
        ring = cell_polygon_lonlat(cell)
        row = {
            "grid_cell_id": cell,
            "centroid_latitude_mean": centroid_lat,
            "centroid_longitude_mean": centroid_lon,
            "latitude_band": band,
            "eligible_temporal_pair_count": pair_counts[cell],
            "eligible_transitions": ";".join(sorted(transitions[cell])),
        }
        csv_rows.append(row)
        features.append(
            {
                "type": "Feature",
                "properties": row,
                "geometry": {"type": "Polygon", "coordinates": [ring]},
            }
        )

    csv_path = outdir / "14C0_eligible_temporal_cells.csv"
    geojson_path = outdir / "14C0_eligible_temporal_cells.geojson"
    write_csv(csv_path, csv_rows)
    geojson_path.write_text(
        json.dumps({"type": "FeatureCollection", "features": features}, indent=2) + "\n",
        encoding="utf-8",
    )
    status = {
        "phase": "14C0",
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETED",
        "eligible_cells": len(eligible_cells),
        "eligible_temporal_pairs": len(eligible_pairs),
        "csv": str(csv_path),
        "geojson": str(geojson_path),
        "grid_geometry": "Exact 25-km square reconstructed from retained spherical LAEA row/column identifier.",
    }
    write_json(outdir / "14C0_run_status.json", status)
    print("PHASE 14C0 — ELIGIBLE TEMPORAL CELLS")
    print(f"Eligible cells: {len(eligible_cells)}")
    print(f"Eligible temporal pairs: {len(eligible_pairs)}")
    print(f"GEOJSON={geojson_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
