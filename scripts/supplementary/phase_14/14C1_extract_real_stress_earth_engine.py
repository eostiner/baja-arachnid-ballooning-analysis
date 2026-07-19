#!/usr/bin/env python3
"""
Phase 14C1 — extract real annual climate and vegetation stress metrics for the
eligible temporal 25-km cells using Google Earth Engine, then calculate within-cell
standardized anomalies.

Primary sources:
  ECMWF/ERA5_LAND/DAILY_AGGR
  MODIS/061/MOD13Q1

Positive anomaly values are oriented to mean MORE environmental stress.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import statistics
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from phase14_common import default_analysis_output_root, read_delimited, write_csv, write_json

SCRIPT_VERSION = "14C1_v0.2.0_2026-07-18"
ERA5_ASSET = "ECMWF/ERA5_LAND/DAILY_AGGR"
MODIS_EVI_ASSET = "MODIS/061/MOD13Q1"
RAW_FIELDS = (
    "mean_temperature_c",
    "annual_max_temperature_c",
    "hot_days_35c",
    "mean_vpd_kpa",
    "annual_precipitation_mm",
    "mean_rootzone_soil_water",
    "mean_evi",
    "evi_p10",
    "mean_evi_observation_count",
)
ANOMALY_SOURCE_FIELDS = (
    "mean_temperature_c",
    "hot_days_35c",
    "mean_vpd_kpa",
    "annual_precipitation_mm",
    "mean_rootzone_soil_water",
    "mean_evi",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract real annual ERA5-Land and MODIS EVI stress anomalies.")
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--ee-project", default=os.environ.get("EARTHENGINE_PROJECT") or os.environ.get("GOOGLE_CLOUD_PROJECT"))
    parser.add_argument("--authenticate", action="store_true", help="Run the interactive Earth Engine authentication flow before initializing.")
    parser.add_argument("--cells-file", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--start-year", type=int, default=2001)
    parser.add_argument("--end-year", type=int, default=2025)
    parser.add_argument("--baseline-start", type=int, default=2001)
    parser.add_argument("--baseline-end", type=int, default=2020)
    parser.add_argument("--hot-threshold-c", type=float, default=35.0)
    parser.add_argument("--era5-scale-m", type=int, default=11132)
    parser.add_argument("--evi-scale-m", type=int, default=250)
    parser.add_argument("--tile-scale", type=int, default=4)
    return parser.parse_args()


def finite_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def initialize_earth_engine(project: str | None, authenticate: bool):
    try:
        import ee  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "The Earth Engine Python API is not installed. Run: "
            "python -m pip install --upgrade earthengine-api"
        ) from exc
    if authenticate:
        ee.Authenticate()
    try:
        if project:
            ee.Initialize(project=project)
        else:
            ee.Initialize()
    except Exception as exc:  # Earth Engine uses multiple exception classes across releases.
        message = (
            "Earth Engine initialization failed. Authenticate once and provide a Cloud project registered for Earth Engine.\n"
            "Example:\n"
            "  earthengine authenticate\n"
            "  python .../14C1_extract_real_stress_earth_engine.py --project-root ... --ee-project YOUR_PROJECT_ID\n"
            f"Original error: {exc}"
        )
        raise RuntimeError(message) from exc
    return ee


def load_cell_features(ee, path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Eligible-cell GeoJSON not found: {path}. Run Phase 14C0 first.")
    import json

    obj = json.loads(path.read_text(encoding="utf-8"))
    features = []
    rows = []
    for feature in obj.get("features", []):
        props = dict(feature.get("properties", {}))
        cell = str(props.get("grid_cell_id", "")).strip()
        coords = feature.get("geometry", {}).get("coordinates")
        if not cell or not coords:
            continue
        geometry = ee.Geometry.Polygon(coords, geodesic=False)
        features.append(ee.Feature(geometry, {"grid_cell_id": cell}))
        rows.append(props)
    if not features:
        raise RuntimeError(f"No valid eligible cell polygons found in {path}")
    return ee.FeatureCollection(features), rows


def vapor_pressure_deficit_image(ee, image):
    temperature_c = image.select("temperature_2m").subtract(273.15)
    dewpoint_c = image.select("dewpoint_temperature_2m").subtract(273.15)
    es = temperature_c.expression(
        "0.6108 * exp((17.27 * t) / (t + 237.3))", {"t": temperature_c}
    )
    ea = dewpoint_c.expression(
        "0.6108 * exp((17.27 * td) / (td + 237.3))", {"td": dewpoint_c}
    )
    return es.subtract(ea).max(0).rename("mean_vpd_kpa")


def annual_era5_image(ee, year: int, hot_threshold_c: float):
    start = ee.Date.fromYMD(year, 1, 1)
    end = start.advance(1, "year")
    daily = ee.ImageCollection(ERA5_ASSET).filterDate(start, end)

    def prep(image):
        mean_temp = image.select("temperature_2m").subtract(273.15).rename("mean_temperature_c")
        max_temp = image.select("temperature_2m_max").subtract(273.15).rename("daily_max_temperature_c")
        hot_day = max_temp.gt(hot_threshold_c).rename("hot_days_35c")
        vpd = vapor_pressure_deficit_image(ee, image)
        precipitation = image.select("total_precipitation_sum").max(ee.Image.constant(0)).multiply(1000).rename("daily_precipitation_mm")
        rootzone = image.select(["volumetric_soil_water_layer_1", "volumetric_soil_water_layer_2"]).reduce(ee.Reducer.mean()).rename("mean_rootzone_soil_water")
        return ee.Image.cat([mean_temp, max_temp, hot_day, vpd, precipitation, rootzone]).copyProperties(image, ["system:time_start"])

    prepared = daily.map(prep)
    annual = ee.Image.cat(
        [
            prepared.select("mean_temperature_c").mean(),
            prepared.select("daily_max_temperature_c").max().rename("annual_max_temperature_c"),
            prepared.select("hot_days_35c").sum(),
            prepared.select("mean_vpd_kpa").mean(),
            prepared.select("daily_precipitation_mm").sum().rename("annual_precipitation_mm"),
            prepared.select("mean_rootzone_soil_water").mean(),
        ]
    )
    return annual.set("year", year)


def annual_evi_image(ee, year: int):
    start = ee.Date.fromYMD(year, 1, 1)
    end = start.advance(1, "year")
    collection = ee.ImageCollection(MODIS_EVI_ASSET).filterDate(start, end)

    def prep(image):
        qa = image.select("SummaryQA").lte(1)
        return image.select("EVI").multiply(0.0001).updateMask(qa).rename("evi")

    evi = collection.map(prep)
    return ee.Image.cat(
        [
            evi.mean().rename("mean_evi"),
            evi.reduce(ee.Reducer.percentile([10])).rename("evi_p10"),
            evi.count().rename("mean_evi_observation_count"),
        ]
    ).set("year", year)


def reduce_to_cells(ee, image, cells_fc, scale: int, tile_scale: int) -> dict[str, dict[str, float | None]]:
    reduced = image.reduceRegions(
        collection=cells_fc,
        reducer=ee.Reducer.mean(),
        scale=scale,
        tileScale=tile_scale,
    )
    info = reduced.getInfo()
    result: dict[str, dict[str, float | None]] = {}
    for feature in info.get("features", []):
        props = feature.get("properties", {})
        cell = str(props.get("grid_cell_id", "")).strip()
        if not cell:
            continue
        result[cell] = {key: finite_float(value) for key, value in props.items() if key != "grid_cell_id"}
    return result


def z_parameters(values: Iterable[float]) -> tuple[float, float] | None:
    vals = [float(value) for value in values if value is not None and math.isfinite(float(value))]
    if len(vals) < 10:
        return None
    mean_value = statistics.fmean(vals)
    sd_value = statistics.stdev(vals)
    if not math.isfinite(sd_value) or sd_value <= 0:
        return None
    return mean_value, sd_value


def average_present(values: Iterable[float | None], minimum: int = 1) -> float | None:
    valid = [float(value) for value in values if value is not None and math.isfinite(float(value))]
    return statistics.fmean(valid) if len(valid) >= minimum else None


def main() -> int:
    args = parse_args()
    if args.start_year > args.end_year:
        raise ValueError("--start-year must be <= --end-year")
    if not (args.start_year <= args.baseline_start <= args.baseline_end <= args.end_year):
        raise ValueError("Baseline years must fall inside the extraction interval.")

    root = args.project_root.expanduser().resolve()
    base = default_analysis_output_root(root)
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else base / "14C_real_stress_anomalies"
    outdir.mkdir(parents=True, exist_ok=True)
    cells_path = args.cells_file.expanduser().resolve() if args.cells_file else outdir / "14C0_eligible_temporal_cells.geojson"

    ee = initialize_earth_engine(args.ee_project, args.authenticate)
    cells_fc, cell_rows = load_cell_features(ee, cells_path)
    cell_ids = sorted({str(row["grid_cell_id"]) for row in cell_rows})

    raw_rows: list[dict[str, Any]] = []
    for year in range(args.start_year, args.end_year + 1):
        print(f"Extracting real stress metrics for {year} ...", flush=True)
        era5 = reduce_to_cells(
            ee,
            annual_era5_image(ee, year, args.hot_threshold_c),
            cells_fc,
            scale=args.era5_scale_m,
            tile_scale=args.tile_scale,
        )
        evi = reduce_to_cells(
            ee,
            annual_evi_image(ee, year),
            cells_fc,
            scale=args.evi_scale_m,
            tile_scale=args.tile_scale,
        )
        for cell in cell_ids:
            row: dict[str, Any] = {"grid_cell_id": cell, "year": year}
            row.update({field: era5.get(cell, {}).get(field) for field in RAW_FIELDS if field in era5.get(cell, {})})
            row.update({field: evi.get(cell, {}).get(field) for field in RAW_FIELDS if field in evi.get(cell, {})})
            for field in RAW_FIELDS:
                row.setdefault(field, None)
            raw_rows.append(row)

    raw_path = outdir / "14C1_cell_year_raw_environment.csv"
    write_csv(raw_path, raw_rows)

    by_cell: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in raw_rows:
        by_cell[str(row["grid_cell_id"])].append(row)

    baseline_rows: list[dict[str, Any]] = []
    params: dict[tuple[str, str], tuple[float, float]] = {}
    for cell, rows in by_cell.items():
        baseline = [row for row in rows if args.baseline_start <= int(row["year"]) <= args.baseline_end]
        for field in ANOMALY_SOURCE_FIELDS:
            values = [finite_float(row.get(field)) for row in baseline]
            parameter = z_parameters(value for value in values if value is not None)
            baseline_rows.append(
                {
                    "grid_cell_id": cell,
                    "source_field": field,
                    "baseline_start": args.baseline_start,
                    "baseline_end": args.baseline_end,
                    "baseline_n_years": sum(value is not None for value in values),
                    "baseline_mean": parameter[0] if parameter else "",
                    "baseline_sd": parameter[1] if parameter else "",
                    "parameter_valid": int(parameter is not None),
                }
            )
            if parameter:
                params[(cell, field)] = parameter

    output_rows: list[dict[str, Any]] = []
    for row in raw_rows:
        cell = str(row["grid_cell_id"])
        output = dict(row)
        zvals: dict[str, float | None] = {}
        for field in ANOMALY_SOURCE_FIELDS:
            value = finite_float(row.get(field))
            parameter = params.get((cell, field))
            zvals[field] = (value - parameter[0]) / parameter[1] if value is not None and parameter else None

        output["mean_temperature_anomaly_z"] = zvals["mean_temperature_c"]
        output["hot_days_35c_anomaly_z"] = zvals["hot_days_35c"]
        output["vpd_anomaly_z"] = zvals["mean_vpd_kpa"]
        output["precipitation_deficit_z"] = -zvals["annual_precipitation_mm"] if zvals["annual_precipitation_mm"] is not None else None
        output["soil_water_deficit_z"] = -zvals["mean_rootzone_soil_water"] if zvals["mean_rootzone_soil_water"] is not None else None
        output["evi_deficit_z"] = -zvals["mean_evi"] if zvals["mean_evi"] is not None else None

        output["thermal_stress_z"] = average_present(
            [output["mean_temperature_anomaly_z"], output["hot_days_35c_anomaly_z"]], minimum=1
        )
        output["moisture_stress_z"] = average_present(
            [output["precipitation_deficit_z"], output["soil_water_deficit_z"]], minimum=1
        )
        output["vegetation_stress_z"] = output["evi_deficit_z"]
        domains = [
            output["thermal_stress_z"],
            output["vpd_anomaly_z"],
            output["moisture_stress_z"],
            output["vegetation_stress_z"],
        ]
        output["stress_domain_count"] = sum(value is not None for value in domains)
        output["stress_composite_z"] = average_present(domains, minimum=3)
        output["source"] = f"{ERA5_ASSET};{MODIS_EVI_ASSET}"
        output["baseline"] = f"within_cell_{args.baseline_start}_{args.baseline_end}"
        output_rows.append(output)

    output_path = outdir / "14C1_cell_year_real_stress_anomalies.csv"
    baseline_path = outdir / "14C1_baseline_parameters.csv"
    write_csv(output_path, output_rows)
    write_csv(baseline_path, baseline_rows)

    complete_composite = sum(finite_float(row.get("stress_composite_z")) is not None for row in output_rows)
    expected = len(cell_ids) * (args.end_year - args.start_year + 1)
    status = {
        "phase": "14C1",
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETED" if complete_composite == expected else "COMPLETED_WITH_MISSING_CELL_YEARS",
        "earth_engine_project": args.ee_project or "default credential project",
        "assets": [ERA5_ASSET, MODIS_EVI_ASSET],
        "eligible_cells": len(cell_ids),
        "years": [args.start_year, args.end_year],
        "expected_cell_year_rows": expected,
        "rows_written": len(output_rows),
        "complete_stress_composite_rows": complete_composite,
        "baseline": [args.baseline_start, args.baseline_end],
        "positive_values_mean": "more environmental stress",
        "hot_day_threshold_c": args.hot_threshold_c,
        "output": str(output_path),
        "important_guardrail": "These are temporally resolved observational anomalies, not experimental treatment assignments.",
    }
    write_json(outdir / "14C1_run_status.json", status)
    print("PHASE 14C1 — REAL STRESS ANOMALIES COMPLETED")
    print(f"Rows written: {len(output_rows)}")
    print(f"Complete composite rows: {complete_composite}/{expected}")
    print(f"OUTPUT={output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
