#!/usr/bin/env python3
"""Phase 14E1 — extract independent and extreme-event stress metrics for H3.

Adds ERA5-Land extremes, Terra+Aqua MODIS land-surface temperature, CHIRPS
rainfall drought/extremes, TerraClimate water balance, and MODIS EVI extremes.
All annual stress variables are standardized within cell against 2001–2020 and
oriented so positive values mean greater stress.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from phase14_common import default_analysis_output_root, write_csv, write_json

SCRIPT_VERSION = "14E1_v0.3.0_2026-07-18"
ERA5_ASSET = "ECMWF/ERA5_LAND/DAILY_AGGR"
MODIS_TERRA_LST = "MODIS/061/MOD11A1"
MODIS_AQUA_LST = "MODIS/061/MYD11A1"
MODIS_EVI_ASSET = "MODIS/061/MOD13Q1"
CHIRPS_ASSET = "UCSB-CHG/CHIRPS/DAILY"
TERRACLIMATE_ASSET = "IDAHO_EPSCOR/TERRACLIMATE"

RAW_ORIENTATION = {
    # field: +1 means high raw values are stress; -1 means low raw values are stress.
    "era5_mean_temperature_c": 1,
    "era5_p95_daily_max_temperature_c": 1,
    "era5_annual_max_temperature_c": 1,
    "era5_hot_days_35c": 1,
    "era5_hot_days_40c": 1,
    "era5_mean_vpd_kpa": 1,
    "era5_p95_vpd_kpa": 1,
    "era5_high_vpd_day_fraction": 1,
    "era5_mean_rootzone_soil_water": -1,
    "modis_lst_day_mean_c": 1,
    "modis_lst_day_p95_c": 1,
    "modis_lst_day_max_c": 1,
    "modis_lst_night_mean_c": 1,
    "modis_lst_hot_observation_fraction_45c": 1,
    "chirps_annual_precipitation_mm": -1,
    "chirps_dry_day_fraction_lt1mm": 1,
    "chirps_very_dry_day_fraction_lt0_1mm": 1,
    "chirps_minimum_monthly_precipitation_mm": -1,
    "chirps_heavy_rain_days_gt20mm": 1,
    "chirps_max_1day_precipitation_mm": 1,
    "terraclimate_climate_water_deficit_mm": 1,
    "terraclimate_pdsi_mean": -1,
    "terraclimate_pdsi_min": -1,
    "terraclimate_soil_moisture_mean_mm": -1,
    "terraclimate_aet_pet_ratio": -1,
    "modis_evi_mean": -1,
    "modis_evi_p10": -1,
}

DOMAIN_MEMBERS = {
    "air_heat_extreme_stress_z": [
        "era5_mean_temperature_c_stress_z",
        "era5_p95_daily_max_temperature_c_stress_z",
        "era5_hot_days_35c_stress_z",
        "era5_hot_days_40c_stress_z",
    ],
    "surface_heat_stress_z": [
        "modis_lst_day_mean_c_stress_z",
        "modis_lst_day_p95_c_stress_z",
        "modis_lst_hot_observation_fraction_45c_stress_z",
        "modis_lst_night_mean_c_stress_z",
    ],
    "vpd_extreme_stress_z": [
        "era5_mean_vpd_kpa_stress_z",
        "era5_p95_vpd_kpa_stress_z",
        "era5_high_vpd_day_fraction_stress_z",
    ],
    "rainfall_drought_stress_z": [
        "chirps_annual_precipitation_mm_stress_z",
        "chirps_dry_day_fraction_lt1mm_stress_z",
        "chirps_minimum_monthly_precipitation_mm_stress_z",
    ],
    "soil_water_stress_z": ["era5_mean_rootzone_soil_water_stress_z"],
    "water_balance_stress_z": [
        "terraclimate_climate_water_deficit_mm_stress_z",
        "terraclimate_pdsi_mean_stress_z",
        "terraclimate_soil_moisture_mean_mm_stress_z",
        "terraclimate_aet_pet_ratio_stress_z",
    ],
    "vegetation_extreme_stress_z": [
        "modis_evi_mean_stress_z",
        "modis_evi_p10_stress_z",
    ],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract extended annual environmental stress metrics.")
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--ee-project", default=os.environ.get("EE_PROJECT") or os.environ.get("EARTHENGINE_PROJECT") or os.environ.get("GOOGLE_CLOUD_PROJECT"))
    parser.add_argument("--authenticate", action="store_true")
    parser.add_argument("--cells-file", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--start-year", type=int, default=2001)
    parser.add_argument("--end-year", type=int, default=2025)
    parser.add_argument("--baseline-start", type=int, default=2001)
    parser.add_argument("--baseline-end", type=int, default=2020)
    parser.add_argument("--era5-scale-m", type=int, default=11132)
    parser.add_argument("--modis-scale-m", type=int, default=1000)
    parser.add_argument("--evi-scale-m", type=int, default=250)
    parser.add_argument("--chirps-scale-m", type=int, default=5566)
    parser.add_argument("--terraclimate-scale-m", type=int, default=4638)
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
        raise RuntimeError("Install Earth Engine with: python -m pip install --upgrade earthengine-api") from exc
    if authenticate:
        ee.Authenticate(auth_mode="localhost")
    if project:
        ee.Initialize(project=project)
    else:
        ee.Initialize()
    return ee


def load_cell_features(ee, path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Eligible-cell GeoJSON not found: {path}. Run 14C0 first.")
    obj = json.loads(path.read_text(encoding="utf-8"))
    features = []
    rows = []
    for feature in obj.get("features", []):
        props = dict(feature.get("properties", {}))
        cell = str(props.get("grid_cell_id", "")).strip()
        coords = feature.get("geometry", {}).get("coordinates")
        if cell and coords:
            features.append(ee.Feature(ee.Geometry.Polygon(coords, geodesic=False), {"grid_cell_id": cell}))
            rows.append(props)
    if not features:
        raise RuntimeError(f"No valid eligible cells in {path}")
    return ee.FeatureCollection(features), rows


def reduce_to_cells(ee, image, cells_fc, scale: int, tile_scale: int) -> dict[str, dict[str, float | None]]:
    reduced = image.reduceRegions(collection=cells_fc, reducer=ee.Reducer.mean(), scale=scale, tileScale=tile_scale)
    info = reduced.getInfo()
    output: dict[str, dict[str, float | None]] = {}
    for feature in info.get("features", []):
        props = feature.get("properties", {})
        cell = str(props.get("grid_cell_id", "")).strip()
        if cell:
            output[cell] = {key: finite_float(value) for key, value in props.items() if key != "grid_cell_id"}
    return output


def vpd_image(ee, image):
    t = image.select("temperature_2m").subtract(273.15)
    td = image.select("dewpoint_temperature_2m").subtract(273.15)
    es = t.expression("0.6108 * exp((17.27 * x) / (x + 237.3))", {"x": t})
    ea = td.expression("0.6108 * exp((17.27 * x) / (x + 237.3))", {"x": td})
    return es.subtract(ea).max(0).rename("vpd")


def annual_era5_image(ee, year: int):
    start = ee.Date.fromYMD(year, 1, 1)
    end = start.advance(1, "year")
    daily = ee.ImageCollection(ERA5_ASSET).filterDate(start, end)

    def prep(image):
        tmean = image.select("temperature_2m").subtract(273.15).rename("tmean")
        tmax = image.select("temperature_2m_max").subtract(273.15).rename("tmax")
        vpd = vpd_image(ee, image)
        soil = image.select(["volumetric_soil_water_layer_1", "volumetric_soil_water_layer_2"]).reduce(ee.Reducer.mean()).rename("soil")
        return ee.Image.cat([
            tmean,
            tmax,
            tmax.gt(35).rename("hot35"),
            tmax.gt(40).rename("hot40"),
            vpd,
            vpd.gt(3).rename("high_vpd"),
            soil,
        ]).copyProperties(image, ["system:time_start"])

    prepared = daily.map(prep)
    return ee.Image.cat([
        prepared.select("tmean").mean().rename("era5_mean_temperature_c"),
        prepared.select("tmax").reduce(ee.Reducer.percentile([95])).rename("era5_p95_daily_max_temperature_c"),
        prepared.select("tmax").max().rename("era5_annual_max_temperature_c"),
        prepared.select("hot35").sum().rename("era5_hot_days_35c"),
        prepared.select("hot40").sum().rename("era5_hot_days_40c"),
        prepared.select("vpd").mean().rename("era5_mean_vpd_kpa"),
        prepared.select("vpd").reduce(ee.Reducer.percentile([95])).rename("era5_p95_vpd_kpa"),
        prepared.select("high_vpd").mean().rename("era5_high_vpd_day_fraction"),
        prepared.select("soil").mean().rename("era5_mean_rootzone_soil_water"),
    ])


def prep_modis_lst(ee, image):
    day_qa = image.select("QC_Day")
    night_qa = image.select("QC_Night")
    day_mask = (
        day_qa.bitwiseAnd(3).lte(1)
        .And(day_qa.rightShift(2).bitwiseAnd(3).eq(0))
        .And(day_qa.rightShift(6).bitwiseAnd(3).lte(1))
    )
    night_mask = (
        night_qa.bitwiseAnd(3).lte(1)
        .And(night_qa.rightShift(2).bitwiseAnd(3).eq(0))
        .And(night_qa.rightShift(6).bitwiseAnd(3).lte(1))
    )
    day = image.select("LST_Day_1km").multiply(0.02).subtract(273.15).updateMask(day_mask).rename("lst_day")
    night = image.select("LST_Night_1km").multiply(0.02).subtract(273.15).updateMask(night_mask).rename("lst_night")
    return ee.Image.cat([
        day,
        night,
        day.gt(45).rename("hot45"),
        day.mask().rename("day_valid"),
        night.mask().rename("night_valid"),
    ]).copyProperties(image, ["system:time_start"])


def annual_modis_lst_image(ee, year: int):
    start = ee.Date.fromYMD(year, 1, 1)
    end = start.advance(1, "year")
    terra = ee.ImageCollection(MODIS_TERRA_LST).filterDate(start, end).map(lambda image: prep_modis_lst(ee, image))
    aqua = ee.ImageCollection(MODIS_AQUA_LST).filterDate(start, end).map(lambda image: prep_modis_lst(ee, image))
    combined = terra.merge(aqua)
    return ee.Image.cat([
        combined.select("lst_day").mean().rename("modis_lst_day_mean_c"),
        combined.select("lst_day").reduce(ee.Reducer.percentile([95])).rename("modis_lst_day_p95_c"),
        combined.select("lst_day").max().rename("modis_lst_day_max_c"),
        combined.select("lst_night").mean().rename("modis_lst_night_mean_c"),
        combined.select("hot45").mean().rename("modis_lst_hot_observation_fraction_45c"),
        combined.select("day_valid").sum().rename("modis_lst_day_observation_count"),
        combined.select("night_valid").sum().rename("modis_lst_night_observation_count"),
    ])


def annual_chirps_image(ee, year: int):
    start = ee.Date.fromYMD(year, 1, 1)
    end = start.advance(1, "year")
    daily = ee.ImageCollection(CHIRPS_ASSET).filterDate(start, end).select("precipitation")
    monthly = []
    for month in range(1, 13):
        mstart = ee.Date.fromYMD(year, month, 1)
        mend = mstart.advance(1, "month")
        monthly.append(daily.filterDate(mstart, mend).sum().rename("monthly_precip"))
    monthly_collection = ee.ImageCollection.fromImages(monthly)
    return ee.Image.cat([
        daily.sum().rename("chirps_annual_precipitation_mm"),
        daily.map(lambda image: image.lt(1).rename("dry")).mean().rename("chirps_dry_day_fraction_lt1mm"),
        daily.map(lambda image: image.lt(0.1).rename("very_dry")).mean().rename("chirps_very_dry_day_fraction_lt0_1mm"),
        monthly_collection.min().rename("chirps_minimum_monthly_precipitation_mm"),
        daily.map(lambda image: image.gt(20).rename("heavy")).sum().rename("chirps_heavy_rain_days_gt20mm"),
        daily.max().rename("chirps_max_1day_precipitation_mm"),
    ])


def annual_evi_image(ee, year: int):
    start = ee.Date.fromYMD(year, 1, 1)
    end = start.advance(1, "year")
    collection = ee.ImageCollection(MODIS_EVI_ASSET).filterDate(start, end)

    def prep(image):
        qa = image.select("SummaryQA").lte(1)
        return image.select("EVI").multiply(0.0001).updateMask(qa).rename("evi")

    evi = collection.map(prep)
    return ee.Image.cat([
        evi.mean().rename("modis_evi_mean"),
        evi.reduce(ee.Reducer.percentile([10])).rename("modis_evi_p10"),
        evi.count().rename("modis_evi_observation_count"),
    ])


def annual_terraclimate_image(ee, year: int):
    start = ee.Date.fromYMD(year, 1, 1)
    end = start.advance(1, "year")
    monthly = ee.ImageCollection(TERRACLIMATE_ASSET).filterDate(start, end)
    deficit = monthly.select("def").sum().multiply(0.1).rename("terraclimate_climate_water_deficit_mm")
    pdsi = monthly.select("pdsi").mean().multiply(0.01).rename("terraclimate_pdsi_mean")
    pdsi_min = monthly.select("pdsi").min().multiply(0.01).rename("terraclimate_pdsi_min")
    soil = monthly.select("soil").mean().multiply(0.1).rename("terraclimate_soil_moisture_mean_mm")
    aet = monthly.select("aet").sum().multiply(0.1)
    pet = monthly.select("pet").sum().multiply(0.1)
    ratio = aet.divide(pet.max(ee.Image.constant(0.001))).rename("terraclimate_aet_pet_ratio")
    return ee.Image.cat([deficit, pdsi, pdsi_min, soil, ratio])


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
    if not (args.start_year <= args.baseline_start <= args.baseline_end <= args.end_year):
        raise ValueError("Baseline years must lie within extraction years.")
    root = args.project_root.expanduser().resolve()
    base = default_analysis_output_root(root)
    c1_dir = base / "14C_real_stress_anomalies"
    cells_path = args.cells_file.expanduser().resolve() if args.cells_file else c1_dir / "14C0_eligible_temporal_cells.geojson"
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else base / "14E_extended_stress_anomalies"
    outdir.mkdir(parents=True, exist_ok=True)

    ee = initialize_earth_engine(args.ee_project, args.authenticate)
    cells_fc, cell_rows = load_cell_features(ee, cells_path)
    cell_ids = sorted({str(row["grid_cell_id"]) for row in cell_rows})
    tc_latest_ms = ee.ImageCollection(TERRACLIMATE_ASSET).aggregate_max("system:time_start").getInfo()
    tc_latest_year = datetime.fromtimestamp(float(tc_latest_ms) / 1000, tz=timezone.utc).year if tc_latest_ms else 0

    raw_rows: list[dict[str, Any]] = []
    for year in range(args.start_year, args.end_year + 1):
        print(f"Extracting extended stress metrics for {year} ...", flush=True)
        datasets = {
            "era5": reduce_to_cells(ee, annual_era5_image(ee, year), cells_fc, args.era5_scale_m, args.tile_scale),
            "modis_lst": reduce_to_cells(ee, annual_modis_lst_image(ee, year), cells_fc, args.modis_scale_m, args.tile_scale),
            "chirps": reduce_to_cells(ee, annual_chirps_image(ee, year), cells_fc, args.chirps_scale_m, args.tile_scale),
            "evi": reduce_to_cells(ee, annual_evi_image(ee, year), cells_fc, args.evi_scale_m, args.tile_scale),
        }
        if year <= tc_latest_year:
            datasets["terraclimate"] = reduce_to_cells(
                ee, annual_terraclimate_image(ee, year), cells_fc, args.terraclimate_scale_m, args.tile_scale
            )
        else:
            datasets["terraclimate"] = {}
        for cell in cell_ids:
            row: dict[str, Any] = {"grid_cell_id": cell, "year": year}
            for values in datasets.values():
                row.update(values.get(cell, {}))
            for field in RAW_ORIENTATION:
                row.setdefault(field, None)
            row.setdefault("modis_lst_day_observation_count", None)
            row.setdefault("modis_lst_night_observation_count", None)
            row.setdefault("modis_evi_observation_count", None)
            raw_rows.append(row)

    raw_path = outdir / "14E1_cell_year_raw_extended_environment.csv"
    write_csv(raw_path, raw_rows)

    by_cell: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in raw_rows:
        by_cell[str(row["grid_cell_id"])].append(row)

    params: dict[tuple[str, str], tuple[float, float]] = {}
    baseline_rows: list[dict[str, Any]] = []
    for cell, rows in by_cell.items():
        baseline = [row for row in rows if args.baseline_start <= int(row["year"]) <= args.baseline_end]
        for field in RAW_ORIENTATION:
            values = [finite_float(row.get(field)) for row in baseline]
            parameter = z_parameters(value for value in values if value is not None)
            baseline_rows.append({
                "grid_cell_id": cell,
                "source_field": field,
                "baseline_start": args.baseline_start,
                "baseline_end": args.baseline_end,
                "baseline_n_years": sum(value is not None for value in values),
                "baseline_mean": parameter[0] if parameter else "",
                "baseline_sd": parameter[1] if parameter else "",
                "orientation": RAW_ORIENTATION[field],
                "parameter_valid": int(parameter is not None),
            })
            if parameter:
                params[(cell, field)] = parameter

    output_rows: list[dict[str, Any]] = []
    for row in raw_rows:
        cell = str(row["grid_cell_id"])
        output = dict(row)
        for field, orientation in RAW_ORIENTATION.items():
            value = finite_float(row.get(field))
            parameter = params.get((cell, field))
            stress = orientation * (value - parameter[0]) / parameter[1] if value is not None and parameter else None
            output[f"{field}_stress_z"] = stress
        for domain, members in DOMAIN_MEMBERS.items():
            minimum = 1 if len(members) == 1 else max(1, len(members) // 2)
            output[domain] = average_present([finite_float(output.get(member)) for member in members], minimum=minimum)
        output["thermal_extreme_stress_z"] = average_present([
            finite_float(output.get("air_heat_extreme_stress_z")),
            finite_float(output.get("surface_heat_stress_z")),
        ], minimum=1)
        output["moisture_extreme_stress_z"] = average_present([
            finite_float(output.get("rainfall_drought_stress_z")),
            finite_float(output.get("soil_water_stress_z")),
            finite_float(output.get("water_balance_stress_z")),
        ], minimum=2)
        composite_domains = [
            finite_float(output.get("thermal_extreme_stress_z")),
            finite_float(output.get("vpd_extreme_stress_z")),
            finite_float(output.get("moisture_extreme_stress_z")),
            finite_float(output.get("vegetation_extreme_stress_z")),
        ]
        output["extended_stress_domain_count"] = sum(value is not None for value in composite_domains)
        output["extended_stress_composite_z"] = average_present(composite_domains, minimum=3)
        output["source"] = ";".join([
            ERA5_ASSET, MODIS_TERRA_LST, MODIS_AQUA_LST, MODIS_EVI_ASSET, CHIRPS_ASSET, TERRACLIMATE_ASSET
        ])
        output["baseline"] = f"within_cell_{args.baseline_start}_{args.baseline_end}"
        output_rows.append(output)

    output_path = outdir / "14E1_cell_year_extended_stress_anomalies.csv"
    write_csv(output_path, output_rows)
    write_csv(outdir / "14E1_baseline_parameters.csv", baseline_rows)

    expected = len(cell_ids) * (args.end_year - args.start_year + 1)
    complete = sum(finite_float(row.get("extended_stress_composite_z")) is not None for row in output_rows)
    status = {
        "phase": "14E1",
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETED" if complete >= expected - len(cell_ids) else "COMPLETED_WITH_MISSING_VALUES",
        "earth_engine_project": args.ee_project or "default",
        "assets": [ERA5_ASSET, MODIS_TERRA_LST, MODIS_AQUA_LST, MODIS_EVI_ASSET, CHIRPS_ASSET, TERRACLIMATE_ASSET],
        "eligible_cells": len(cell_ids),
        "years": [args.start_year, args.end_year],
        "expected_rows": expected,
        "rows_written": len(output_rows),
        "complete_extended_composite_rows": complete,
        "terraclimate_latest_available_year_detected": tc_latest_year,
        "guardrail": "Independent products and extremes improve stress measurement but do not increase biological replication.",
        "output": str(output_path),
    }
    write_json(outdir / "14E1_run_status.json", status)
    print("PHASE 14E1 — EXTENDED STRESS EXTRACTION COMPLETED")
    print(f"Rows written: {len(output_rows)}")
    print(f"Complete extended composite rows: {complete}/{expected}")
    print(f"OUTPUT={output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
