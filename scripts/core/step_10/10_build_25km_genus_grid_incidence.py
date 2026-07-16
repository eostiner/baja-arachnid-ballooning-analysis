#!/usr/bin/env python3
"""
Build aligned 25-km genus × grid-cell incidence matrices for the regenerated
Baja Ballooning final and sensitivity datasets.

The grid uses a spherical Lambert azimuthal equal-area projection centered on
Baja California (27.5 N, 113.5 W). Grid squares are exactly 25,000 m on each
side in projected space. The biodiversity-final dataset defines the common
row and column universe so all matrices are directly comparable.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

SCRIPT_VERSION = "10_v1_2026-07-13"
EARTH_RADIUS_M = 6_371_008.8
DEFAULT_CENTER_LAT = 27.5
DEFAULT_CENTER_LON = -113.5
DEFAULT_CELL_SIZE_M = 25_000.0
FLAGGED_FESA_GBIF_ID = "3095884508"

DATASET_SPECS = [
    {
        "key": "biodiversity_final",
        "required": True,
        "candidates": ["05_biodiversity_final_records.tsv"],
    },
    {
        "key": "ballooning_final",
        "required": True,
        "candidates": ["05_ballooning_final_records.tsv"],
    },
    {
        "key": "environmental_strict",
        "required": False,
        "candidates": [
            "05_environmental_model_final_records.tsv",
            "05_environmental_model_strict_records.tsv",
            "05_environmental_strict_records.tsv",
        ],
    },
    {
        "key": "environmental_inclusive",
        "required": False,
        "candidates": [
            "05_environmental_model_inclusive_sensitivity.tsv",
            "05_environmental_inclusive_sensitivity.tsv",
        ],
    },
    {
        "key": "environmental_temporal_2001_2020",
        "required": False,
        "candidates": [
            "05_environmental_temporal_sensitivity_2001_2020.tsv",
            "05_environmental_model_temporal_sensitivity_2001_2020.tsv",
            "05_temporal_sensitivity_2001_2020.tsv",
        ],
    },
    {
        "key": "biodiversity_taxonomy_strict",
        "required": False,
        "candidates": [
            "05_biodiversity_taxonomy_strict_sensitivity.tsv",
            "05_biodiversity_final_taxonomy_strict_sensitivity.tsv",
        ],
    },
    {
        "key": "ballooning_taxonomy_strict",
        "required": False,
        "candidates": [
            "05_ballooning_taxonomy_strict_sensitivity.tsv",
            "05_ballooning_final_taxonomy_strict_sensitivity.tsv",
        ],
    },
]

GENUS_FIELDS = ["analysis_genus", "genus", "trait_review_genus"]
LATITUDE_FIELDS = ["decimalLatitude", "latitude", "decimal_latitude"]
LONGITUDE_FIELDS = ["decimalLongitude", "longitude", "decimal_longitude"]
ID_FIELDS = ["gbifID", "gbifId", "key", "occurrenceID"]
ORDER_FIELDS = ["order", "trait_order", "analysis_order"]
FAMILY_FIELDS = ["family", "trait_family", "analysis_family"]


@dataclass(frozen=True, order=True)
class CellKey:
    row: int
    col: int

    @property
    def cell_id(self) -> str:
        return f"BJA25K_C{self.col:+05d}_R{self.row:+05d}"


@dataclass
class DatasetResult:
    key: str
    source_path: Path
    record_count: int
    unique_ids: set[str]
    genera: set[str]
    cells: set[CellKey]
    genus_cell_record_counts: Counter[tuple[str, CellKey]]
    genus_cell_ids: dict[tuple[str, CellKey], set[str]]
    genus_record_counts: Counter[str]
    genus_ids: dict[str, set[str]]
    cell_record_counts: Counter[CellKey]
    cell_genera: dict[CellKey, set[str]]
    genus_latitudes: dict[str, list[float]]
    genus_longitudes: dict[str, list[float]]
    genus_orders: dict[str, set[str]]
    genus_families: dict[str, set[str]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build 25-km genus-by-grid-cell incidence matrices."
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path.home() / "Desktop" / "Baja_Ballooning_Pipeline",
    )
    parser.add_argument("--cell-size-km", type=float, default=25.0)
    parser.add_argument("--center-lat", type=float, default=DEFAULT_CENTER_LAT)
    parser.add_argument("--center-lon", type=float, default=DEFAULT_CENTER_LON)
    return parser.parse_args()


def normalize(value: object) -> str:
    return str(value or "").strip()


def normalize_genus(value: object) -> str:
    value = normalize(value)
    if not value:
        return ""
    return value[0].upper() + value[1:]


def choose_field(fields: list[str], candidates: Iterable[str], required: bool = True) -> str | None:
    for candidate in candidates:
        if candidate in fields:
            return candidate
    if required:
        raise RuntimeError("None of these required fields were found: " + ", ".join(candidates))
    return None


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def write_csv(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    temporary.replace(path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def forward_laea(lat_deg: float, lon_deg: float, center_lat: float, center_lon: float) -> tuple[float, float]:
    phi = math.radians(lat_deg)
    lam = math.radians(lon_deg)
    phi0 = math.radians(center_lat)
    lam0 = math.radians(center_lon)
    dlam = lam - lam0
    denominator = 1.0 + math.sin(phi0) * math.sin(phi) + math.cos(phi0) * math.cos(phi) * math.cos(dlam)
    if denominator <= 0:
        raise ValueError("Coordinate is antipodal to the projection center.")
    k = math.sqrt(2.0 / denominator)
    x = EARTH_RADIUS_M * k * math.cos(phi) * math.sin(dlam)
    y = EARTH_RADIUS_M * k * (
        math.cos(phi0) * math.sin(phi)
        - math.sin(phi0) * math.cos(phi) * math.cos(dlam)
    )
    return x, y


def inverse_laea(x: float, y: float, center_lat: float, center_lon: float) -> tuple[float, float]:
    phi0 = math.radians(center_lat)
    lam0 = math.radians(center_lon)
    rho = math.hypot(x, y)
    if rho == 0:
        return center_lat, center_lon
    ratio = min(1.0, rho / (2.0 * EARTH_RADIUS_M))
    c = 2.0 * math.asin(ratio)
    sin_c = math.sin(c)
    cos_c = math.cos(c)
    lat = math.asin(
        cos_c * math.sin(phi0)
        + (y * sin_c * math.cos(phi0) / rho)
    )
    lon = lam0 + math.atan2(
        x * sin_c,
        rho * math.cos(phi0) * cos_c - y * math.sin(phi0) * sin_c,
    )
    lon_deg = (math.degrees(lon) + 540.0) % 360.0 - 180.0
    return math.degrees(lat), lon_deg


def cell_for_coordinate(lat: float, lon: float, cell_size_m: float, center_lat: float, center_lon: float) -> CellKey:
    x, y = forward_laea(lat, lon, center_lat, center_lon)
    return CellKey(row=math.floor(y / cell_size_m), col=math.floor(x / cell_size_m))


def cell_bounds(cell: CellKey, cell_size_m: float) -> tuple[float, float, float, float]:
    x_min = cell.col * cell_size_m
    x_max = (cell.col + 1) * cell_size_m
    y_min = cell.row * cell_size_m
    y_max = (cell.row + 1) * cell_size_m
    return x_min, y_min, x_max, y_max


def cell_centroid(cell: CellKey, cell_size_m: float, center_lat: float, center_lon: float) -> tuple[float, float, float, float]:
    x_min, y_min, x_max, y_max = cell_bounds(cell, cell_size_m)
    x = (x_min + x_max) / 2.0
    y = (y_min + y_max) / 2.0
    lat, lon = inverse_laea(x, y, center_lat, center_lon)
    return x, y, lat, lon


def latitude_band(lat: float) -> str:
    if 23 <= lat < 24:
        return "23-24N"
    if 24 <= lat < 26:
        return "24-26N"
    if 26 <= lat < 28:
        return "26-28N"
    if 28 <= lat < 30:
        return "28-30N"
    if 30 <= lat <= 32.6:
        return "30-32N"
    return "outside_study_bands"


def resolve_dataset(input_dir: Path, spec: dict[str, object]) -> Path | None:
    for filename in spec["candidates"]:
        candidate = input_dir / str(filename)
        if candidate.exists():
            return candidate
    if spec["required"]:
        raise FileNotFoundError(
            f"Required dataset for {spec['key']} not found. Tried:\n"
            + "\n".join(str(input_dir / str(name)) for name in spec["candidates"])
        )
    return None


def load_dataset(path: Path, key: str, cell_size_m: float, center_lat: float, center_lon: float) -> DatasetResult:
    fields, rows = read_tsv(path)
    genus_field = choose_field(fields, GENUS_FIELDS)
    lat_field = choose_field(fields, LATITUDE_FIELDS)
    lon_field = choose_field(fields, LONGITUDE_FIELDS)
    id_field = choose_field(fields, ID_FIELDS, required=False)
    order_field = choose_field(fields, ORDER_FIELDS, required=False)
    family_field = choose_field(fields, FAMILY_FIELDS, required=False)

    unique_ids: set[str] = set()
    genera: set[str] = set()
    cells: set[CellKey] = set()
    genus_cell_record_counts: Counter[tuple[str, CellKey]] = Counter()
    genus_cell_ids: dict[tuple[str, CellKey], set[str]] = defaultdict(set)
    genus_record_counts: Counter[str] = Counter()
    genus_ids: dict[str, set[str]] = defaultdict(set)
    cell_record_counts: Counter[CellKey] = Counter()
    cell_genera: dict[CellKey, set[str]] = defaultdict(set)
    genus_latitudes: dict[str, list[float]] = defaultdict(list)
    genus_longitudes: dict[str, list[float]] = defaultdict(list)
    genus_orders: dict[str, set[str]] = defaultdict(set)
    genus_families: dict[str, set[str]] = defaultdict(set)

    problems: list[str] = []
    for index, row in enumerate(rows, start=2):
        genus = normalize_genus(row.get(genus_field))
        if not genus:
            problems.append(f"row {index}: missing genus")
            continue
        try:
            lat = float(normalize(row.get(lat_field)))
            lon = float(normalize(row.get(lon_field)))
        except ValueError:
            problems.append(f"row {index}: invalid coordinates")
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            problems.append(f"row {index}: coordinates outside valid ranges")
            continue

        occurrence_id = normalize(row.get(id_field)) if id_field else f"ROW_{index}"
        cell = cell_for_coordinate(lat, lon, cell_size_m, center_lat, center_lon)
        pair = (genus, cell)

        unique_ids.add(occurrence_id)
        genera.add(genus)
        cells.add(cell)
        genus_cell_record_counts[pair] += 1
        genus_cell_ids[pair].add(occurrence_id)
        genus_record_counts[genus] += 1
        genus_ids[genus].add(occurrence_id)
        cell_record_counts[cell] += 1
        cell_genera[cell].add(genus)
        genus_latitudes[genus].append(lat)
        genus_longitudes[genus].append(lon)

        if order_field and normalize(row.get(order_field)):
            genus_orders[genus].add(normalize(row.get(order_field)))
        if family_field and normalize(row.get(family_field)):
            genus_families[genus].add(normalize(row.get(family_field)))

    if problems:
        preview = "\n".join(problems[:20])
        raise RuntimeError(
            f"{path.name} contains {len(problems)} unusable rows. First issues:\n{preview}"
        )

    if FLAGGED_FESA_GBIF_ID in unique_ids:
        raise RuntimeError(f"Flagged Fesa gbifID {FLAGGED_FESA_GBIF_ID} remains in {path.name}.")

    return DatasetResult(
        key=key,
        source_path=path,
        record_count=len(rows),
        unique_ids=unique_ids,
        genera=genera,
        cells=cells,
        genus_cell_record_counts=genus_cell_record_counts,
        genus_cell_ids=genus_cell_ids,
        genus_record_counts=genus_record_counts,
        genus_ids=genus_ids,
        cell_record_counts=cell_record_counts,
        cell_genera=cell_genera,
        genus_latitudes=genus_latitudes,
        genus_longitudes=genus_longitudes,
        genus_orders=genus_orders,
        genus_families=genus_families,
    )


def sort_cells(cells: Iterable[CellKey]) -> list[CellKey]:
    return sorted(cells, key=lambda cell: (cell.row, cell.col))


def write_common_grid(
    output_dir: Path,
    reference: DatasetResult,
    common_cells: list[CellKey],
    cell_size_m: float,
    center_lat: float,
    center_lon: float,
) -> tuple[Path, Path]:
    lookup_path = output_dir / "10_common_grid25km_cell_lookup.csv"
    geojson_path = output_dir / "10_common_grid25km_cells.geojson"
    lookup_rows = []
    features = []

    for position, cell in enumerate(common_cells, start=1):
        x_min, y_min, x_max, y_max = cell_bounds(cell, cell_size_m)
        centroid_x, centroid_y, centroid_lat, centroid_lon = cell_centroid(
            cell, cell_size_m, center_lat, center_lon
        )
        lookup_rows.append(
            {
                "grid_cell_order": position,
                "grid_cell_id": cell.cell_id,
                "grid_row": cell.row,
                "grid_column": cell.col,
                "x_min_m": round(x_min, 3),
                "y_min_m": round(y_min, 3),
                "x_max_m": round(x_max, 3),
                "y_max_m": round(y_max, 3),
                "centroid_x_m": round(centroid_x, 3),
                "centroid_y_m": round(centroid_y, 3),
                "centroid_latitude": round(centroid_lat, 7),
                "centroid_longitude": round(centroid_lon, 7),
                "centroid_latitude_band": latitude_band(centroid_lat),
                "biodiversity_record_count": reference.cell_record_counts[cell],
                "biodiversity_genus_richness": len(reference.cell_genera[cell]),
            }
        )

        corners_xy = [
            (x_min, y_min),
            (x_max, y_min),
            (x_max, y_max),
            (x_min, y_max),
            (x_min, y_min),
        ]
        ring = []
        for x, y in corners_xy:
            lat, lon = inverse_laea(x, y, center_lat, center_lon)
            ring.append([round(lon, 7), round(lat, 7)])
        features.append(
            {
                "type": "Feature",
                "properties": {
                    "grid_cell_order": position,
                    "grid_cell_id": cell.cell_id,
                    "grid_row": cell.row,
                    "grid_column": cell.col,
                    "cell_size_m": cell_size_m,
                    "centroid_latitude": round(centroid_lat, 7),
                    "centroid_longitude": round(centroid_lon, 7),
                    "centroid_latitude_band": latitude_band(centroid_lat),
                    "biodiversity_record_count": reference.cell_record_counts[cell],
                    "biodiversity_genus_richness": len(reference.cell_genera[cell]),
                },
                "geometry": {"type": "Polygon", "coordinates": [ring]},
            }
        )

    fields = list(lookup_rows[0].keys()) if lookup_rows else []
    write_csv(lookup_path, fields, lookup_rows)
    geojson_path.write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "name": "Baja_25km_common_occupied_grid",
                "crs_note": (
                    "Polygons are GeoJSON WGS84 coordinates derived from a spherical "
                    "Lambert azimuthal equal-area grid centered at 27.5N, -113.5E/W."
                ),
                "features": features,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return lookup_path, geojson_path


def write_common_genera(output_dir: Path, reference: DatasetResult, common_genera: list[str]) -> Path:
    path = output_dir / "10_common_genus_lookup.csv"
    rows = []
    for genus in common_genera:
        orders = sorted(reference.genus_orders.get(genus, set()))
        families = sorted(reference.genus_families.get(genus, set()))
        rows.append(
            {
                "genus_order": len(rows) + 1,
                "genus": genus,
                "order": " | ".join(orders),
                "family": " | ".join(families),
                "biodiversity_record_count": reference.genus_record_counts[genus],
                "biodiversity_unique_occurrence_ids": len(reference.genus_ids[genus]),
                "biodiversity_occupied_grid_cells": sum(
                    1 for candidate in reference.cells if (genus, candidate) in reference.genus_cell_record_counts
                ),
                "minimum_latitude": round(min(reference.genus_latitudes[genus]), 7),
                "maximum_latitude": round(max(reference.genus_latitudes[genus]), 7),
                "minimum_longitude": round(min(reference.genus_longitudes[genus]), 7),
                "maximum_longitude": round(max(reference.genus_longitudes[genus]), 7),
            }
        )
    write_csv(path, list(rows[0].keys()) if rows else [], rows)
    return path


def write_dataset_outputs(
    output_dir: Path,
    result: DatasetResult,
    common_genera: list[str],
    common_cells: list[CellKey],
    cell_size_m: float,
    center_lat: float,
    center_lon: float,
) -> dict[str, object]:
    matrix_path = output_dir / f"10_{result.key}_genus_by_grid25km_incidence.csv"
    long_path = output_dir / f"10_{result.key}_genus_grid25km_incidence_long.csv"
    genus_summary_path = output_dir / f"10_{result.key}_genus_summary.csv"
    cell_summary_path = output_dir / f"10_{result.key}_cell_summary.csv"

    matrix_fields = ["genus"] + [cell.cell_id for cell in common_cells]
    matrix_rows = []
    presence_count = 0
    for genus in common_genera:
        row: dict[str, object] = {"genus": genus}
        for cell in common_cells:
            value = 1 if (genus, cell) in result.genus_cell_record_counts else 0
            row[cell.cell_id] = value
            presence_count += value
        matrix_rows.append(row)
    write_csv(matrix_path, matrix_fields, matrix_rows)

    long_rows = []
    for (genus, cell), record_count in sorted(
        result.genus_cell_record_counts.items(), key=lambda item: (item[0][0].casefold(), item[0][1].row, item[0][1].col)
    ):
        _, _, centroid_lat, centroid_lon = cell_centroid(cell, cell_size_m, center_lat, center_lon)
        long_rows.append(
            {
                "genus": genus,
                "grid_cell_id": cell.cell_id,
                "incidence": 1,
                "occurrence_record_count": record_count,
                "unique_occurrence_ids": len(result.genus_cell_ids[(genus, cell)]),
                "grid_row": cell.row,
                "grid_column": cell.col,
                "centroid_latitude": round(centroid_lat, 7),
                "centroid_longitude": round(centroid_lon, 7),
                "centroid_latitude_band": latitude_band(centroid_lat),
            }
        )
    write_csv(long_path, list(long_rows[0].keys()) if long_rows else [], long_rows)

    genus_rows = []
    for genus in common_genera:
        occupied_cells = [cell for cell in common_cells if (genus, cell) in result.genus_cell_record_counts]
        latitudes = result.genus_latitudes.get(genus, [])
        longitudes = result.genus_longitudes.get(genus, [])
        genus_rows.append(
            {
                "genus": genus,
                "present_in_dataset": int(genus in result.genera),
                "occurrence_record_count": result.genus_record_counts.get(genus, 0),
                "unique_occurrence_ids": len(result.genus_ids.get(genus, set())),
                "occupied_grid_cells": len(occupied_cells),
                "minimum_latitude": round(min(latitudes), 7) if latitudes else "",
                "maximum_latitude": round(max(latitudes), 7) if latitudes else "",
                "minimum_longitude": round(min(longitudes), 7) if longitudes else "",
                "maximum_longitude": round(max(longitudes), 7) if longitudes else "",
            }
        )
    write_csv(genus_summary_path, list(genus_rows[0].keys()) if genus_rows else [], genus_rows)

    cell_rows = []
    for cell in common_cells:
        _, _, centroid_lat, centroid_lon = cell_centroid(cell, cell_size_m, center_lat, center_lon)
        cell_rows.append(
            {
                "grid_cell_id": cell.cell_id,
                "present_in_dataset": int(cell in result.cells),
                "occurrence_record_count": result.cell_record_counts.get(cell, 0),
                "genus_richness": len(result.cell_genera.get(cell, set())),
                "grid_row": cell.row,
                "grid_column": cell.col,
                "centroid_latitude": round(centroid_lat, 7),
                "centroid_longitude": round(centroid_lon, 7),
                "centroid_latitude_band": latitude_band(centroid_lat),
            }
        )
    write_csv(cell_summary_path, list(cell_rows[0].keys()) if cell_rows else [], cell_rows)

    return {
        "matrix_path": matrix_path,
        "long_path": long_path,
        "genus_summary_path": genus_summary_path,
        "cell_summary_path": cell_summary_path,
        "records": result.record_count,
        "unique_ids": len(result.unique_ids),
        "observed_genera": len(result.genera),
        "observed_cells": len(result.cells),
        "matrix_rows": len(common_genera),
        "matrix_columns": len(common_cells),
        "presence_count": presence_count,
        "matrix_density": presence_count / (len(common_genera) * len(common_cells)),
    }


def main() -> int:
    args = parse_args()
    project_root = args.project_root.expanduser().resolve()
    input_dir = project_root / "02_data_clean" / "05_final_qc_flags"
    output_dir = project_root / "02_data_clean" / "08_grid25km_incidence"
    archive_root = project_root / "08_archive"
    logs_dir = project_root / "06_logs"

    if not input_dir.exists():
        raise FileNotFoundError(f"Final QC folder not found:\n{input_dir}")

    cell_size_m = args.cell_size_km * 1000.0
    if cell_size_m <= 0:
        raise ValueError("Cell size must be positive.")

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    if output_dir.exists() and any(output_dir.iterdir()):
        archive_dir = archive_root / f"10_grid25km_incidence_{timestamp}"
        archive_dir.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(output_dir, archive_dir)
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_path = logs_dir / f"10_build_25km_incidence_{timestamp}.log"

    log_lines: list[str] = []
    def log(message: str = "") -> None:
        print(message)
        log_lines.append(message)

    log("Loading regenerated final datasets...")
    resolved: dict[str, Path] = {}
    results: dict[str, DatasetResult] = {}
    for spec in DATASET_SPECS:
        path = resolve_dataset(input_dir, spec)
        if path is None:
            log(f"Skipping optional dataset not found: {spec['key']}")
            continue
        resolved[str(spec["key"])] = path
        result = load_dataset(path, str(spec["key"]), cell_size_m, args.center_lat, args.center_lon)
        results[result.key] = result
        log(f"  {result.key}: {result.record_count:,} records, {len(result.genera):,} genera, {len(result.cells):,} occupied cells")

    reference = results["biodiversity_final"]
    ballooning = results["ballooning_final"]

    if reference.unique_ids != ballooning.unique_ids:
        raise RuntimeError("Biodiversity-final and ballooning-final occurrence ID sets differ.")
    if reference.genus_cell_record_counts != ballooning.genus_cell_record_counts:
        raise RuntimeError("Biodiversity-final and ballooning-final 25-km incidence structures differ.")
    if "Fesa" in reference.genera:
        raise RuntimeError("Fesa remains in the biodiversity-final genus set.")

    common_genera = sorted(reference.genera, key=str.casefold)
    common_cells = sort_cells(reference.cells)

    for key, result in results.items():
        if not result.unique_ids.issubset(reference.unique_ids):
            raise RuntimeError(f"{key} contains occurrence IDs absent from biodiversity_final.")
        if not result.genera.issubset(reference.genera):
            raise RuntimeError(f"{key} contains genera absent from biodiversity_final: {sorted(result.genera - reference.genera)}")
        if not result.cells.issubset(reference.cells):
            raise RuntimeError(f"{key} contains grid cells absent from biodiversity_final.")

    common_grid_lookup, common_grid_geojson = write_common_grid(
        output_dir, reference, common_cells, cell_size_m, args.center_lat, args.center_lon
    )
    common_genus_lookup = write_common_genera(output_dir, reference, common_genera)

    manifest_rows = []
    output_metadata: dict[str, object] = {}
    for key, result in results.items():
        metadata = write_dataset_outputs(
            output_dir, result, common_genera, common_cells, cell_size_m, args.center_lat, args.center_lon
        )
        output_metadata[key] = {
            **metadata,
            "source_path": str(result.source_path),
            "source_sha256": sha256_file(result.source_path),
        }
        manifest_rows.append(
            {
                "dataset": key,
                "source_file": result.source_path.name,
                "source_records": metadata["records"],
                "source_unique_ids": metadata["unique_ids"],
                "observed_genera": metadata["observed_genera"],
                "observed_grid_cells": metadata["observed_cells"],
                "aligned_matrix_rows": metadata["matrix_rows"],
                "aligned_matrix_grid_columns": metadata["matrix_columns"],
                "genus_cell_presences": metadata["presence_count"],
                "matrix_density": round(float(metadata["matrix_density"]), 8),
                "matrix_file": Path(metadata["matrix_path"]).name,
                "long_file": Path(metadata["long_path"]).name,
                "genus_summary_file": Path(metadata["genus_summary_path"]).name,
                "cell_summary_file": Path(metadata["cell_summary_path"]).name,
            }
        )

    manifest_path = output_dir / "10_grid25km_matrix_manifest.csv"
    write_csv(manifest_path, list(manifest_rows[0].keys()), manifest_rows)

    validation_rows = [
        {"check": "biodiversity_and_ballooning_id_sets_equal", "status": "PASS", "detail": str(reference.unique_ids == ballooning.unique_ids)},
        {"check": "biodiversity_and_ballooning_incidence_equal", "status": "PASS", "detail": str(reference.genus_cell_record_counts == ballooning.genus_cell_record_counts)},
        {"check": "flagged_fesa_gbifID_absent", "status": "PASS", "detail": FLAGGED_FESA_GBIF_ID},
        {"check": "fesa_genus_absent", "status": "PASS", "detail": str("Fesa" not in reference.genera)},
        {"check": "common_genus_count", "status": "PASS", "detail": len(common_genera)},
        {"check": "common_grid_cell_count", "status": "PASS", "detail": len(common_cells)},
        {"check": "biodiversity_record_count", "status": "PASS", "detail": reference.record_count},
        {"check": "ballooning_record_count", "status": "PASS", "detail": ballooning.record_count},
    ]
    validation_path = output_dir / "10_grid25km_validation.csv"
    write_csv(validation_path, ["check", "status", "detail"], validation_rows)

    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script_version": SCRIPT_VERSION,
        "script_path": str(Path(__file__).resolve()),
        "project_root": str(project_root),
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "grid_definition": {
            "cell_size_m": cell_size_m,
            "cell_size_km": args.cell_size_km,
            "projection": "Spherical Lambert azimuthal equal-area",
            "earth_radius_m": EARTH_RADIUS_M,
            "center_latitude": args.center_lat,
            "center_longitude": args.center_lon,
            "grid_origin_x_m": 0.0,
            "grid_origin_y_m": 0.0,
            "cell_index_rule": "col=floor(x/cell_size); row=floor(y/cell_size)",
            "cell_id_rule": "BJA25K_C{signed_col}_R{signed_row}",
            "common_cell_universe": "Cells occupied in biodiversity_final",
            "common_genus_universe": "Genera present in biodiversity_final",
        },
        "reference_dimensions": {
            "records": reference.record_count,
            "genera": len(common_genera),
            "occupied_grid_cells": len(common_cells),
        },
        "dataset_outputs": output_metadata,
        "common_outputs": {
            "cell_lookup": {"path": str(common_grid_lookup), "sha256": sha256_file(common_grid_lookup)},
            "cell_geojson": {"path": str(common_grid_geojson), "sha256": sha256_file(common_grid_geojson)},
            "genus_lookup": {"path": str(common_genus_lookup), "sha256": sha256_file(common_genus_lookup)},
            "manifest": {"path": str(manifest_path), "sha256": sha256_file(manifest_path)},
            "validation": {"path": str(validation_path), "sha256": sha256_file(validation_path)},
        },
        "incidence_definition": (
            "A matrix cell equals 1 when at least one retained occurrence of a genus falls "
            "within the corresponding 25-km grid cell; otherwise it equals 0. Multiple "
            "records in the same genus-cell combination are collapsed to one incidence."
        ),
    }
    provenance_path = output_dir / "10_grid25km_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2, default=str) + "\n", encoding="utf-8")

    log_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")

    log("\n25-km incidence matrices completed successfully.")
    log(f"Common matrix dimensions: {len(common_genera):,} genera × {len(common_cells):,} grid cells")
    log(f"Biodiversity genus-cell presences: {len(reference.genus_cell_record_counts):,}")
    log(f"Outputs: {output_dir}")
    log(f"Manifest: {manifest_path}")
    log(f"Provenance: {provenance_path}")
    log_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"\nERROR: {error}", file=sys.stderr)
        raise
