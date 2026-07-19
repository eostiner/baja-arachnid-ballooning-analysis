#!/usr/bin/env python3
"""
Phase 14A — audit whether retained final-QC occurrences can support a direct
within-cell temporal test of recent environmental-stress tracking (H3).

This script does not fit the H3 environmental model. It:
  * reads the retained final-QC occurrence table;
  * recovers collection year/date precision;
  * assigns the exact retained 25-km grid cell;
  * joins authoritative C3/N0 genus classifications;
  * quantifies records, events, genera, datasets, cells, and periods;
  * applies predeclared cell-period eligibility thresholds;
  * reports a transparent PASS / CONDITIONAL / FAIL feasibility result.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from phase14_common import (
    DATE_FIELDS,
    DATASET_FIELDS,
    DAY_FIELDS,
    DEFAULT_CELL_SIZE_M,
    DEFAULT_CENTER_LAT,
    DEFAULT_CENTER_LON,
    EVENT_ID_FIELDS,
    GENUS_FIELDS,
    ID_FIELDS,
    LAT_FIELDS,
    LON_FIELDS,
    MONTH_FIELDS,
    YEAR_FIELDS,
    cell_for_coordinate,
    clean,
    default_analysis_output_root,
    first_field,
    latitude_band,
    load_periods,
    load_trait_lookup,
    normalize_genus,
    parse_occurrence_date,
    period_for_year,
    read_delimited,
    resolve_occurrence_file,
    resolve_trait_file,
    sha256_file,
    write_csv,
    write_json,
)

SCRIPT_VERSION = "14A_v0.1.0_2026-07-18"


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Audit temporal coverage for Phase 14 H3.")
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--occurrence-file", type=Path)
    parser.add_argument("--trait-file", type=Path)
    parser.add_argument(
        "--period-config",
        type=Path,
        default=here / "configs" / "phase_14_temporal_windows_frozen.csv",
    )
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--min-records", type=int, default=5)
    parser.add_argument("--min-events", type=int, default=3)
    parser.add_argument("--min-genera", type=int, default=3)
    parser.add_argument("--min-direct-cells", type=int, default=20)
    parser.add_argument("--min-direct-comparisons", type=int, default=30)
    parser.add_argument("--min-conditional-cells", type=int, default=10)
    parser.add_argument("--min-conditional-comparisons", type=int, default=15)
    parser.add_argument("--cell-size-km", type=float, default=25.0)
    parser.add_argument("--center-lat", type=float, default=DEFAULT_CENTER_LAT)
    parser.add_argument("--center-lon", type=float, default=DEFAULT_CENTER_LON)
    return parser.parse_args()


def safe_float(value: Any) -> float | None:
    text = clean(value)
    if not text:
        return None
    try:
        result = float(text)
    except ValueError:
        return None
    return result if math.isfinite(result) else None


def event_key(
    row: dict[str, str],
    event_field: str | None,
    dataset_value: str,
    year: int,
    month: int | None,
    day: int | None,
    lat: float,
    lon: float,
    occurrence_id: str,
) -> str:
    explicit = clean(row.get(event_field)) if event_field else ""
    if explicit:
        return "EVENTID::" + explicit
    date_token = f"{year:04d}-{month or 0:02d}-{day or 0:02d}"
    # Conservative synthetic event: same dataset, date precision and locality.
    # Exact occurrence ID is retained only when no usable date-locality grouping exists.
    if month is not None:
        return f"SYNTH::{dataset_value}::{date_token}::{lat:.4f}::{lon:.4f}"
    return f"SYNTH_YEAR::{dataset_value}::{year:04d}::{lat:.4f}::{lon:.4f}"


def count_summary(records: list[dict[str, Any]], keys: tuple[str, ...]) -> list[dict[str, Any]]:
    grouped: dict[tuple[Any, ...], dict[str, Any]] = {}
    for record in records:
        group_key = tuple(record[key] for key in keys)
        entry = grouped.setdefault(
            group_key,
            {
                **{key: record[key] for key in keys},
                "occurrence_records": 0,
                "unique_events_set": set(),
                "unique_genera_set": set(),
                "unique_cells_set": set(),
                "unique_datasets_set": set(),
            },
        )
        entry["occurrence_records"] += 1
        entry["unique_events_set"].add(record["event_key"])
        entry["unique_genera_set"].add(record["genus"])
        entry["unique_cells_set"].add(record["grid_cell_id"])
        entry["unique_datasets_set"].add(record["dataset_key"])
    rows = []
    for entry in grouped.values():
        rows.append(
            {
                **{key: entry[key] for key in keys},
                "occurrence_records": entry["occurrence_records"],
                "unique_events": len(entry["unique_events_set"]),
                "unique_genera": len(entry["unique_genera_set"]),
                "unique_cells": len(entry["unique_cells_set"]),
                "unique_datasets": len(entry["unique_datasets_set"]),
            }
        )
    return sorted(rows, key=lambda row: tuple(str(row[key]) for key in keys))


def make_plots(outdir: Path, year_rows: list[dict[str, Any]], cell_period_rows: list[dict[str, Any]]) -> list[str]:
    outputs: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except ImportError:
        return outputs

    groups = ("C3", "N0")
    years = sorted({int(row["year"]) for row in year_rows if row["trait_group"] in groups})
    by_group = {(row["trait_group"], int(row["year"])): int(row["occurrence_records"]) for row in year_rows}
    if years:
        fig, ax = plt.subplots(figsize=(11, 5.5))
        for group in groups:
            ax.plot(years, [by_group.get((group, year), 0) for year in years], marker="o", markersize=2, linewidth=1, label=group)
        ax.set_xlabel("Collection year")
        ax.set_ylabel("Retained occurrence records")
        ax.set_title("Phase 14A temporal coverage by trait group")
        ax.legend()
        fig.tight_layout()
        for suffix in ("png", "svg"):
            path = outdir / f"14A_records_by_year.{suffix}"
            fig.savefig(path, dpi=300 if suffix == "png" else None)
            outputs.append(str(path))
        plt.close(fig)

    eligible = [row for row in cell_period_rows if int(row["paired_eligible"]) == 1]
    if eligible:
        fig, ax = plt.subplots(figsize=(7, 8))
        x = [float(row["centroid_longitude_mean"]) for row in eligible]
        y = [float(row["centroid_latitude_mean"]) for row in eligible]
        counts = Counter(row["grid_cell_id"] for row in eligible)
        sizes = [25 + 25 * counts[row["grid_cell_id"]] for row in eligible]
        scatter = ax.scatter(x, y, s=sizes, c=[counts[row["grid_cell_id"]] for row in eligible])
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title("Cell-periods eligible for paired C3/N0 temporal analysis")
        fig.colorbar(scatter, ax=ax, label="Eligible periods in cell")
        fig.tight_layout()
        for suffix in ("png", "svg"):
            path = outdir / f"14A_eligible_temporal_cells.{suffix}"
            fig.savefig(path, dpi=300 if suffix == "png" else None)
            outputs.append(str(path))
        plt.close(fig)
    return outputs


def main() -> int:
    args = parse_args()
    project_root = args.project_root.expanduser().resolve()
    if not project_root.is_dir():
        raise FileNotFoundError(f"Project root not found: {project_root}")
    occurrence_path = resolve_occurrence_file(project_root, args.occurrence_file)
    trait_path = resolve_trait_file(project_root, args.trait_file)
    period_path = args.period_config.expanduser().resolve()
    periods = load_periods(period_path)
    analysis_start = min(period["start_year"] for period in periods)
    analysis_end = max(period["end_year"] for period in periods)
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else (
        default_analysis_output_root(project_root) / "14A_temporal_feasibility_audit"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    fields, rows = read_delimited(occurrence_path)
    genus_field = first_field(fields, GENUS_FIELDS, required=True)
    lat_field = first_field(fields, LAT_FIELDS, required=True)
    lon_field = first_field(fields, LON_FIELDS, required=True)
    id_field = first_field(fields, ID_FIELDS, required=False)
    date_field = first_field(fields, DATE_FIELDS, required=False)
    year_field = first_field(fields, YEAR_FIELDS, required=False)
    month_field = first_field(fields, MONTH_FIELDS, required=False)
    day_field = first_field(fields, DAY_FIELDS, required=False)
    event_field = first_field(fields, EVENT_ID_FIELDS, required=False)
    dataset_field = first_field(fields, DATASET_FIELDS, required=False)
    if not date_field and not year_field:
        raise RuntimeError(
            "The final-QC occurrence table contains neither a recognized eventDate field nor a year field. "
            "H3 cannot be reconstructed from this table without recovering collection dates upstream."
        )

    traits = load_trait_lookup(trait_path)
    counters: Counter[str] = Counter()
    date_precision: Counter[str] = Counter()
    analysis_records: list[dict[str, Any]] = []
    all_valid_years: list[int] = []
    cell_coords: dict[str, list[tuple[float, float]]] = defaultdict(list)

    current_year = datetime.now().year
    for row_number, row in enumerate(rows, start=2):
        counters["input_rows"] += 1
        genus = normalize_genus(row.get(genus_field))
        lat = safe_float(row.get(lat_field))
        lon = safe_float(row.get(lon_field))
        if not genus:
            counters["missing_genus"] += 1
            continue
        if lat is None or lon is None or not (-90 <= lat <= 90 and -180 <= lon <= 180):
            counters["invalid_coordinates"] += 1
            continue
        trait = traits.get(genus, {"trait_group": "UNRESOLVED", "evidence_class": ""})
        trait_group = trait["trait_group"]
        counters[f"trait_{trait_group}"] += 1
        parsed = parse_occurrence_date(
            row, date_field, year_field, month_field, day_field, current_year=current_year
        )
        date_precision[parsed.precision] += 1
        if parsed.year is None:
            counters["missing_or_invalid_date"] += 1
            continue
        all_valid_years.append(parsed.year)
        if not analysis_start <= parsed.year <= analysis_end:
            counters["outside_primary_window"] += 1
            continue
        period_id = period_for_year(parsed.year, periods)
        if period_id is None:
            counters["no_period_assignment"] += 1
            continue
        if trait_group not in {"C3", "N0"}:
            counters["excluded_nonprimary_trait"] += 1
            continue
        occurrence_id = clean(row.get(id_field)) if id_field else f"ROW_{row_number}"
        if not occurrence_id:
            occurrence_id = f"ROW_{row_number}"
        dataset_value = clean(row.get(dataset_field)) if dataset_field else "UNSPECIFIED_DATASET"
        if not dataset_value:
            dataset_value = "UNSPECIFIED_DATASET"
        cell = cell_for_coordinate(
            lat,
            lon,
            cell_size_m=args.cell_size_km * 1000.0,
            center_lat=args.center_lat,
            center_lon=args.center_lon,
        )
        key = event_key(
            row, event_field, dataset_value, parsed.year, parsed.month, parsed.day,
            lat, lon, occurrence_id,
        )
        record = {
            "occurrence_id": occurrence_id,
            "genus": genus,
            "trait_group": trait_group,
            "evidence_class": trait["evidence_class"],
            "year": parsed.year,
            "month": parsed.month or "",
            "day": parsed.day or "",
            "date_precision": parsed.precision,
            "date_source_field": parsed.source,
            "period_id": period_id,
            "grid_cell_id": cell.cell_id,
            "grid_row": cell.row,
            "grid_column": cell.col,
            "decimalLatitude": round(lat, 7),
            "decimalLongitude": round(lon, 7),
            "latitude_band": latitude_band(lat),
            "event_key": key,
            "dataset_key": dataset_value,
        }
        analysis_records.append(record)
        cell_coords[cell.cell_id].append((lat, lon))
        counters["primary_window_records"] += 1

    # Validate exact compatibility with the retained Step 10 common grid when available.
    grid_lookup_candidates = [
        project_root / "02_data_clean" / "08_grid25km_incidence" / "10_common_grid25km_cell_lookup.csv",
        project_root / "ANALYSIS_READY_INPUTS" / "02_incidence_matrices_25km" / "10_common_grid25km_cell_lookup.csv",
    ]
    grid_lookup_path = next((candidate for candidate in grid_lookup_candidates if candidate.exists()), None)
    grid_validation = {
        "status": "SKIPPED_LOOKUP_NOT_FOUND",
        "lookup_path": "",
        "computed_unique_cells": len({record["grid_cell_id"] for record in analysis_records}),
        "cells_absent_from_retained_lookup": [],
    }
    if grid_lookup_path is not None:
        lookup_fields, lookup_rows = read_delimited(grid_lookup_path)
        if "grid_cell_id" not in lookup_fields:
            raise RuntimeError(f"Retained grid lookup lacks grid_cell_id: {grid_lookup_path}")
        retained_cells = {clean(row.get("grid_cell_id")) for row in lookup_rows}
        computed_cells = {record["grid_cell_id"] for record in analysis_records}
        absent = sorted(computed_cells - retained_cells)
        grid_validation = {
            "status": "PASS" if not absent else "FAIL",
            "lookup_path": str(grid_lookup_path),
            "lookup_sha256": sha256_file(grid_lookup_path),
            "retained_lookup_cells": len(retained_cells),
            "computed_unique_cells": len(computed_cells),
            "cells_absent_from_retained_lookup": absent,
        }
        if absent:
            raise RuntimeError(
                "Phase 14 computed cells absent from the retained Step 10 grid lookup: "
                + ", ".join(absent[:20])
            )

    record_index_path = outdir / "14A_record_temporal_index.tsv"
    record_fields = list(analysis_records[0].keys()) if analysis_records else [
        "occurrence_id", "genus", "trait_group", "evidence_class", "year", "month", "day",
        "date_precision", "date_source_field", "period_id", "grid_cell_id", "grid_row",
        "grid_column", "decimalLatitude", "decimalLongitude", "latitude_band", "event_key",
        "dataset_key",
    ]
    # TSV is used because dataset titles may contain commas.
    with record_index_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=record_fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(analysis_records)

    year_rows = count_summary(analysis_records, ("year", "trait_group"))
    write_csv(outdir / "14A_year_trait_coverage.csv", year_rows)
    dataset_rows = count_summary(analysis_records, ("dataset_key", "trait_group"))
    write_csv(outdir / "14A_dataset_trait_coverage.csv", dataset_rows)

    group_summary: dict[tuple[str, str, str], dict[str, Any]] = {}
    for record in analysis_records:
        key = (record["grid_cell_id"], record["period_id"], record["trait_group"])
        entry = group_summary.setdefault(
            key,
            {
                "records": 0,
                "events": set(),
                "genera": set(),
                "datasets": set(),
                "bands": set(),
            },
        )
        entry["records"] += 1
        entry["events"].add(record["event_key"])
        entry["genera"].add(record["genus"])
        entry["datasets"].add(record["dataset_key"])
        entry["bands"].add(record["latitude_band"])

    all_cell_periods = sorted({(cell, period) for cell, period, _ in group_summary})
    cell_period_rows = []
    for cell_id, period_id in all_cell_periods:
        c3 = group_summary.get((cell_id, period_id, "C3"), {})
        n0 = group_summary.get((cell_id, period_id, "N0"), {})
        coords = cell_coords[cell_id]
        c3_records = int(c3.get("records", 0))
        n0_records = int(n0.get("records", 0))
        c3_events = len(c3.get("events", set()))
        n0_events = len(n0.get("events", set()))
        c3_genera = len(c3.get("genera", set()))
        n0_genera = len(n0.get("genera", set()))
        c3_eligible = (
            c3_records >= args.min_records
            and c3_events >= args.min_events
            and c3_genera >= args.min_genera
        )
        n0_eligible = (
            n0_records >= args.min_records
            and n0_events >= args.min_events
            and n0_genera >= args.min_genera
        )
        cell_period_rows.append(
            {
                "grid_cell_id": cell_id,
                "period_id": period_id,
                "centroid_latitude_mean": round(sum(x[0] for x in coords) / len(coords), 7),
                "centroid_longitude_mean": round(sum(x[1] for x in coords) / len(coords), 7),
                "latitude_band": latitude_band(sum(x[0] for x in coords) / len(coords)),
                "c3_records": c3_records,
                "c3_unique_events": c3_events,
                "c3_genera": c3_genera,
                "c3_datasets": len(c3.get("datasets", set())),
                "n0_records": n0_records,
                "n0_unique_events": n0_events,
                "n0_genera": n0_genera,
                "n0_datasets": len(n0.get("datasets", set())),
                "c3_eligible": int(c3_eligible),
                "n0_eligible": int(n0_eligible),
                "paired_eligible": int(c3_eligible and n0_eligible),
            }
        )
    write_csv(outdir / "14A_cell_period_coverage.csv", cell_period_rows)

    eligibility = {(row["grid_cell_id"], row["period_id"]): int(row["paired_eligible"]) == 1 for row in cell_period_rows}
    period_ids = [str(period["period_id"]) for period in periods]
    pair_rows = []
    for cell_id in sorted({row["grid_cell_id"] for row in cell_period_rows}):
        row_lookup = {(row["grid_cell_id"], row["period_id"]): row for row in cell_period_rows}
        for left, right in zip(period_ids, period_ids[1:]):
            left_ok = eligibility.get((cell_id, left), False)
            right_ok = eligibility.get((cell_id, right), False)
            representative = row_lookup.get((cell_id, left)) or row_lookup.get((cell_id, right))
            pair_rows.append(
                {
                    "grid_cell_id": cell_id,
                    "period_1": left,
                    "period_2": right,
                    "latitude_band": representative["latitude_band"] if representative else "",
                    "period_1_paired_eligible": int(left_ok),
                    "period_2_paired_eligible": int(right_ok),
                    "temporal_pair_eligible": int(left_ok and right_ok),
                }
            )
    write_csv(outdir / "14A_adjacent_period_pair_eligibility.csv", pair_rows)

    eligible_pairs = [row for row in pair_rows if int(row["temporal_pair_eligible"]) == 1]
    eligible_cells = sorted({row["grid_cell_id"] for row in eligible_pairs})
    eligible_bands = sorted({row["latitude_band"] for row in eligible_pairs if row["latitude_band"] != "outside_study_bands"})
    direct_pass = (
        len(eligible_cells) >= args.min_direct_cells
        and len(eligible_pairs) >= args.min_direct_comparisons
        and len(eligible_bands) >= 4
    )
    conditional_pass = (
        len(eligible_cells) >= args.min_conditional_cells
        and len(eligible_pairs) >= args.min_conditional_comparisons
        and len(eligible_bands) >= 3
    )
    if direct_pass:
        status = "PASS_DIRECT_H3_TEMPORAL_TEST"
        recommendation = "Proceed to Phase 14B temporal turnover and Phase 14C stress-anomaly construction."
    elif conditional_pass:
        status = "CONDITIONAL_LIMITED_H3_TEST"
        recommendation = (
            "Proceed only as an explicitly limited/exploratory H3 test, with equal-event resampling, "
            "dataset-continuity sensitivity analyses, and cautious geographic scope."
        )
    else:
        status = "FAIL_INSUFFICIENT_REPEATED_TEMPORAL_COVERAGE"
        recommendation = (
            "Do not fit a peninsula-wide H3 model from the present opportunistic records. "
            "Consider broader periods, dataset-restricted tests, or standardized resurveys."
        )

    sensitivity_rows = []
    for min_events in (1, 2, 3, 5):
        for min_genera in (2, 3, 5):
            paired_ok = {}
            for row in cell_period_rows:
                c3_ok = int(row["c3_records"]) >= args.min_records and int(row["c3_unique_events"]) >= min_events and int(row["c3_genera"]) >= min_genera
                n0_ok = int(row["n0_records"]) >= args.min_records and int(row["n0_unique_events"]) >= min_events and int(row["n0_genera"]) >= min_genera
                paired_ok[(row["grid_cell_id"], row["period_id"])] = c3_ok and n0_ok
            candidate_pairs = [
                row for row in pair_rows
                if paired_ok.get((row["grid_cell_id"], row["period_1"]), False)
                and paired_ok.get((row["grid_cell_id"], row["period_2"]), False)
            ]
            sensitivity_rows.append(
                {
                    "min_records_fixed": args.min_records,
                    "min_events_per_group_cell_period": min_events,
                    "min_genera_per_group_cell_period": min_genera,
                    "eligible_temporal_pairs": len(candidate_pairs),
                    "eligible_cells": len({row["grid_cell_id"] for row in candidate_pairs}),
                    "latitude_bands": len({row["latitude_band"] for row in candidate_pairs if row["latitude_band"] != "outside_study_bands"}),
                }
            )
    write_csv(outdir / "14A_threshold_sensitivity.csv", sensitivity_rows)

    go_no_go_rows = [
        {"criterion": "primary_window", "value": f"{analysis_start}-{analysis_end}", "threshold": "frozen config", "pass": 1},
        {"criterion": "valid_primary_C3_N0_records", "value": len(analysis_records), "threshold": ">0", "pass": int(bool(analysis_records))},
        {"criterion": "eligible_temporal_pairs", "value": len(eligible_pairs), "threshold": args.min_direct_comparisons, "pass": int(len(eligible_pairs) >= args.min_direct_comparisons)},
        {"criterion": "eligible_cells", "value": len(eligible_cells), "threshold": args.min_direct_cells, "pass": int(len(eligible_cells) >= args.min_direct_cells)},
        {"criterion": "eligible_latitude_bands", "value": len(eligible_bands), "threshold": 4, "pass": int(len(eligible_bands) >= 4)},
        {"criterion": "overall_status", "value": status, "threshold": "operational feasibility rule", "pass": int(direct_pass)},
    ]
    write_csv(outdir / "14A_go_no_go_summary.csv", go_no_go_rows)

    field_audit = [
        {"role": "genus", "selected_field": genus_field, "required": 1},
        {"role": "latitude", "selected_field": lat_field, "required": 1},
        {"role": "longitude", "selected_field": lon_field, "required": 1},
        {"role": "occurrence_id", "selected_field": id_field or "synthetic row id", "required": 0},
        {"role": "event_date", "selected_field": date_field or "", "required": 0},
        {"role": "year", "selected_field": year_field or "", "required": 0},
        {"role": "month", "selected_field": month_field or "", "required": 0},
        {"role": "day", "selected_field": day_field or "", "required": 0},
        {"role": "event_id", "selected_field": event_field or "synthetic event key", "required": 0},
        {"role": "dataset", "selected_field": dataset_field or "UNSPECIFIED_DATASET", "required": 0},
    ]
    write_csv(outdir / "14A_field_audit.csv", field_audit)

    plot_outputs = make_plots(outdir, year_rows, cell_period_rows)
    valid_date_fraction = len(all_valid_years) / len(rows) if rows else 0.0
    summary = {
        "phase": "14A",
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "recommendation": recommendation,
        "project_root": str(project_root),
        "inputs": {
            "occurrence_file": str(occurrence_path),
            "occurrence_sha256": sha256_file(occurrence_path),
            "trait_file": str(trait_path),
            "trait_sha256": sha256_file(trait_path),
            "period_config": str(period_path),
            "period_config_sha256": sha256_file(period_path),
        },
        "selected_fields": {row["role"]: row["selected_field"] for row in field_audit},
        "grid_validation": grid_validation,
        "date_audit": {
            "valid_date_fraction_all_input_rows": valid_date_fraction,
            "minimum_valid_year": min(all_valid_years) if all_valid_years else None,
            "maximum_valid_year": max(all_valid_years) if all_valid_years else None,
            "precision_counts": dict(date_precision),
        },
        "primary_window": {"start_year": analysis_start, "end_year": analysis_end, "periods": periods},
        "eligibility_thresholds": {
            "min_records_per_trait_cell_period": args.min_records,
            "min_unique_events_per_trait_cell_period": args.min_events,
            "min_genera_per_trait_cell_period": args.min_genera,
            "direct_min_cells": args.min_direct_cells,
            "direct_min_adjacent_comparisons": args.min_direct_comparisons,
            "direct_min_latitude_bands": 4,
            "conditional_min_cells": args.min_conditional_cells,
            "conditional_min_adjacent_comparisons": args.min_conditional_comparisons,
            "conditional_min_latitude_bands": 3,
        },
        "counts": {
            **dict(counters),
            "analysis_records_written": len(analysis_records),
            "eligible_cell_periods": sum(int(row["paired_eligible"]) for row in cell_period_rows),
            "eligible_adjacent_temporal_pairs": len(eligible_pairs),
            "eligible_cells": len(eligible_cells),
            "eligible_latitude_bands": eligible_bands,
        },
        "plot_outputs": plot_outputs,
        "important_interpretive_rule": (
            "PASS indicates operational feasibility for an observational temporal test; it does not establish "
            "causation or guarantee statistical power. FAIL prohibits a peninsula-wide confirmatory H3 model "
            "under the frozen thresholds but does not prove the mechanism absent."
        ),
    }
    write_json(outdir / "14A_temporal_feasibility_summary.json", summary)
    readme = f"""PHASE 14A — TEMPORAL FEASIBILITY AUDIT
======================================
Status: {status}

Input occurrence table: {occurrence_path}
Input trait table: {trait_path}
Primary analysis window: {analysis_start}-{analysis_end}
Retained C3/N0 records in primary window: {len(analysis_records):,}
Eligible adjacent cell-period comparisons: {len(eligible_pairs):,}
Eligible cells: {len(eligible_cells):,}
Eligible latitude bands: {', '.join(eligible_bands) if eligible_bands else 'none'}

Eligibility per trait group within each cell-period:
  records >= {args.min_records}
  unique events >= {args.min_events}
  genera >= {args.min_genera}

Recommendation:
{recommendation}

Interpretation:
This audit determines whether the existing occurrence data can support a repeated,
within-cell temporal comparison. It does not itself test stress tracking and does not
convert opportunistic GBIF records into standardized surveys.
"""
    (outdir / "14A_README.txt").write_text(readme, encoding="utf-8")
    print(readme)
    print(f"OUTPUT_DIR={outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
