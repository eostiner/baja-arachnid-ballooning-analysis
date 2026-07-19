#!/usr/bin/env python3
"""
Phase 14C2 — validate and join temporally resolved cell-year stress anomalies to
the paired temporal turnover response from Phase 14B.

Only predictors frozen in phase_14_stress_predictors_frozen.csv are analyzed.
The primary complete-case gate depends only on predictors marked core=1, preventing
optional domain variables from deleting otherwise valid H3 pairs.
"""
from __future__ import annotations

import argparse
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean
from typing import Any

from phase14_common import default_analysis_output_root, load_periods, read_delimited, write_csv, write_json

SCRIPT_VERSION = "14C2_v0.2.0_2026-07-18"


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Join cell-year stress anomalies to Phase 14B turnover.")
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--cell-year-stress", type=Path)
    parser.add_argument("--turnover-file", type=Path)
    parser.add_argument("--period-config", type=Path, default=here / "configs" / "phase_14_temporal_windows_frozen.csv")
    parser.add_argument("--predictor-config", type=Path, default=here / "configs" / "phase_14_stress_predictors_frozen.csv")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--min-years-per-period", type=int, default=3)
    return parser.parse_args()


def to_float(value: Any) -> float | None:
    try:
        result = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def main() -> int:
    args = parse_args()
    root = args.project_root.expanduser().resolve()
    base = default_analysis_output_root(root)
    turnover_path = args.turnover_file.expanduser().resolve() if args.turnover_file else (
        base / "14B_temporal_community_turnover" / "14B_paired_temporal_turnover.csv"
    )
    stress_path = args.cell_year_stress.expanduser().resolve() if args.cell_year_stress else (
        base / "14C_real_stress_anomalies" / "14C1_cell_year_real_stress_anomalies.csv"
    )
    period_path = args.period_config.expanduser().resolve()
    predictor_path = args.predictor_config.expanduser().resolve()
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else base / "14C_temporal_stress_join"
    outdir.mkdir(parents=True, exist_ok=True)

    for path, label in ((turnover_path, "Phase 14B turnover"), (stress_path, "cell-year stress"), (predictor_path, "predictor config")):
        if not path.exists():
            raise FileNotFoundError(f"{label} file not found: {path}")

    _, turnover_rows = read_delimited(turnover_path)
    stress_fields, stress_rows = read_delimited(stress_path)
    _, predictor_rows = read_delimited(predictor_path)
    required = {"grid_cell_id", "year"}
    if not required.issubset(stress_fields):
        raise RuntimeError(f"Stress table missing required columns: {sorted(required - set(stress_fields))}")

    configured = [str(row.get("predictor", "")).strip() for row in predictor_rows if str(row.get("predictor", "")).strip()]
    missing_configured = [field for field in configured if field not in stress_fields]
    if missing_configured:
        raise RuntimeError(f"Stress table is missing configured predictors: {missing_configured}")
    numeric_fields = [field for field in configured if any(to_float(row.get(field)) is not None for row in stress_rows)]
    if not numeric_fields:
        raise RuntimeError("No numeric configured stress predictors were found.")
    core_fields = [
        str(row["predictor"]).strip()
        for row in predictor_rows
        if str(row.get("predictor", "")).strip() in numeric_fields and int(float(str(row.get("core", 0)))) == 1
    ]
    if not core_fields:
        raise RuntimeError("Predictor config contains no available core predictor.")

    periods = load_periods(period_path)
    period_lookup = {period["period_id"]: period for period in periods}
    annual: dict[tuple[str, int], dict[str, float]] = {}
    duplicates = 0
    for row in stress_rows:
        cell = str(row.get("grid_cell_id", "")).strip()
        try:
            year = int(float(str(row.get("year", "")).strip()))
        except ValueError:
            continue
        key = (cell, year)
        if key in annual:
            duplicates += 1
        annual[key] = {field: value for field in numeric_fields if (value := to_float(row.get(field))) is not None}
    if duplicates:
        raise RuntimeError(f"Stress table contains {duplicates} duplicate cell-year rows.")

    joined: list[dict[str, Any]] = []
    missing_rows: list[dict[str, Any]] = []
    coverage_rows: list[dict[str, Any]] = []
    for row in turnover_rows:
        cell = row["grid_cell_id"]
        p1 = period_lookup[row["period_1"]]
        p2 = period_lookup[row["period_2"]]
        output: dict[str, Any] = dict(row)
        failed_core: list[str] = []
        for field in numeric_fields:
            values1 = [
                annual[(cell, year)][field]
                for year in range(p1["start_year"], p1["end_year"] + 1)
                if (cell, year) in annual and field in annual[(cell, year)]
            ]
            values2 = [
                annual[(cell, year)][field]
                for year in range(p2["start_year"], p2["end_year"] + 1)
                if (cell, year) in annual and field in annual[(cell, year)]
            ]
            output[f"{field}_period1_mean"] = mean(values1) if values1 else ""
            output[f"{field}_period2_mean"] = mean(values2) if values2 else ""
            output[f"{field}_period1_n_years"] = len(values1)
            output[f"{field}_period2_n_years"] = len(values2)
            sufficient = len(values1) >= args.min_years_per_period and len(values2) >= args.min_years_per_period
            delta = mean(values2) - mean(values1) if sufficient else None
            output[f"delta_{field}_period2_minus_period1"] = delta if delta is not None else ""
            output[f"abs_delta_{field}_period2_minus_period1"] = abs(delta) if delta is not None else ""
            if field in core_fields and not sufficient:
                failed_core.append(field)
            coverage_rows.append(
                {
                    "pair_id": row.get("pair_id", ""),
                    "grid_cell_id": cell,
                    "period_1": row["period_1"],
                    "period_2": row["period_2"],
                    "predictor": field,
                    "core": int(field in core_fields),
                    "period1_years": len(values1),
                    "period2_years": len(values2),
                    "sufficient": int(sufficient),
                }
            )
        output["stress_join_complete"] = int(not failed_core)
        output["failed_core_predictors"] = ";".join(failed_core)
        joined.append(output)
        if failed_core:
            missing_rows.append(
                {
                    "pair_id": row.get("pair_id", ""),
                    "grid_cell_id": cell,
                    "period_1": row["period_1"],
                    "period_2": row["period_2"],
                    "failed_core_predictors": ";".join(failed_core),
                }
            )

    write_csv(outdir / "14C_H3_model_input_all_pairs.csv", joined)
    complete = [row for row in joined if int(row["stress_join_complete"]) == 1]
    write_csv(outdir / "14C_H3_model_input_complete_cases.csv", complete)
    write_csv(outdir / "14C_predictor_coverage_by_pair.csv", coverage_rows)
    write_csv(outdir / "14C_missing_core_stress_coverage.csv", missing_rows)
    status = {
        "phase": "14C2",
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETED" if complete else "NO_COMPLETE_CORE_STRESS_JOIN_CASES",
        "stress_columns": numeric_fields,
        "core_stress_columns": core_fields,
        "turnover_rows": len(turnover_rows),
        "complete_join_rows": len(complete),
        "minimum_years_per_period_per_predictor": args.min_years_per_period,
        "primary_predictor": "delta_stress_composite_z_period2_minus_period1",
        "secondary_predictor": "abs_delta_stress_composite_z_period2_minus_period1",
        "next_step": "Run Phase 14D only as a limited exploratory H3 test under the Phase 14A conditional status.",
    }
    write_json(outdir / "14C_run_status.json", status)
    print("PHASE 14C2 — TEMPORAL STRESS JOIN")
    print(f"Turnover rows: {len(turnover_rows)}")
    print(f"Complete core-stress joins: {len(complete)}")
    print(f"OUTPUT={outdir / '14C_H3_model_input_complete_cases.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
