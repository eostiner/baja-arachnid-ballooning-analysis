#!/usr/bin/env python3
"""Phase 14E2 — summarize annual extended stress within each five-year period and join to turnover."""
from __future__ import annotations

import argparse
import math
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean
from typing import Any

from phase14_common import default_analysis_output_root, load_periods, read_delimited, write_csv, write_json

SCRIPT_VERSION = "14E2_v0.3.0_2026-07-18"


def to_float(value: Any) -> float | None:
    try:
        result = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def summarize(values: list[float]) -> dict[str, float]:
    positive = [max(0.0, value) for value in values]
    return {
        "mean": mean(values),
        "max": max(values),
        "burden": mean(positive),
        "extreme_fraction": sum(value >= 1.0 for value in values) / len(values),
    }


def main() -> int:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Join extended period stress summaries to temporal turnover.")
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--cell-year-stress", type=Path)
    parser.add_argument("--turnover-file", type=Path)
    parser.add_argument("--period-config", type=Path, default=here / "configs" / "phase_14_temporal_windows_frozen.csv")
    parser.add_argument("--predictor-config", type=Path, default=here / "configs" / "phase_14_extended_predictors_frozen.csv")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--min-years-per-period", type=int, default=3)
    args = parser.parse_args()

    root = args.project_root.expanduser().resolve()
    base = default_analysis_output_root(root)
    stress_path = args.cell_year_stress.expanduser().resolve() if args.cell_year_stress else (
        base / "14E_extended_stress_anomalies" / "14E1_cell_year_extended_stress_anomalies.csv"
    )
    turnover_path = args.turnover_file.expanduser().resolve() if args.turnover_file else (
        base / "14B_temporal_community_turnover" / "14B_paired_temporal_turnover.csv"
    )
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else base / "14E_extended_stress_join"
    outdir.mkdir(parents=True, exist_ok=True)

    for path, label in ((stress_path, "extended annual stress"), (turnover_path, "temporal turnover"), (args.predictor_config, "predictor config")):
        if not path.exists():
            raise FileNotFoundError(f"{label} file not found: {path}")

    stress_fields, stress_rows = read_delimited(stress_path)
    _, turnover_rows = read_delimited(turnover_path)
    _, predictor_rows = read_delimited(args.predictor_config.expanduser().resolve())
    periods = {period["period_id"]: period for period in load_periods(args.period_config.expanduser().resolve())}

    predictors = [str(row.get("predictor", "")).strip() for row in predictor_rows if str(row.get("predictor", "")).strip()]
    missing = [predictor for predictor in predictors if predictor not in stress_fields]
    if missing:
        raise RuntimeError(f"Extended stress table is missing configured predictors: {missing}")
    core_predictors = [
        str(row["predictor"]).strip() for row in predictor_rows
        if int(float(str(row.get("core", 0)))) == 1
    ]

    annual: dict[tuple[str, int], dict[str, float]] = {}
    for row in stress_rows:
        cell = str(row.get("grid_cell_id", "")).strip()
        try:
            year = int(float(str(row.get("year", "")).strip()))
        except ValueError:
            continue
        key = (cell, year)
        if key in annual:
            raise RuntimeError(f"Duplicate extended stress row for {cell}, {year}")
        annual[key] = {
            predictor: value
            for predictor in predictors
            if (value := to_float(row.get(predictor))) is not None
        }

    joined: list[dict[str, Any]] = []
    coverage_rows: list[dict[str, Any]] = []
    for turnover in turnover_rows:
        cell = str(turnover["grid_cell_id"])
        p1 = periods[str(turnover["period_1"])]
        p2 = periods[str(turnover["period_2"])]
        output: dict[str, Any] = dict(turnover)
        failed_core: list[str] = []
        for predictor in predictors:
            values_by_period: list[list[float]] = []
            for period in (p1, p2):
                values = [
                    annual[(cell, year)][predictor]
                    for year in range(int(period["start_year"]), int(period["end_year"]) + 1)
                    if (cell, year) in annual and predictor in annual[(cell, year)]
                ]
                values_by_period.append(values)
            sufficient = all(len(values) >= args.min_years_per_period for values in values_by_period)
            coverage_rows.append({
                "pair_id": turnover.get("pair_id", ""),
                "grid_cell_id": cell,
                "period_1": turnover["period_1"],
                "period_2": turnover["period_2"],
                "predictor": predictor,
                "core": int(predictor in core_predictors),
                "period1_years": len(values_by_period[0]),
                "period2_years": len(values_by_period[1]),
                "sufficient": int(sufficient),
            })
            if not sufficient:
                if predictor in core_predictors:
                    failed_core.append(predictor)
                for summary_name in ("mean", "max", "burden", "extreme_fraction"):
                    output[f"{predictor}_period1_{summary_name}"] = ""
                    output[f"{predictor}_period2_{summary_name}"] = ""
                    output[f"delta_{summary_name}_{predictor}_period2_minus_period1"] = ""
                continue
            summary1 = summarize(values_by_period[0])
            summary2 = summarize(values_by_period[1])
            for summary_name in ("mean", "max", "burden", "extreme_fraction"):
                output[f"{predictor}_period1_{summary_name}"] = summary1[summary_name]
                output[f"{predictor}_period2_{summary_name}"] = summary2[summary_name]
                output[f"delta_{summary_name}_{predictor}_period2_minus_period1"] = summary2[summary_name] - summary1[summary_name]
        output["extended_stress_join_complete"] = int(not failed_core)
        output["failed_core_predictors"] = ";".join(failed_core)
        joined.append(output)

    complete = [row for row in joined if int(row["extended_stress_join_complete"]) == 1]
    write_csv(outdir / "14E2_H3_extended_model_input_all_pairs.csv", joined)
    write_csv(outdir / "14E2_H3_extended_model_input_primary_complete.csv", complete)
    write_csv(outdir / "14E2_predictor_coverage_by_pair.csv", coverage_rows)
    status = {
        "phase": "14E2",
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETED" if complete else "NO_PRIMARY_COMPLETE_CASES",
        "turnover_rows": len(turnover_rows),
        "primary_complete_rows": len(complete),
        "predictors": predictors,
        "summary_modes": ["mean", "max", "burden", "extreme_fraction"],
        "primary_predictor": "delta_mean_extended_stress_composite_z_period2_minus_period1",
        "minimum_years_per_period": args.min_years_per_period,
        "guardrail": "Alternative summaries are prespecified sensitivities and must not be cherry-picked.",
    }
    write_json(outdir / "14E2_run_status.json", status)
    print("PHASE 14E2 — EXTENDED STRESS JOIN COMPLETED")
    print(f"Turnover rows: {len(turnover_rows)}")
    print(f"Primary complete rows: {len(complete)}")
    print(f"OUTPUT={outdir / '14E2_H3_extended_model_input_primary_complete.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
