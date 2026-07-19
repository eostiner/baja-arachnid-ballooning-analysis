#!/usr/bin/env python3
"""Phase 15A — freeze and audit the paired Bayesian H3 model input."""
from __future__ import annotations

import argparse
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from phase15_common import analysis_root, as_float, period_midpoint, phase15_root, read_csv, write_csv, write_json

SCRIPT_VERSION = "15A_v0.1.0_2026-07-18"
RESPONSE = "resampled_delta_simpson_C3_minus_N0_median"
Q025 = "resampled_delta_simpson_C3_minus_N0_q025"
Q975 = "resampled_delta_simpson_C3_minus_N0_q975"
PREDICTOR = "delta_stress_composite_z_period2_minus_period1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare Phase 15 Bayesian H3 input.")
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--input-file", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--min-pairs", type=int, default=12)
    parser.add_argument("--min-cells", type=int, default=10)
    parser.add_argument("--minimum-observation-sd", type=float, default=0.01)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.project_root.expanduser().resolve()
    if not root.exists():
        raise FileNotFoundError(f"Project root not found: {root}")
    source = args.input_file.expanduser().resolve() if args.input_file else (
        analysis_root(root) / "14C_temporal_stress_join" / "14C_H3_model_input_complete_cases.csv"
    )
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else phase15_root(root) / "15A_input_audit"
    outdir.mkdir(parents=True, exist_ok=True)
    if not source.exists():
        raise FileNotFoundError(
            f"Pre-specified Phase 14C primary input not found: {source}\n"
            "Phase 15 deliberately uses the original Phase 14 composite as primary rather than selecting an extended predictor after seeing results."
        )
    rows = read_csv(source)
    required = {
        "pair_id", "grid_cell_id", "period_1", "period_2", "common_event_resample_n",
        RESPONSE, Q025, Q975, PREDICTOR,
    }
    missing = required - set(rows[0] if rows else {})
    if missing:
        raise RuntimeError(f"Input missing required columns: {sorted(missing)}")

    provisional = []
    raw_ses = []
    excluded = []
    for row in rows:
        y = as_float(row.get(RESPONSE))
        lo = as_float(row.get(Q025))
        hi = as_float(row.get(Q975))
        stress = as_float(row.get(PREDICTOR))
        events = as_float(row.get("common_event_resample_n"))
        if None in (y, lo, hi, stress, events) or events is None or events < 1:
            excluded.append({"pair_id": row.get("pair_id", ""), "reason": "missing_or_invalid_required_value"})
            continue
        if not (-1.000001 <= y <= 1.000001):
            excluded.append({"pair_id": row.get("pair_id", ""), "reason": "response_outside_minus1_plus1"})
            continue
        se = max(0.0, (float(hi) - float(lo)) / (2.0 * 1.959963984540054))
        if se > 0:
            raw_ses.append(se)
        provisional.append((row, float(y), float(lo), float(hi), float(stress), float(events), se))

    adaptive_floor = args.minimum_observation_sd
    if raw_ses:
        adaptive_floor = max(adaptive_floor, 0.10 * float(np.median(raw_ses)))
    prepared = []
    floor_count = 0
    for row, y, lo, hi, stress, events, se in provisional:
        used = max(se, adaptive_floor)
        floor_count += int(se < adaptive_floor)
        prepared.append({
            "pair_id": row["pair_id"],
            "grid_cell_id": row["grid_cell_id"],
            "period_1": row["period_1"],
            "period_2": row["period_2"],
            "latitude_band": row.get("latitude_band", ""),
            "y_observed": y,
            "resampling_q025": lo,
            "resampling_q975": hi,
            "observation_sd_raw_approx": se,
            "observation_sd_approx": used,
            "observation_sd_floor_applied": int(se < adaptive_floor),
            "stress_change_raw": stress,
            "common_event_resample_n": events,
            "log_common_events_raw": math.log1p(events),
            "transition_midpoint_raw": period_midpoint(row["period_1"], row["period_2"]),
        })

    cells = {row["grid_cell_id"] for row in prepared}
    status = "PASS" if len(prepared) >= args.min_pairs and len(cells) >= args.min_cells else "FAIL_INSUFFICIENT_DATA"
    write_csv(outdir / "15A_bayesian_model_input.csv", prepared)
    write_csv(outdir / "15A_excluded_rows.csv", excluded)
    audit = [
        {"item": "status", "value": status},
        {"item": "source_file", "value": str(source)},
        {"item": "eligible_pairs", "value": len(prepared)},
        {"item": "eligible_cells", "value": len(cells)},
        {"item": "observation_sd_floor", "value": adaptive_floor},
        {"item": "pairs_using_sd_floor", "value": floor_count},
        {"item": "primary_response", "value": RESPONSE},
        {"item": "primary_predictor", "value": PREDICTOR},
        {"item": "measurement_caveat", "value": "Resampling quantiles approximate pair-level observation dispersion; they are not treated as independent biological replicates."},
        {"item": "selection_guardrail", "value": "The original Phase 14 primary composite remains primary; extended predictors are not substituted post hoc."},
    ]
    write_csv(outdir / "15A_input_audit.csv", audit)
    payload = {
        "phase": "15A", "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": status, "pairs": len(prepared), "cells": len(cells),
        "observation_sd_floor": adaptive_floor, "pairs_using_floor": floor_count,
        "source": str(source),
    }
    write_json(outdir / "15A_run_status.json", payload)
    readme = f"""PHASE 15A — BAYESIAN H3 INPUT AUDIT
=====================================
Status: {status}
Pairs retained: {len(prepared)}
Cells represented: {len(cells)}
Approximate observation-SD floor: {adaptive_floor:.6f}
Pairs using floor: {floor_count}

The primary predictor remains the original Phase 14 composite stress change. This
prevents choosing a different environmental metric after seeing the Phase 14 results.
The Phase 14 equal-event resampling interval is propagated as approximate observation
dispersion in Phase 15; it does not create additional biological replication.
"""
    (outdir / "15A_README.txt").write_text(readme, encoding="utf-8")
    print(readme)
    print(f"OUTPUT={outdir / '15A_bayesian_model_input.csv'}")
    if status != "PASS":
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
