#!/usr/bin/env python3
"""Phase 14G — sampling-continuity sensitivity for opportunistic temporal records.

Rebuilds C3/N0 temporal turnover after restricting records to datasets and/or
calendar quarters represented in both periods and both trait groups.
"""
from __future__ import annotations

import argparse
import random
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from phase14_common import default_analysis_output_root, load_periods, quantile, read_delimited, write_csv, write_json
from phase14_model_utils import as_float, classify, run_model

SCRIPT_VERSION = "14G_v0.3.0_2026-07-18"
PRIMARY_RESPONSE = "resampled_delta_simpson_C3_minus_N0_median"
PRIMARY_PREDICTOR = "delta_mean_extended_stress_composite_z_period2_minus_period1"


def beta_metrics(left: set[str], right: set[str]) -> tuple[float, float]:
    shared = len(left & right)
    left_only = len(left - right)
    right_only = len(right - left)
    union = shared + left_only + right_only
    jaccard = (left_only + right_only) / union if union else 0.0
    minimum = min(left_only, right_only)
    simpson = minimum / (shared + minimum) if (shared + minimum) else 0.0
    return jaccard, simpson


def q(values: list[float], probability: float) -> float:
    return quantile(values, probability)


def main() -> int:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Test dataset and seasonal sampling-continuity sensitivities.")
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--record-index", type=Path)
    parser.add_argument("--eligibility-file", type=Path)
    parser.add_argument("--stress-join", type=Path)
    parser.add_argument("--period-config", type=Path, default=here / "configs" / "phase_14_temporal_windows_frozen.csv")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument("--bootstrap-iterations", type=int, default=2500)
    parser.add_argument("--permutation-iterations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260718)
    args = parser.parse_args()

    root = args.project_root.expanduser().resolve()
    base = default_analysis_output_root(root)
    audit_dir = base / "14A_temporal_feasibility_audit"
    record_path = args.record_index.expanduser().resolve() if args.record_index else audit_dir / "14A_record_temporal_index.tsv"
    eligibility_path = args.eligibility_file.expanduser().resolve() if args.eligibility_file else audit_dir / "14A_adjacent_period_pair_eligibility.csv"
    stress_path = args.stress_join.expanduser().resolve() if args.stress_join else base / "14E_extended_stress_join" / "14E2_H3_extended_model_input_all_pairs.csv"
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else base / "14G_sampling_continuity_sensitivity"
    outdir.mkdir(parents=True, exist_ok=True)

    for path in (record_path, eligibility_path, stress_path):
        if not path.exists():
            raise FileNotFoundError(path)
    _, records = read_delimited(record_path)
    _, eligibility = read_delimited(eligibility_path)
    _, stress_rows = read_delimited(stress_path)
    periods = {period["period_id"]: period for period in load_periods(args.period_config.expanduser().resolve())}
    eligible = [row for row in eligibility if int(row["temporal_pair_eligible"]) == 1]
    stress_lookup = {
        (str(row["grid_cell_id"]), str(row["period_1"]), str(row["period_2"])): row
        for row in stress_rows
    }

    # event_info[(cell, period, trait)][event] = metadata and genera
    event_info: dict[tuple[str, str, str], dict[str, dict[str, Any]]] = defaultdict(dict)
    for record in records:
        key = (str(record["grid_cell_id"]), str(record["period_id"]), str(record["trait_group"]))
        event = str(record["event_key"])
        month_text = str(record.get("month", "")).strip()
        quarter = ""
        if month_text:
            try:
                month = int(float(month_text))
                quarter = f"Q{(month - 1) // 3 + 1}" if 1 <= month <= 12 else ""
            except ValueError:
                quarter = ""
        info = event_info[key].setdefault(event, {
            "genera": set(),
            "dataset": str(record.get("dataset_key", "UNSPECIFIED_DATASET")),
            "quarter": quarter,
        })
        info["genera"].add(str(record["genus"]))

    schemes = (
        "dataset_continuity_strict",
        "quarter_continuity_strict",
        "dataset_and_quarter_continuity_strict",
    )
    rng = random.Random(args.seed)
    turnover_rows: list[dict[str, Any]] = []
    coverage_rows: list[dict[str, Any]] = []

    for pair_number, pair in enumerate(eligible, start=1):
        cell, p1, p2 = str(pair["grid_cell_id"]), str(pair["period_1"]), str(pair["period_2"])
        maps = {(period, trait): event_info[(cell, period, trait)] for period in (p1, p2) for trait in ("C3", "N0")}
        dataset_sets = {
            key: {info["dataset"] for info in events.values() if info["dataset"]}
            for key, events in maps.items()
        }
        quarter_sets = {
            key: {info["quarter"] for info in events.values() if info["quarter"]}
            for key, events in maps.items()
        }
        shared_datasets = set.intersection(*(dataset_sets[key] for key in dataset_sets)) if dataset_sets else set()
        shared_quarters = set.intersection(*(quarter_sets[key] for key in quarter_sets)) if quarter_sets else set()

        for scheme in schemes:
            filtered: dict[tuple[str, str], dict[str, dict[str, Any]]] = {}
            for key, events in maps.items():
                selected = {}
                for event, info in events.items():
                    dataset_ok = info["dataset"] in shared_datasets
                    quarter_ok = info["quarter"] in shared_quarters and bool(info["quarter"])
                    keep = (
                        dataset_ok if scheme == "dataset_continuity_strict" else
                        quarter_ok if scheme == "quarter_continuity_strict" else
                        dataset_ok and quarter_ok
                    )
                    if keep:
                        selected[event] = info
                filtered[key] = selected
            event_counts = {key: len(events) for key, events in filtered.items()}
            common_events = min(event_counts.values()) if event_counts else 0
            coverage_rows.append({
                "pair_id": f"TPAIR_{pair_number:04d}",
                "grid_cell_id": cell,
                "period_1": p1,
                "period_2": p2,
                "scheme": scheme,
                "shared_dataset_count": len(shared_datasets),
                "shared_quarter_count": len(shared_quarters),
                "common_event_resample_n": common_events,
                "usable": int(common_events >= 1),
            })
            if common_events < 1:
                continue

            delta_simpson: list[float] = []
            delta_jaccard: list[float] = []
            for _ in range(args.iterations):
                metrics: dict[str, tuple[float, float]] = {}
                for trait in ("C3", "N0"):
                    period_genera = []
                    for period in (p1, p2):
                        events = filtered[(period, trait)]
                        selected_events = rng.sample(list(events), common_events)
                        genera = set().union(*(events[event]["genera"] for event in selected_events))
                        period_genera.append(genera)
                    metrics[trait] = beta_metrics(period_genera[0], period_genera[1])
                delta_jaccard.append(metrics["C3"][0] - metrics["N0"][0])
                delta_simpson.append(metrics["C3"][1] - metrics["N0"][1])

            row: dict[str, Any] = {
                "pair_id": f"TPAIR_{pair_number:04d}",
                "grid_cell_id": cell,
                "period_1": p1,
                "period_2": p2,
                "latitude_band": pair.get("latitude_band", ""),
                "scheme": scheme,
                "common_event_resample_n": common_events,
                "shared_dataset_count": len(shared_datasets),
                "shared_quarter_count": len(shared_quarters),
                "resampled_delta_simpson_C3_minus_N0_median": q(delta_simpson, 0.5),
                "resampled_delta_simpson_C3_minus_N0_q025": q(delta_simpson, 0.025),
                "resampled_delta_simpson_C3_minus_N0_q975": q(delta_simpson, 0.975),
                "resampled_delta_jaccard_C3_minus_N0_median": q(delta_jaccard, 0.5),
            }
            stress = stress_lookup.get((cell, p1, p2), {})
            row[PRIMARY_PREDICTOR] = stress.get(PRIMARY_PREDICTOR, "")
            turnover_rows.append(row)

    write_csv(outdir / "14G_sampling_continuity_turnover.csv", turnover_rows)
    write_csv(outdir / "14G_sampling_continuity_coverage.csv", coverage_rows)

    model_rows: list[dict[str, Any]] = []
    for index, scheme in enumerate(schemes, start=1):
        subset = [row for row in turnover_rows if row["scheme"] == scheme and as_float(row.get(PRIMARY_PREDICTOR)) is not None]
        cells = {row["grid_cell_id"] for row in subset}
        if len(subset) < 6 or len(cells) < 5:
            model_rows.append({
                "scheme": scheme,
                "status": "SKIPPED_INSUFFICIENT_PAIRS",
                "n_pairs": len(subset),
                "n_cells": len(cells),
            })
            continue
        fitted = None
        control_set = "effort_time"
        error = ""
        try:
            fitted = run_model(
                subset, PRIMARY_PREDICTOR, PRIMARY_RESPONSE, periods,
                args.bootstrap_iterations, args.permutation_iterations, args.seed + index * 1000,
                include_band=False, include_time=True, include_effort=True,
            )
        except Exception as exc:
            error = str(exc)
            control_set = "effort_only"
            try:
                fitted = run_model(
                    subset, PRIMARY_PREDICTOR, PRIMARY_RESPONSE, periods,
                    args.bootstrap_iterations, args.permutation_iterations, args.seed + index * 1000 + 100,
                    include_band=False, include_time=False, include_effort=True,
                )
            except Exception as exc2:
                error += f"; fallback failed: {exc2}"
        if fitted is None:
            model_rows.append({
                "scheme": scheme,
                "status": f"MODEL_FAILED: {error}",
                "n_pairs": len(subset),
                "n_cells": len(cells),
            })
            continue
        conclusion = classify(fitted["beta"], fitted["ci_low"], fitted["ci_high"], fitted["pvalue"])
        model_rows.append({
            "scheme": scheme,
            "status": "OK",
            "control_set": control_set,
            "n_pairs": len(fitted["built"]["retained"]),
            "n_cells": len({item["cell"] for item in fitted["built"]["retained"]}),
            "coefficient_per_predictor_sd": fitted["beta"],
            "cluster_bootstrap_q025": fitted["ci_low"],
            "cluster_bootstrap_q975": fitted["ci_high"],
            "wild_cluster_p_two_sided": fitted["pvalue"],
            "conclusion": conclusion,
        })
    write_csv(outdir / "14G_sampling_continuity_models.csv", model_rows)

    readme = f"""PHASE 14G — SAMPLING-CONTINUITY SENSITIVITY
================================================
Eligible original temporal pairs: {len(eligible)}
Continuity turnover rows produced: {len(turnover_rows)}
Dataset-continuity usable pairs: {sum(row['scheme'] == 'dataset_continuity_strict' for row in turnover_rows)}
Quarter-continuity usable pairs: {sum(row['scheme'] == 'quarter_continuity_strict' for row in turnover_rows)}
Dataset+quarter usable pairs: {sum(row['scheme'] == 'dataset_and_quarter_continuity_strict' for row in turnover_rows)}

These restrictions test whether apparent temporal turnover could be caused by changes
in contributing datasets or seasonal collecting. They generally reduce sample size,
so a failed or uncertain model is not evidence that sampling bias is absent.
"""
    (outdir / "14G_README.txt").write_text(readme, encoding="utf-8")
    write_json(outdir / "14G_run_status.json", {
        "phase": "14G",
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETED",
        "original_eligible_pairs": len(eligible),
        "turnover_rows": len(turnover_rows),
        "models": model_rows,
    })
    print(readme)
    print(f"OUTPUT_DIR={outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
