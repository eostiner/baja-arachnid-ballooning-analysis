#!/usr/bin/env python3
"""Phase 14H — temporal H3 sensitivity to C1/C2/C3 ballooning evidence thresholds."""
from __future__ import annotations

import argparse
import random
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from phase14_common import default_analysis_output_root, load_periods, quantile, read_delimited, write_csv, write_json
from phase14_model_utils import as_float, classify, run_model

SCRIPT_VERSION = "14H_v0.3.0_2026-07-18"
PRIMARY_RESPONSE = "resampled_delta_simpson_C3_minus_N0_median"
PRIMARY_PREDICTOR = "delta_mean_extended_stress_composite_z_period2_minus_period1"
THRESHOLDS = {
    "C1_D1_only": {"D1"},
    "C2_D1_D2": {"D1", "D2"},
    "C3_D1_D3_primary": {"D1", "D2", "D3"},
}


def beta_metrics(left: set[str], right: set[str]) -> tuple[float, float]:
    shared = len(left & right)
    left_only = len(left - right)
    right_only = len(right - left)
    union = shared + left_only + right_only
    jaccard = (left_only + right_only) / union if union else 0.0
    minimum = min(left_only, right_only)
    simpson = minimum / (shared + minimum) if (shared + minimum) else 0.0
    return jaccard, simpson


def main() -> int:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Test temporal H3 sensitivity across C1/C2/C3 definitions.")
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--record-index", type=Path)
    parser.add_argument("--stress-join", type=Path)
    parser.add_argument("--period-config", type=Path, default=here / "configs" / "phase_14_temporal_windows_frozen.csv")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument("--bootstrap-iterations", type=int, default=2500)
    parser.add_argument("--permutation-iterations", type=int, default=10000)
    parser.add_argument("--min-records", type=int, default=5)
    parser.add_argument("--min-events", type=int, default=3)
    parser.add_argument("--min-genera", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20260718)
    args = parser.parse_args()

    root = args.project_root.expanduser().resolve()
    base = default_analysis_output_root(root)
    record_path = args.record_index.expanduser().resolve() if args.record_index else base / "14A_temporal_feasibility_audit" / "14A_record_temporal_index.tsv"
    stress_path = args.stress_join.expanduser().resolve() if args.stress_join else base / "14E_extended_stress_join" / "14E2_H3_extended_model_input_all_pairs.csv"
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else base / "14H_trait_threshold_temporal_sensitivity"
    outdir.mkdir(parents=True, exist_ok=True)
    if not record_path.exists() or not stress_path.exists():
        raise FileNotFoundError("Required Phase 14A or 14E2 input is missing.")

    _, records = read_delimited(record_path)
    _, stress_rows = read_delimited(stress_path)
    periods_list = load_periods(args.period_config.expanduser().resolve())
    periods = {period["period_id"]: period for period in periods_list}
    period_ids = [period["period_id"] for period in periods_list]
    stress_lookup = {
        (str(row["grid_cell_id"]), str(row["period_1"]), str(row["period_2"])): row
        for row in stress_rows
    }
    rng = random.Random(args.seed)

    coverage_rows: list[dict[str, Any]] = []
    turnover_rows: list[dict[str, Any]] = []
    model_rows: list[dict[str, Any]] = []

    for threshold_index, (threshold, allowed) in enumerate(THRESHOLDS.items(), start=1):
        event_genera: dict[tuple[str, str, str], dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
        metadata: dict[tuple[str, str], str] = {}
        record_counts: dict[tuple[str, str, str], int] = defaultdict(int)
        for record in records:
            trait = str(record.get("trait_group", ""))
            evidence = str(record.get("evidence_class", "")).upper()
            group = "N0" if trait == "N0" else ("B" if evidence in allowed else "")
            if not group:
                continue
            key = (str(record["grid_cell_id"]), str(record["period_id"]), group)
            event_genera[key][str(record["event_key"])].add(str(record["genus"]))
            record_counts[key] += 1
            metadata[(str(record["grid_cell_id"]), str(record["period_id"]))] = str(record.get("latitude_band", ""))

        cell_period_eligible: dict[tuple[str, str], bool] = {}
        all_cells = sorted({key[0] for key in event_genera})
        for cell in all_cells:
            for period in period_ids:
                stats = {}
                for group in ("B", "N0"):
                    events = event_genera[(cell, period, group)]
                    genera = set().union(*events.values()) if events else set()
                    stats[group] = {
                        "records": record_counts[(cell, period, group)],
                        "events": len(events),
                        "genera": len(genera),
                    }
                eligible = all(
                    stats[group]["records"] >= args.min_records
                    and stats[group]["events"] >= args.min_events
                    and stats[group]["genera"] >= args.min_genera
                    for group in ("B", "N0")
                )
                cell_period_eligible[(cell, period)] = eligible
                coverage_rows.append({
                    "threshold": threshold,
                    "grid_cell_id": cell,
                    "period_id": period,
                    "latitude_band": metadata.get((cell, period), ""),
                    "ballooning_records": stats["B"]["records"],
                    "ballooning_events": stats["B"]["events"],
                    "ballooning_genera": stats["B"]["genera"],
                    "n0_records": stats["N0"]["records"],
                    "n0_events": stats["N0"]["events"],
                    "n0_genera": stats["N0"]["genera"],
                    "paired_eligible": int(eligible),
                })

        threshold_pairs: list[dict[str, Any]] = []
        pair_counter = 0
        for cell in all_cells:
            for p1, p2 in zip(period_ids, period_ids[1:]):
                if not (cell_period_eligible.get((cell, p1), False) and cell_period_eligible.get((cell, p2), False)):
                    continue
                pair_counter += 1
                samples = {(period, group): event_genera[(cell, period, group)] for period in (p1, p2) for group in ("B", "N0")}
                common_events = min(len(value) for value in samples.values())
                if common_events < 1:
                    continue
                delta_simpson: list[float] = []
                delta_jaccard: list[float] = []
                for _ in range(args.iterations):
                    metrics = {}
                    for group in ("B", "N0"):
                        period_sets = []
                        for period in (p1, p2):
                            event_map = samples[(period, group)]
                            selected = rng.sample(list(event_map), common_events)
                            period_sets.append(set().union(*(event_map[event] for event in selected)))
                        metrics[group] = beta_metrics(period_sets[0], period_sets[1])
                    delta_jaccard.append(metrics["B"][0] - metrics["N0"][0])
                    delta_simpson.append(metrics["B"][1] - metrics["N0"][1])
                stress = stress_lookup.get((cell, p1, p2), {})
                row = {
                    "threshold": threshold,
                    "pair_id": f"{threshold}_TPAIR_{pair_counter:04d}",
                    "grid_cell_id": cell,
                    "period_1": p1,
                    "period_2": p2,
                    "latitude_band": metadata.get((cell, p1), metadata.get((cell, p2), "")),
                    "common_event_resample_n": common_events,
                    "resampled_delta_simpson_C3_minus_N0_median": quantile(delta_simpson, 0.5),
                    "resampled_delta_simpson_C3_minus_N0_q025": quantile(delta_simpson, 0.025),
                    "resampled_delta_simpson_C3_minus_N0_q975": quantile(delta_simpson, 0.975),
                    "resampled_delta_jaccard_C3_minus_N0_median": quantile(delta_jaccard, 0.5),
                    PRIMARY_PREDICTOR: stress.get(PRIMARY_PREDICTOR, ""),
                }
                turnover_rows.append(row)
                threshold_pairs.append(row)

        complete = [row for row in threshold_pairs if as_float(row.get(PRIMARY_PREDICTOR)) is not None]
        cells = {row["grid_cell_id"] for row in complete}
        if len(complete) < 6 or len(cells) < 5:
            model_rows.append({
                "threshold": threshold,
                "status": "SKIPPED_INSUFFICIENT_PAIRS",
                "n_pairs": len(complete),
                "n_cells": len(cells),
            })
            continue
        try:
            model = run_model(
                complete, PRIMARY_PREDICTOR, PRIMARY_RESPONSE, periods,
                args.bootstrap_iterations, args.permutation_iterations,
                args.seed + threshold_index * 5000,
                include_band=False, include_time=True, include_effort=True,
            )
            model_rows.append({
                "threshold": threshold,
                "status": "OK",
                "n_pairs": len(model["built"]["retained"]),
                "n_cells": len({item["cell"] for item in model["built"]["retained"]}),
                "coefficient_per_predictor_sd": model["beta"],
                "cluster_bootstrap_q025": model["ci_low"],
                "cluster_bootstrap_q975": model["ci_high"],
                "wild_cluster_p_two_sided": model["pvalue"],
                "conclusion": classify(model["beta"], model["ci_low"], model["ci_high"], model["pvalue"]),
            })
        except Exception as exc:
            model_rows.append({
                "threshold": threshold,
                "status": f"MODEL_FAILED: {exc}",
                "n_pairs": len(complete),
                "n_cells": len(cells),
            })

    write_csv(outdir / "14H_trait_threshold_cell_period_coverage.csv", coverage_rows)
    write_csv(outdir / "14H_trait_threshold_turnover.csv", turnover_rows)
    write_csv(outdir / "14H_trait_threshold_models.csv", model_rows)
    readme = f"""PHASE 14H — TRAIT-THRESHOLD TEMPORAL SENSITIVITY
====================================================
Thresholds tested: C1 (D1 only), C2 (D1+D2), and C3 (D1+D2+D3 primary)
Turnover rows produced across thresholds: {len(turnover_rows)}
Models fit successfully: {sum(row.get('status') == 'OK' for row in model_rows)}

C4 is not reconstructed here because D4 records were intentionally excluded from
the Phase 14A primary temporal index. The C1/C2 analyses ask whether any H3 signal
is restricted to the broader juvenile-ballooning definition or persists under
more conservative evidence thresholds.
"""
    (outdir / "14H_README.txt").write_text(readme, encoding="utf-8")
    write_json(outdir / "14H_run_status.json", {
        "phase": "14H",
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETED",
        "thresholds": list(THRESHOLDS),
        "models": model_rows,
    })
    print(readme)
    print(f"OUTPUT_DIR={outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
