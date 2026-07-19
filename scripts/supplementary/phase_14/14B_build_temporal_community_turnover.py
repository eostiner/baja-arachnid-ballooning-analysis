#!/usr/bin/env python3
"""
Phase 14B — construct paired within-cell temporal C3/N0 community turnover.

Requires the record-level temporal index and eligibility table from Phase 14A.
For each eligible adjacent period comparison, it computes raw and equal-event
resampled Jaccard dissimilarity and Simpson replacement for C3 and N0, then writes
the paired C3-minus-N0 response needed for an H3 stress-anomaly model.
"""
from __future__ import annotations

import argparse
import csv
import random
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from phase14_common import default_analysis_output_root, quantile, read_delimited, write_csv, write_json

SCRIPT_VERSION = "14B_v0.1.0_2026-07-18"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build temporal C3/N0 turnover response table.")
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--audit-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260718)
    return parser.parse_args()


def beta_metrics(left: set[str], right: set[str]) -> dict[str, float | int]:
    shared = len(left & right)
    left_only = len(left - right)
    right_only = len(right - left)
    union = shared + left_only + right_only
    jaccard = (left_only + right_only) / union if union else 0.0
    minimum = min(left_only, right_only)
    simpson = minimum / (shared + minimum) if (shared + minimum) else 0.0
    return {
        "shared_genera": shared,
        "period1_only_genera": left_only,
        "period2_only_genera": right_only,
        "jaccard_dissimilarity": jaccard,
        "simpson_replacement": simpson,
    }


def summarize_draws(values: list[float], prefix: str) -> dict[str, float]:
    return {
        f"{prefix}_mean": sum(values) / len(values) if values else float("nan"),
        f"{prefix}_median": quantile(values, 0.5),
        f"{prefix}_q025": quantile(values, 0.025),
        f"{prefix}_q975": quantile(values, 0.975),
    }


def main() -> int:
    args = parse_args()
    if args.iterations < 1:
        raise ValueError("--iterations must be at least 1")
    root = args.project_root.expanduser().resolve()
    audit_dir = args.audit_dir.expanduser().resolve() if args.audit_dir else (
        default_analysis_output_root(root) / "14A_temporal_feasibility_audit"
    )
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else (
        default_analysis_output_root(root) / "14B_temporal_community_turnover"
    )
    outdir.mkdir(parents=True, exist_ok=True)
    record_path = audit_dir / "14A_record_temporal_index.tsv"
    eligibility_path = audit_dir / "14A_adjacent_period_pair_eligibility.csv"
    if not record_path.exists() or not eligibility_path.exists():
        raise FileNotFoundError("Phase 14A outputs not found. Run 14A before 14B.")

    _, records = read_delimited(record_path)
    _, eligibility_rows = read_delimited(eligibility_path)
    eligible = [row for row in eligibility_rows if int(row["temporal_pair_eligible"]) == 1]
    if not eligible:
        status = {
            "phase": "14B",
            "script_version": SCRIPT_VERSION,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "status": "NOT_RUN_NO_ELIGIBLE_TEMPORAL_PAIRS",
            "message": "Phase 14A found no adjacent cell-period pairs meeting the frozen paired C3/N0 thresholds.",
        }
        write_json(outdir / "14B_run_status.json", status)
        (outdir / "14B_README.txt").write_text(status["message"] + "\n", encoding="utf-8")
        print(status["message"])
        return 0

    # event_genera[(cell, period, trait)][event] = set(genera)
    event_genera: dict[tuple[str, str, str], dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    metadata: dict[tuple[str, str], dict[str, str]] = {}
    for row in records:
        key = (row["grid_cell_id"], row["period_id"], row["trait_group"])
        event_genera[key][row["event_key"]].add(row["genus"])
        metadata[(row["grid_cell_id"], row["period_id"])] = {
            "latitude_band": row["latitude_band"],
            "decimalLatitude": row["decimalLatitude"],
            "decimalLongitude": row["decimalLongitude"],
        }

    rng = random.Random(args.seed)
    output_rows: list[dict[str, Any]] = []
    for pair_number, pair in enumerate(eligible, start=1):
        cell = pair["grid_cell_id"]
        p1 = pair["period_1"]
        p2 = pair["period_2"]
        samples = {
            (period, trait): event_genera[(cell, period, trait)]
            for period in (p1, p2)
            for trait in ("C3", "N0")
        }
        event_counts = {key: len(value) for key, value in samples.items()}
        common_events = min(event_counts.values())
        if common_events < 1:
            continue

        raw = {}
        for trait in ("C3", "N0"):
            left = set().union(*samples[(p1, trait)].values()) if samples[(p1, trait)] else set()
            right = set().union(*samples[(p2, trait)].values()) if samples[(p2, trait)] else set()
            metrics = beta_metrics(left, right)
            for name, value in metrics.items():
                raw[f"raw_{trait.lower()}_{name}"] = value
        raw["raw_delta_jaccard_C3_minus_N0"] = float(raw["raw_c3_jaccard_dissimilarity"]) - float(raw["raw_n0_jaccard_dissimilarity"])
        raw["raw_delta_simpson_C3_minus_N0"] = float(raw["raw_c3_simpson_replacement"]) - float(raw["raw_n0_simpson_replacement"])

        draws: dict[str, list[float]] = defaultdict(list)
        for _ in range(args.iterations):
            draw_metrics = {}
            for trait in ("C3", "N0"):
                period_sets = []
                for period in (p1, p2):
                    event_map = samples[(period, trait)]
                    selected = rng.sample(list(event_map), common_events)
                    genera = set().union(*(event_map[event] for event in selected))
                    period_sets.append(genera)
                metrics = beta_metrics(period_sets[0], period_sets[1])
                draw_metrics[trait] = metrics
                draws[f"{trait.lower()}_jaccard"].append(float(metrics["jaccard_dissimilarity"]))
                draws[f"{trait.lower()}_simpson"].append(float(metrics["simpson_replacement"]))
            draws["delta_jaccard_C3_minus_N0"].append(
                float(draw_metrics["C3"]["jaccard_dissimilarity"]) - float(draw_metrics["N0"]["jaccard_dissimilarity"])
            )
            draws["delta_simpson_C3_minus_N0"].append(
                float(draw_metrics["C3"]["simpson_replacement"]) - float(draw_metrics["N0"]["simpson_replacement"])
            )

        representative = metadata.get((cell, p1), metadata.get((cell, p2), {}))
        row: dict[str, Any] = {
            "pair_id": f"TPAIR_{pair_number:04d}",
            "grid_cell_id": cell,
            "period_1": p1,
            "period_2": p2,
            "latitude_band": pair.get("latitude_band", representative.get("latitude_band", "")),
            "common_event_resample_n": common_events,
            "c3_period1_events": event_counts[(p1, "C3")],
            "c3_period2_events": event_counts[(p2, "C3")],
            "n0_period1_events": event_counts[(p1, "N0")],
            "n0_period2_events": event_counts[(p2, "N0")],
            **raw,
        }
        for key, values in draws.items():
            row.update(summarize_draws(values, f"resampled_{key}"))
        output_rows.append(row)

    output_path = outdir / "14B_paired_temporal_turnover.csv"
    write_csv(output_path, output_rows)
    status_name = "COMPLETED" if output_rows else "NO_VALID_ROWS_AFTER_EVENT_RESAMPLING"
    status = {
        "phase": "14B",
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": status_name,
        "iterations": args.iterations,
        "seed": args.seed,
        "eligible_pairs_from_14A": len(eligible),
        "turnover_rows_written": len(output_rows),
        "primary_future_H3_response": "resampled_delta_simpson_C3_minus_N0_median",
        "secondary_future_H3_response": "resampled_delta_jaccard_C3_minus_N0_median",
        "interpretation": (
            "A positive stress coefficient for the paired C3-minus-N0 response would support H3. "
            "This table alone contains temporal turnover, not environmental anomalies."
        ),
    }
    write_json(outdir / "14B_run_status.json", status)
    readme = f"""PHASE 14B — PAIRED TEMPORAL COMMUNITY TURNOVER
================================================
Status: {status_name}
Eligible Phase 14A temporal pairs: {len(eligible):,}
Rows written: {len(output_rows):,}
Equal-event resampling iterations: {args.iterations:,}
Random seed: {args.seed}

Primary future H3 response:
  resampled_delta_simpson_C3_minus_N0_median

Secondary future H3 response:
  resampled_delta_jaccard_C3_minus_N0_median

The environmental stress-anomaly predictors have not yet been joined. A positive
association between stress change and the C3-minus-N0 response is the directional
prediction for H3.
"""
    (outdir / "14B_README.txt").write_text(readme, encoding="utf-8")
    print(readme)
    print(f"OUTPUT={output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
