#!/usr/bin/env python3
"""Phase 14I — synthesize all available non-field H3 evidence without cherry-picking."""
from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from phase14_common import default_analysis_output_root, read_delimited, write_csv, write_json

SCRIPT_VERSION = "14I_v0.3.0_2026-07-18"


def num(value: Any) -> float | None:
    try:
        result = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def read_rows(path: Path) -> list[dict[str, str]]:
    return read_delimited(path)[1] if path.exists() else []


def evidence_row(source: str, model: str, row: dict[str, Any]) -> dict[str, Any] | None:
    beta = num(row.get("coefficient_per_predictor_sd", row.get("coefficient")))
    low = num(row.get("cluster_bootstrap_q025"))
    high = num(row.get("cluster_bootstrap_q975"))
    pvalue = num(row.get("wild_cluster_p_two_sided"))
    if beta is None:
        return None
    return {
        "source": source,
        "model": model,
        "status": row.get("status", "OK"),
        "n_pairs": row.get("n_pairs", ""),
        "n_cells": row.get("n_cells", ""),
        "coefficient_per_predictor_sd": beta,
        "cluster_bootstrap_q025": low if low is not None else "",
        "cluster_bootstrap_q975": high if high is not None else "",
        "wild_cluster_p_two_sided": pvalue if pvalue is not None else "",
        "positive_direction": int(beta > 0),
        "interval_excludes_zero_positive": int(low is not None and low > 0),
        "conclusion": row.get("conclusion", ""),
    }


def write_summary_plot(outdir: Path, rows: list[dict[str, Any]]) -> list[str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    plotted = [row for row in rows if num(row.get("cluster_bootstrap_q025")) is not None and num(row.get("cluster_bootstrap_q975")) is not None]
    if not plotted:
        return []
    plotted = list(reversed(plotted))
    y = list(range(len(plotted)))
    estimates = [float(row["coefficient_per_predictor_sd"]) for row in plotted]
    lower = [float(row["cluster_bootstrap_q025"]) for row in plotted]
    upper = [float(row["cluster_bootstrap_q975"]) for row in plotted]
    labels = [f"{row['source']}: {row['model']}" for row in plotted]
    fig, ax = plt.subplots(figsize=(10.5, max(5.5, 0.42 * len(plotted) + 1.8)))
    ax.errorbar(estimates, y, xerr=[
        [estimate - lo for estimate, lo in zip(estimates, lower)],
        [hi - estimate for estimate, hi in zip(estimates, upper)],
    ], fmt="o", capsize=3)
    ax.axvline(0, linewidth=1)
    ax.set_yticks(y, labels)
    ax.set_xlabel("C3 − N0 turnover response per 1 SD worsening stress")
    ax.set_title("Phase 14 H3 evidence synthesis")
    fig.tight_layout()
    outputs = []
    for suffix in ("png", "svg"):
        path = outdir / f"14I_H3_evidence_synthesis.{suffix}"
        fig.savefig(path, dpi=300 if suffix == "png" else None)
        outputs.append(str(path))
    plt.close(fig)
    return outputs


def main() -> int:
    parser = argparse.ArgumentParser(description="Synthesize Phase 14 H3 evidence.")
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    root = args.project_root.expanduser().resolve()
    base = default_analysis_output_root(root)
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else base / "14I_H3_evidence_synthesis"
    outdir.mkdir(parents=True, exist_ok=True)

    evidence: list[dict[str, Any]] = []
    d_rows = read_rows(base / "14D_recent_stress_tracking_test" / "14D_primary_model_summary.csv")
    if d_rows:
        row = evidence_row("14D", "original composite primary", d_rows[0])
        if row:
            evidence.append(row)

    f_rows = read_rows(base / "14F_extended_stress_sensitivity" / "14F_extended_model_results.csv")
    for item in f_rows:
        if item.get("status") != "OK":
            continue
        if item.get("model_id") in {
            "E_PRIMARY_COMPOSITE_MEAN",
            "E_SECONDARY_THERMAL_MEAN",
            "E_SECONDARY_VPD_MEAN",
            "E_SECONDARY_MOISTURE_MEAN",
            "E_SECONDARY_VEGETATION_MEAN",
        }:
            row = evidence_row("14F", str(item.get("model_id", "")), item)
            if row:
                evidence.append(row)

    g_rows = read_rows(base / "14G_sampling_continuity_sensitivity" / "14G_sampling_continuity_models.csv")
    for item in g_rows:
        if item.get("status") == "OK":
            row = evidence_row("14G", str(item.get("scheme", "sampling continuity")), item)
            if row:
                evidence.append(row)

    h_rows = read_rows(base / "14H_trait_threshold_temporal_sensitivity" / "14H_trait_threshold_models.csv")
    for item in h_rows:
        if item.get("status") == "OK":
            row = evidence_row("14H", str(item.get("threshold", "trait threshold")), item)
            if row:
                evidence.append(row)

    write_csv(outdir / "14I_H3_evidence_table.csv", evidence)
    positive = [row for row in evidence if int(row["positive_direction"]) == 1]
    positive_excluding_zero = [row for row in evidence if int(row["interval_excludes_zero_positive"]) == 1]
    original = next((row for row in evidence if row["source"] == "14D"), None)
    expanded = next((row for row in evidence if row["source"] == "14F" and row["model"] == "E_PRIMARY_COMPOSITE_MEAN"), None)

    if original and expanded and int(original["interval_excludes_zero_positive"]) and int(expanded["interval_excludes_zero_positive"]):
        overall = "REPEATED_EXPLORATORY_DIRECTIONAL_SUPPORT_NOT_CAUSAL"
        wording = (
            "Both the original and expanded exploratory models estimated stronger C3 temporal turnover under worsening stress, "
            "with uncertainty intervals above zero. Because the analysis remains based on few opportunistically sampled cells, "
            "this is repeated exploratory support rather than confirmatory or causal evidence."
        )
    elif original and expanded and float(original["coefficient_per_predictor_sd"]) > 0 and float(expanded["coefficient_per_predictor_sd"]) > 0:
        overall = "POSITIVE_BUT_UNCERTAIN_ACROSS_PRIMARY_MODELS"
        wording = (
            "Both primary analyses pointed in the H3-predicted direction, but at least one uncertainty interval included zero. "
            "The expanded environmental data therefore preserve a plausible signal without demonstrating that it is real."
        )
    elif positive and len(positive) >= max(2, len(evidence) // 2):
        overall = "MIXED_WITH_MORE_POSITIVE_THAN_NEGATIVE_ESTIMATES"
        wording = (
            "More analyses pointed in the predicted direction than in the opposite direction, but the pattern was not sufficiently "
            "consistent or precise to support H3."
        )
    else:
        overall = "NO_CONSISTENT_DIRECTIONAL_SUPPORT"
        wording = (
            "The expanded analyses did not produce a consistent directional relationship between worsening recent stress and stronger C3 turnover."
        )

    summary_rows = [
        {"item": "overall_status", "value": overall},
        {"item": "models_with_numeric_estimates", "value": len(evidence)},
        {"item": "positive_estimates", "value": len(positive)},
        {"item": "positive_intervals_excluding_zero", "value": len(positive_excluding_zero)},
        {"item": "manuscript_wording", "value": wording},
        {"item": "scope_guardrail", "value": "Inference is restricted to repeated eligible cells and is not peninsula-wide."},
        {"item": "causal_guardrail", "value": "Observational turnover-stress association cannot demonstrate stress-caused redistribution."},
        {"item": "power_guardrail", "value": "Additional environmental products improve measurement but do not add biological replicates."},
    ]
    write_csv(outdir / "14I_conclusion_summary.csv", summary_rows)
    plot_outputs = write_summary_plot(outdir, evidence)
    readme = f"""PHASE 14I — H3 EVIDENCE SYNTHESIS
===================================
Overall status: {overall}
Models with numeric estimates: {len(evidence)}
Positive estimates: {len(positive)}
Positive estimates whose 95% interval excluded zero: {len(positive_excluding_zero)}

Interpretation:
{wording}

This synthesis is the result to use when deciding manuscript language. It prevents
selection of whichever individual environmental variable happened to look strongest.
All Phase 14 analyses remain exploratory because the biological dataset contains few
repeated cell-period comparisons and was assembled opportunistically.
"""
    (outdir / "14I_README.txt").write_text(readme, encoding="utf-8")
    write_json(outdir / "14I_run_status.json", {
        "phase": "14I",
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETED",
        "overall_status": overall,
        "models_synthesized": len(evidence),
        "plot_outputs": plot_outputs,
    })
    print(readme)
    print(f"OUTPUT_DIR={outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
