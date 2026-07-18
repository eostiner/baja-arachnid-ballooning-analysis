#!/usr/bin/env python3
"""
Phase 13A — conservative input/provenance audit.

Purpose:
- Discover candidate retained Phase 10/11/12 inputs beneath a frozen project root.
- Inventory likely C3/N0 incidence, cell-coordinate, environmental, boundary/ecoregion,
  and trait files.
- Refuse to silently choose among ambiguous candidates.
- Write machine-readable and human-readable audit outputs.

This script DOES NOT run Phase 13 biological models.
"""

from __future__ import annotations
import argparse
import csv
import json
import re
from pathlib import Path
from typing import Iterable

TEXT_EXT = {".csv", ".tsv", ".txt", ".json", ".xlsx", ".xls", ".parquet", ".gpkg", ".geojson", ".shp"}

PATTERNS = {
    "incidence_or_genus_cell": [
        r"genus.*cell", r"cell.*genus", r"incidence", r"25km", r"25_km"
    ],
    "trait_c3_n0": [
        r"\bC3\b", r"\bN0\b", r"trait", r"balloon"
    ],
    "environment_cell_table": [
        r"12C", r"environment", r"predictor", r"cell.*table", r"spatial"
    ],
    "boundary_or_ecoregion": [
        r"10[A-I]", r"break", r"boundary", r"ecoregion", r"biogeograph"
    ],
    "step11_outputs": [
        r"11G", r"11H", r"turnover", r"jaccard", r"simpson"
    ],
}

EXPECTED_REFERENCE = {
    "note": "Reference expectations from retained C3 workflow; 13A reports but does not force-match if upstream files document a newer retained state.",
    "total_genera_reference": 267,
    "occupied_25km_cells_reference": 205,
    "c3_genera_step11h_support": 87,
    "n0_genera_step11h_support": 140,
    "c3_occupied_cells_step11h": 148,
    "n0_occupied_cells_step11h": 176,
}

def score(path: Path, pats: list[str]) -> int:
    s = str(path).lower()
    return sum(bool(re.search(p.lower(), s, flags=re.I)) for p in pats)

def discover(root: Path) -> dict[str, list[dict]]:
    files = [p for p in root.rglob("*") if p.is_file() and p.suffix.lower() in TEXT_EXT]
    out = {}
    for category, pats in PATTERNS.items():
        ranked = []
        for p in files:
            sc = score(p, pats)
            if sc:
                ranked.append({
                    "path": str(p),
                    "score": sc,
                    "size_bytes": p.stat().st_size,
                    "suffix": p.suffix.lower(),
                })
        ranked.sort(key=lambda x: (-x["score"], x["size_bytes"], x["path"]))
        out[category] = ranked[:30]
    return out

def sniff_csv(path: Path) -> dict:
    result = {"path": str(path), "readable": False}
    try:
        if path.suffix.lower() not in {".csv", ".tsv"}:
            return result
        delim = "\t" if path.suffix.lower() == ".tsv" else ","
        with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as f:
            reader = csv.reader(f, delimiter=delim)
            header = next(reader, [])
            rows = []
            for i, row in enumerate(reader):
                rows.append(row)
                if i >= 4:
                    break
        result.update({
            "readable": True,
            "columns": header,
            "n_preview_rows": len(rows),
        })
    except Exception as e:
        result["error"] = repr(e)
    return result

def ambiguity_status(items: list[dict]) -> str:
    if not items:
        return "MISSING"
    if len(items) == 1:
        return "ONE_CANDIDATE"
    if items[0]["score"] > items[1]["score"]:
        return "TOP_CANDIDATE_BUT_REVIEW"
    return "AMBIGUOUS_REVIEW_REQUIRED"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    ap.add_argument("--output-dir", default=None)
    args = ap.parse_args()

    root = Path(args.project_root).expanduser().resolve()
    if not root.exists():
        raise SystemExit(f"FAIL: project root does not exist: {root}")

    outdir = Path(args.output_dir).expanduser().resolve() if args.output_dir else (
        root / "04_analysis_USE _THIS" / "13_historical_vs_contemporary" / "13A_input_audit"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    candidates = discover(root)
    report = {
        "phase": "13A",
        "project_root": str(root),
        "expected_reference": EXPECTED_REFERENCE,
        "categories": {},
        "overall_status": "REVIEW_REQUIRED",
        "important_rule": "No candidate is automatically accepted as the authoritative Phase 13 input solely by filename score.",
    }

    for cat, items in candidates.items():
        previews = []
        for item in items[:5]:
            previews.append(sniff_csv(Path(item["path"])))
        report["categories"][cat] = {
            "status": ambiguity_status(items),
            "candidates": items,
            "top_csv_tsv_previews": previews,
        }

    missing = [k for k, v in report["categories"].items() if v["status"] == "MISSING"]
    report["missing_categories"] = missing
    if missing:
        report["overall_status"] = "FAIL_MISSING_INPUT_CATEGORIES"

    json_path = outdir / "13A_input_audit.json"
    json_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    lines = [
        "PHASE 13A — INPUT / PROVENANCE AUDIT",
        "=" * 44,
        f"Project root: {root}",
        f"Status: {report['overall_status']}",
        "",
        "REFERENCE EXPECTATIONS (review only; not blindly enforced):",
    ]
    for k, v in EXPECTED_REFERENCE.items():
        lines.append(f"  {k}: {v}")
    lines += ["", "CANDIDATE INPUT CATEGORIES:"]
    for cat, info in report["categories"].items():
        lines.append(f"\n[{cat}] {info['status']}")
        for item in info["candidates"][:10]:
            lines.append(f"  score={item['score']} size={item['size_bytes']}  {item['path']}")
    lines += [
        "",
        "NEXT ACTION:",
        "Manually confirm authoritative files for incidence/traits, cell coordinates,",
        "Phase 12 environment, and a priori boundary configuration before Phase 13B.",
        "Do not run inferential Phase 13 models while any required category is missing",
        "or while authoritative inputs are ambiguous.",
    ]
    (outdir / "13A_INPUT_AUDIT_SUMMARY.txt").write_text("\n".join(lines), encoding="utf-8")

    print("\n".join(lines))
    print(f"\nJSON={json_path}")
    print(f"OUTPUT_DIR={outdir}")

if __name__ == "__main__":
    main()
