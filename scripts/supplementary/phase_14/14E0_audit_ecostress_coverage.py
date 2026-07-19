#!/usr/bin/env python3
"""Phase 14E0 — audit whether Earth Engine ECOSTRESS currently covers eligible Baja cells."""
from __future__ import annotations

import argparse
import os
from datetime import datetime, timezone
from pathlib import Path

from phase14_common import default_analysis_output_root, write_csv, write_json
from importlib.util import spec_from_file_location, module_from_spec

SCRIPT_VERSION = "14E0_v0.3.0_2026-07-18"
ASSET = "NASA/ECOSTRESS/L2T_LSTE/V2"


def load_c1_helpers():
    path = Path(__file__).resolve().parent / "14C1_extract_real_stress_earth_engine.py"
    spec = spec_from_file_location("phase14_c1_helpers", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import Earth Engine helpers from {path}")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.initialize_earth_engine, module.load_cell_features


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit ECOSTRESS coverage for eligible Phase 14 cells.")
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--ee-project", default=os.environ.get("EE_PROJECT") or os.environ.get("EARTHENGINE_PROJECT") or os.environ.get("GOOGLE_CLOUD_PROJECT"))
    parser.add_argument("--authenticate", action="store_true")
    parser.add_argument("--cells-file", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--start-date", default="2018-07-09")
    parser.add_argument("--end-date", default="2026-01-01")
    args = parser.parse_args()

    root = args.project_root.expanduser().resolve()
    base = default_analysis_output_root(root)
    c1_dir = base / "14C_real_stress_anomalies"
    cells_path = args.cells_file.expanduser().resolve() if args.cells_file else c1_dir / "14C0_eligible_temporal_cells.geojson"
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else base / "14E0_ecostress_coverage_audit"
    outdir.mkdir(parents=True, exist_ok=True)

    initialize_ee, load_cells = load_c1_helpers()
    ee = initialize_ee(args.ee_project, args.authenticate)
    cells_fc, cell_rows = load_cells(ee, cells_path)
    cell_ids = sorted({str(row["grid_cell_id"]) for row in cell_rows})

    collection = (
        ee.ImageCollection(ASSET)
        .filterDate(args.start_date, args.end_date)
        .filterBounds(cells_fc.geometry())
    )
    intersecting_scenes = int(collection.size().getInfo())
    rows = []
    if intersecting_scenes == 0:
        rows = [{"grid_cell_id": cell, "intersecting_scene_count": 0, "usable_for_phase14": 0} for cell in cell_ids]
        status_name = "NO_EARTH_ENGINE_ECOSTRESS_COVERAGE_FOR_ELIGIBLE_BAJA_CELLS"
    else:
        # Scene count per cell is a coverage audit only; it is not a valid temperature summary.
        for row in cell_rows:
            cell = str(row["grid_cell_id"])
            feature = cells_fc.filter(ee.Filter.eq("grid_cell_id", cell)).first()
            count = int(collection.filterBounds(feature.geometry()).size().getInfo())
            rows.append({
                "grid_cell_id": cell,
                "intersecting_scene_count": count,
                "usable_for_phase14": int(count >= 10),
            })
        status_name = "COVERAGE_PRESENT_REQUIRES_SEPARATE_QUALITY_AND_TEMPORAL_AUDIT"

    write_csv(outdir / "14E0_ecostress_scene_coverage_by_cell.csv", rows)
    status = {
        "phase": "14E0",
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": status_name,
        "asset": ASSET,
        "date_range": [args.start_date, args.end_date],
        "eligible_cells": len(cell_ids),
        "intersecting_scenes_total": intersecting_scenes,
        "decision": (
            "Do not add ECOSTRESS to H3 unless coverage is spatially and temporally adequate for repeated cell-period comparisons. "
            "A scene intersection alone is not sufficient."
        ),
    }
    write_json(outdir / "14E0_run_status.json", status)
    (outdir / "14E0_README.txt").write_text(
        "PHASE 14E0 — ECOSTRESS COVERAGE AUDIT\n"
        "=======================================\n"
        f"Status: {status_name}\n"
        f"Eligible cells: {len(cell_ids)}\n"
        f"Intersecting scenes: {intersecting_scenes}\n\n"
        "ECOSTRESS is not admitted into the H3 model unless this audit and a subsequent quality audit pass.\n",
        encoding="utf-8",
    )
    print("PHASE 14E0 — ECOSTRESS COVERAGE AUDIT")
    print(f"Status: {status_name}")
    print(f"Intersecting scenes: {intersecting_scenes}")
    print(f"OUTPUT_DIR={outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
