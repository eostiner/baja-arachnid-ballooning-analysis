#!/usr/bin/env python3
"""
13H1_audit_mulege_ecotone_inputs.py

Purpose
-------
Audit retained inputs for a pre-specified local-scale Phase 13H test of the
apparent central Gulf / Mulegé bump and ecoregion-junction hypothesis.

This script performs NO hypothesis testing and does NOT select a Mulegé zone
from the observed C3/N0 outcome.

It audits:
1. Retained 25-km cell-level C3/N0 richness/proportion inputs.
2. Step 10 ecoregion assignment / polygon / boundary candidates.
3. Existing cell-level ecoregion metrics, if any.
4. Candidate geographic fields needed to distinguish Gulf vs Pacific sides.

Outputs
-------
13H1_candidate_files.csv
13H1_cell_input_audit.csv
13H1_audit_summary.txt
"""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd
import numpy as np


SEARCH_HINTS = (
    "ecoreg", "ecoregion", "gonzalez", "abraham", "10i", "10h", "10g",
    "biogeog", "boundary", "junction", "transition", "mulege", "mulegé",
    "ratio", "cell_ecoreg", "ecotone"
)

TABULAR_EXT = {".csv", ".tsv", ".txt", ".xlsx", ".xls"}
SPATIAL_EXT = {".shp", ".geojson", ".gpkg", ".json"}


def parse_args():
    p = argparse.ArgumentParser(
        description="Audit retained inputs for Phase 13H Mulegé/ecotone test."
    )
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--output-dir", type=Path, default=None)
    return p.parse_args()


def safe_read_table(path: Path):
    try:
        if path.suffix.lower() == ".csv":
            return pd.read_csv(path, low_memory=False)
        if path.suffix.lower() == ".tsv":
            return pd.read_csv(path, sep="\t", low_memory=False)
        if path.suffix.lower() in {".xlsx", ".xls"}:
            return pd.read_excel(path)
        if path.suffix.lower() == ".txt":
            # Try comma, tab, then whitespace-ish fallback.
            for sep in [",", "\t"]:
                try:
                    df = pd.read_csv(path, sep=sep, low_memory=False)
                    if df.shape[1] > 1:
                        return df
                except Exception:
                    pass
    except Exception:
        return None
    return None


def main():
    args = parse_args()
    root = args.project_root.expanduser().resolve()
    outdir = args.output_dir or (
        root / "04_analysis_USE _THIS" / "13_historical_vs_contemporary"
        / "13H_mulege_ecotone" / "13H1_input_audit"
    )
    outdir = outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print("PHASE 13H1 — MULEGÉ / ECOREGION-JUNCTION INPUT AUDIT")
    print("=" * 80)
    print(f"PROJECT ROOT : {root}")
    print(f"OUTPUT DIR   : {outdir}")

    # Core retained cell data.
    richness_path = (
        root / "04_analysis_USE _THIS" / "13_historical_vs_contemporary"
        / "13D_paired_C3_N0_community_dissimilarity"
        / "13D_cell_trait_richness.csv"
    )
    env_path = (
        root / "04_analysis_USE _THIS" / "12C_cell_environment_model_table"
        / "12C_cell_environment_model_table.csv"
    )
    lookup_path = (
        root / "02_data_clean" / "08_grid25km_incidence"
        / "10_common_grid25km_cell_lookup.csv"
    )

    for p in [richness_path, env_path, lookup_path]:
        if not p.exists():
            raise FileNotFoundError(f"Required retained input not found: {p}")

    richness = pd.read_csv(richness_path)
    env = pd.read_csv(env_path)
    lookup = pd.read_csv(lookup_path)

    if "grid_cell_id" not in richness.columns:
        raise ValueError("13D richness table missing grid_cell_id.")

    richness["classified_richness_C3_plus_N0"] = (
        richness["C3_richness"] + richness["N0_richness"]
    )
    richness["C3_fraction"] = np.where(
        richness["classified_richness_C3_plus_N0"] > 0,
        richness["C3_richness"] / richness["classified_richness_C3_plus_N0"],
        np.nan
    )

    print("\nCORE CELL INPUT")
    print(f"  13D cells: {len(richness)}")
    print(f"  cells with classified richness > 0: "
          f"{int((richness['classified_richness_C3_plus_N0'] > 0).sum())}")
    print(f"  latitude range: "
          f"{richness['centroid_latitude'].min():.3f} to "
          f"{richness['centroid_latitude'].max():.3f}")
    print(f"  longitude range: "
          f"{richness['centroid_longitude'].min():.3f} to "
          f"{richness['centroid_longitude'].max():.3f}")

    cell_audit = richness[[
        "grid_cell_id", "grid_cell_order", "centroid_latitude",
        "centroid_longitude", "latitude_band",
        "C3_richness", "N0_richness",
        "classified_richness_C3_plus_N0", "C3_fraction"
    ]].copy()
    cell_audit.to_csv(outdir / "13H1_cell_input_audit.csv", index=False)

    # Search likely retained Step 10 / ecoregion files.
    candidates = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue

        suffix = path.suffix.lower()
        if suffix not in TABULAR_EXT | SPATIAL_EXT:
            continue

        rel = str(path.relative_to(root))
        low = rel.lower()

        # Prioritize retained analysis and obvious ecoregion hints.
        hint_score = sum(h in low for h in SEARCH_HINTS)
        retained_bonus = 2 if "04_analysis_use _this" in low else 0
        step10_bonus = 2 if "/10" in low or "\\10" in low else 0
        score = hint_score + retained_bonus + step10_bonus

        if score <= 0:
            continue

        rec = {
            "path": str(path),
            "relative_path": rel,
            "suffix": suffix,
            "score": score,
            "size_bytes": path.stat().st_size,
            "readable_tabular": False,
            "n_rows": np.nan,
            "n_cols": np.nan,
            "columns": "",
            "has_grid_cell_id": False,
            "has_ecoregion_like_column": False,
            "has_latitude": False,
            "has_longitude": False,
            "notes": "",
        }

        if suffix in TABULAR_EXT:
            df = safe_read_table(path)
            if df is not None:
                cols = [str(c) for c in df.columns]
                lcols = [c.lower() for c in cols]
                rec["readable_tabular"] = True
                rec["n_rows"] = len(df)
                rec["n_cols"] = len(cols)
                rec["columns"] = " | ".join(cols[:80])
                rec["has_grid_cell_id"] = any(c == "grid_cell_id" for c in lcols)
                rec["has_ecoregion_like_column"] = any(
                    ("ecoreg" in c) or ("region" in c) or ("ecotone" in c)
                    for c in lcols
                )
                rec["has_latitude"] = any("lat" in c for c in lcols)
                rec["has_longitude"] = any(("lon" in c) or ("long" in c) for c in lcols)

                # Flag especially useful tables.
                if rec["has_grid_cell_id"] and rec["has_ecoregion_like_column"]:
                    rec["notes"] = "HIGH_VALUE_CELL_ECOREGION_CANDIDATE"
                elif rec["has_ecoregion_like_column"]:
                    rec["notes"] = "ECOREGION_TABLE_CANDIDATE"
                elif rec["has_grid_cell_id"]:
                    rec["notes"] = "CELL_LEVEL_CANDIDATE"

        elif suffix in SPATIAL_EXT:
            rec["notes"] = "SPATIAL_GEOMETRY_CANDIDATE"

        candidates.append(rec)

    cand = pd.DataFrame(candidates)
    if len(cand):
        cand = cand.sort_values(
            ["score", "notes", "size_bytes"],
            ascending=[False, True, False]
        ).reset_index(drop=True)
    cand_path = outdir / "13H1_candidate_files.csv"
    cand.to_csv(cand_path, index=False)

    print("\nTOP CANDIDATE FILES")
    if len(cand):
        show_cols = [
            "score", "notes", "n_rows", "n_cols", "relative_path"
        ]
        print(cand[show_cols].head(40).to_string(index=False))
    else:
        print("  No ecoregion/Step10 candidates found.")

    # Explicitly inspect top high-value tabular candidates.
    high = cand.loc[
        cand["notes"].isin([
            "HIGH_VALUE_CELL_ECOREGION_CANDIDATE",
            "ECOREGION_TABLE_CANDIDATE"
        ])
    ].head(15)

    print("\nHIGH-VALUE TABLE DETAILS")
    if len(high) == 0:
        print("  None automatically identified.")
    else:
        for _, r in high.iterrows():
            print("\n" + "-" * 80)
            print(r["path"])
            df = safe_read_table(Path(r["path"]))
            if df is None:
                print("Could not read.")
                continue
            print(f"SHAPE: {df.shape}")
            print("COLUMNS:")
            print(list(df.columns))
            print("FIRST 5 ROWS:")
            print(df.head(5).to_string(index=False))

    # Audit whether likely side/coast fields already exist.
    combined_cols = sorted(set(map(str, env.columns)) | set(map(str, lookup.columns)))
    coast_like = [
        c for c in combined_cols
        if any(k in c.lower() for k in [
            "coast", "gulf", "pacific", "side", "ecoreg", "region", "province"
        ])
    ]

    print("\nEXISTING COAST / SIDE / ECOREGION-LIKE FIELDS IN CORE TABLES")
    if coast_like:
        for c in coast_like:
            print(f"  {c}")
    else:
        print("  None detected automatically.")

    # Write summary.
    summary_lines = [
        "PHASE 13H1 INPUT AUDIT SUMMARY",
        f"project_root={root}",
        f"n_13D_cells={len(richness)}",
        f"n_candidate_files={len(cand)}",
        f"n_high_value_cell_ecoregion_candidates="
        f"{int((cand['notes'] == 'HIGH_VALUE_CELL_ECOREGION_CANDIDATE').sum()) if len(cand) else 0}",
        f"n_ecoregion_table_candidates="
        f"{int((cand['notes'] == 'ECOREGION_TABLE_CANDIDATE').sum()) if len(cand) else 0}",
        "",
        "NEXT DECISION:",
        "Use independently defined retained ecoregion geometry/assignments to build:",
        "1) local ecoregion-junction complexity around each 25-km cell,",
        "2) an independently fixed Mulege/central-Gulf focal zone,",
        "3) a Gulf-side nonlinear north-south test.",
        "",
        "NO OUTCOME-BASED ZONE SELECTION WAS PERFORMED IN THIS AUDIT."
    ]
    (outdir / "13H1_audit_summary.txt").write_text(
        "\n".join(summary_lines) + "\n", encoding="utf-8"
    )

    print("\nFILES WRITTEN")
    print(f"  {cand_path}")
    print(f"  {outdir / '13H1_cell_input_audit.csv'}")
    print(f"  {outdir / '13H1_audit_summary.txt'}")

    print("\nSTATUS: PASS")
    print("13H1 audited Mulegé/ecotone inputs without testing the outcome.")
    print("Next: review the top retained ecoregion candidates, then build 13H2.")


if __name__ == "__main__":
    main()
