#!/usr/bin/env python3
"""
13H2_build_mulege_ecotone_predictors.py

Construct frozen local-scale predictors for the exploratory Mulegé /
central-Gulf ecotone analysis. This step performs NO C3/N0 hypothesis test.

Inputs:
- retained Step 10B 205-cell ecoregion crosswalk
- retained Step 10B all-cell ecoregion overlaps
- retained Step 10D rook-neighbor edges
- retained 13H1 cell table

Outputs:
- 205-cell predictor table
- frozen Mulegé focal definition
- predictor summary and join audit
- QC map
"""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

MULEGE_LAT = 26.89
MULEGE_LON = -111.98
PRIMARY_RADIUS_KM = 75.0
SENSITIVITY_RADII_KM = [50.0, 100.0, 125.0]

COSTA_GULFO = "Costa Central del Golfo"
SIERRA_GIGANTA = "Sierra de la Giganta"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0088
    p1 = math.radians(float(lat1))
    p2 = math.radians(float(lat2))
    dp = math.radians(float(lat2) - float(lat1))
    dl = math.radians(float(lon2) - float(lon1))
    a = (
        math.sin(dp / 2.0) ** 2
        + math.cos(p1) * math.cos(p2) * math.sin(dl / 2.0) ** 2
    )
    return 2.0 * r * math.asin(math.sqrt(a))


def safe_z(series):
    x = pd.to_numeric(series, errors="coerce")
    sd = x.std(ddof=1)
    if not np.isfinite(sd) or sd <= 0:
        return pd.Series(np.nan, index=x.index)
    return (x - x.mean()) / sd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--crosswalk", type=Path, default=None)
    p.add_argument("--overlaps", type=Path, default=None)
    p.add_argument("--neighbors", type=Path, default=None)
    p.add_argument("--cell-input", type=Path, default=None)
    p.add_argument("--output-dir", type=Path, default=None)
    return p.parse_args()


def main():
    args = parse_args()
    root = args.project_root.expanduser().resolve()
    step10 = (
        root / "04_analysis_USE _THIS" / "C3_pipeline_rebuild"
        / "09_C3_biogeographic_concordance"
    )

    cross_path = args.crosswalk or (
        step10 / "10B_cell_ecoregion_crosswalk"
        / "10B_cell_ecoregion_crosswalk.csv"
    )
    overlaps_path = args.overlaps or (
        step10 / "10B_cell_ecoregion_crosswalk"
        / "10B_all_cell_ecoregion_overlaps.csv"
    )
    neighbors_path = args.neighbors or (
        step10 / "10D_ecoregion_boundary_turnover"
        / "10D_neighbor_pair_metrics.csv"
    )
    cell_path = args.cell_input or (
        root / "04_analysis_USE _THIS" / "13_historical_vs_contemporary"
        / "13H_mulege_ecotone" / "13H1_input_audit"
        / "13H1_cell_input_audit.csv"
    )
    outdir = args.output_dir or (
        root / "04_analysis_USE _THIS" / "13_historical_vs_contemporary"
        / "13H_mulege_ecotone" / "13H2_ecotone_predictors"
    )

    cross_path = cross_path.expanduser().resolve()
    overlaps_path = overlaps_path.expanduser().resolve()
    neighbors_path = neighbors_path.expanduser().resolve()
    cell_path = cell_path.expanduser().resolve()
    outdir = outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print("PHASE 13H2 — BUILD MULEGÉ / ECOREGION-JUNCTION PREDICTORS")
    print("=" * 80)
    print(f"CROSSWALK  : {cross_path}")
    print(f"OVERLAPS   : {overlaps_path}")
    print(f"NEIGHBORS  : {neighbors_path}")
    print(f"CELL INPUT : {cell_path}")
    print(f"OUTPUT DIR : {outdir}")

    for p in [cross_path, overlaps_path, neighbors_path, cell_path]:
        if not p.exists():
            raise FileNotFoundError(p)

    cross = pd.read_csv(cross_path, low_memory=False)
    overlaps = pd.read_csv(overlaps_path, low_memory=False)
    edges = pd.read_csv(neighbors_path, low_memory=False)
    cells = pd.read_csv(cell_path, low_memory=False)

    cross = cross.rename(columns={"step10b_cell_id": "grid_cell_id"})
    overlaps = overlaps.rename(columns={"step10b_cell_id": "grid_cell_id"})
    for d in [cross, overlaps, cells]:
        d["grid_cell_id"] = d["grid_cell_id"].astype(str)
    edges["cell_i"] = edges["cell_i"].astype(str)
    edges["cell_j"] = edges["cell_j"].astype(str)

    n_analysis_cells = cells["grid_cell_id"].nunique()

    if len(cells) != n_analysis_cells:
        raise ValueError(
            "13H1 input contains duplicated grid_cell_id rows: "
            f"{len(cells)} rows but {n_analysis_cells} unique cells."
        )
    if n_analysis_cells < 150:
        raise ValueError(
            f"13H1 analysis set is unexpectedly small ({n_analysis_cells} cells)."
        )

    if len(cross) != 205 or cross["grid_cell_id"].nunique() != 205:
        raise ValueError(
            "Step 10B reference crosswalk must contain exactly 205 unique cells."
        )

    analysis_set = set(cells["grid_cell_id"])
    cross_set = set(cross["grid_cell_id"])
    missing_from_crosswalk = sorted(analysis_set - cross_set)
    if missing_from_crosswalk:
        raise ValueError(
            f"{len(missing_from_crosswalk)} 13H1 cells are absent from the "
            f"Step 10B crosswalk. First examples: {missing_from_crosswalk[:20]}"
        )

    # Step 10B contains all 205 occupied cells, whereas 13H1 inherits the
    # Phase 13C-eligible analysis set (expected 189). Restrict all predictor
    # construction and neighbor calculations to the exact 13H1 cell set.
    cross_analysis = cross.loc[
        cross["grid_cell_id"].isin(analysis_set)
    ].copy()

    print("\nCELL-SET ALIGNMENT")
    print(f"  Step 10B reference cells : {len(cross)}")
    print(f"  13H1 analysis cells      : {n_analysis_cells}")
    print(f"  Crosswalk cells retained : {len(cross_analysis)}")

    # Check centroids match.
    chk = cells[
        ["grid_cell_id", "centroid_latitude", "centroid_longitude"]
    ].merge(
        cross_analysis[["grid_cell_id", "centroid_latitude", "centroid_longitude"]],
        on="grid_cell_id", suffixes=("_h1", "_10b"), validate="one_to_one"
    )
    if (
        (chk["centroid_latitude_h1"] - chk["centroid_latitude_10b"]).abs().max()
        > 1e-6
        or
        (chk["centroid_longitude_h1"] - chk["centroid_longitude_10b"]).abs().max()
        > 1e-6
    ):
        raise ValueError("13H1 and Step 10B centroid coordinates differ.")

    out = cells.copy()

    cross_cols = [
        "grid_cell_id", "dominant_ecoregion", "second_ecoregion",
        "dominant_overlap_km2", "second_overlap_km2",
        "ecoregion_covered_km2", "full_cell_area_km2",
        "dominant_fraction_of_covered_land",
        "second_fraction_of_covered_land",
        "dominant_minus_second_fraction",
        "dominant_fraction_of_full_cell",
        "ecoregion_coverage_fraction",
        "centroid_ecoregion",
        "centroid_on_ecoregion_boundary",
        "centroid_outside_mapped_mainland",
        "centroid_agrees_with_dominant",
        "ambiguous_dominant_assignment",
        "low_mapped_land_coverage",
        "primary_assignment_eligible",
        "sensitivity_unambiguous_eligible",
    ]
    cross_cols = [c for c in cross_cols if c in cross_analysis.columns]
    out = out.merge(
        cross_analysis[cross_cols],
        on="grid_cell_id",
        validate="one_to_one"
    )

    # Within-cell overlap diversity.
    overlaps["overlap_km2"] = pd.to_numeric(overlaps["overlap_km2"], errors="coerce")
    overlaps = overlaps.loc[overlaps["overlap_km2"].fillna(0) > 0].copy()

    overlaps = overlaps.loc[
        overlaps["grid_cell_id"].isin(analysis_set)
    ].copy()

    covered = cross_analysis.set_index("grid_cell_id")["ecoregion_covered_km2"]
    full_area = cross_analysis.set_index("grid_cell_id")["full_cell_area_km2"]
    overlaps["covered_area"] = overlaps["grid_cell_id"].map(covered)
    overlaps["full_area"] = overlaps["grid_cell_id"].map(full_area)
    overlaps["frac_mapped"] = overlaps["overlap_km2"] / overlaps["covered_area"]
    overlaps["frac_full"] = overlaps["overlap_km2"] / overlaps["full_area"]

    mix_rows = []
    for cid, g in overlaps.groupby("grid_cell_id"):
        p = g["frac_mapped"].to_numpy(float)
        p = p[np.isfinite(p) & (p > 0)]
        if len(p):
            p = p / p.sum()
            shannon = float(-(p * np.log(p)).sum())
            effective = float(np.exp(shannon))
            simpson = float(1.0 - np.square(p).sum())
        else:
            shannon = effective = simpson = np.nan

        regions = set(g["ecoregion_label"].dropna().astype(str).str.strip())
        mix_rows.append({
            "grid_cell_id": cid,
            "ecoregion_overlap_count_any": int(len(g)),
            "ecoregion_overlap_count_ge1pct_mapped": int((g["frac_mapped"] >= 0.01).sum()),
            "ecoregion_overlap_count_ge5pct_mapped": int((g["frac_mapped"] >= 0.05).sum()),
            "ecoregion_overlap_count_ge1pct_full_cell": int((g["frac_full"] >= 0.01).sum()),
            "ecoregion_overlap_shannon": shannon,
            "ecoregion_effective_number": effective,
            "ecoregion_overlap_simpson_diversity": simpson,
            "overlaps_Costa_Central_del_Golfo": COSTA_GULFO in regions,
            "overlaps_Sierra_de_la_Giganta": SIERRA_GIGANTA in regions,
            "overlaps_Costa_Golfo_and_Sierra_Giganta": (
                COSTA_GULFO in regions and SIERRA_GIGANTA in regions
            ),
        })

    out = out.merge(pd.DataFrame(mix_rows), on="grid_cell_id", how="left")
    for c in [
        "ecoregion_overlap_count_any",
        "ecoregion_overlap_count_ge1pct_mapped",
        "ecoregion_overlap_count_ge5pct_mapped",
        "ecoregion_overlap_count_ge1pct_full_cell",
    ]:
        out[c] = out[c].fillna(0).astype(int)
    for c in [
        "overlaps_Costa_Central_del_Golfo",
        "overlaps_Sierra_de_la_Giganta",
        "overlaps_Costa_Golfo_and_Sierra_Giganta",
    ]:
        out[c] = out[c].fillna(False).astype(bool)

    out["dominant_second_Costa_Golfo_Sierra_Giganta_interface"] = out[
        ["dominant_ecoregion", "second_ecoregion"]
    ].apply(
        lambda r: (
            COSTA_GULFO in {str(x).strip() for x in r if pd.notna(x)}
            and SIERRA_GIGANTA in {str(x).strip() for x in r if pd.notna(x)}
        ),
        axis=1,
    )

    # Rook-neighbor junction structure.
    dominant = out.set_index("grid_cell_id")["dominant_ecoregion"].to_dict()
    cell_set = set(out["grid_cell_id"])
    neighbor_sets = defaultdict(set)
    edge_n = defaultdict(int)
    cross_n = defaultdict(int)
    interface_sets = defaultdict(set)

    for _, r in edges.iterrows():
        i, j = r["cell_i"], r["cell_j"]
        if i not in cell_set or j not in cell_set:
            continue
        neighbor_sets[i].add(j)
        neighbor_sets[j].add(i)
        edge_n[i] += 1
        edge_n[j] += 1
        if bool(r["crosses_ecoregion_boundary"]):
            cross_n[i] += 1
            cross_n[j] += 1
            interface_sets[i].add(str(r["ecoregion_pair"]))
            interface_sets[j].add(str(r["ecoregion_pair"]))

    local_rows = []
    for cid in out["grid_cell_id"]:
        local_regions = set()
        own = dominant.get(cid)
        if pd.notna(own):
            local_regions.add(str(own))
        for nid in neighbor_sets.get(cid, set()):
            reg = dominant.get(nid)
            if pd.notna(reg):
                local_regions.add(str(reg))
        deg = edge_n.get(cid, 0)
        cn = cross_n.get(cid, 0)
        local_rows.append({
            "grid_cell_id": cid,
            "rook_neighbor_degree": int(deg),
            "local_unique_dominant_ecoregions_focal_plus_neighbors": int(len(local_regions)),
            "adjacent_cross_boundary_edges": int(cn),
            "adjacent_cross_boundary_fraction": (cn / deg if deg else np.nan),
            "local_distinct_ecoregion_interfaces": int(len(interface_sets.get(cid, set()))),
        })
    out = out.merge(pd.DataFrame(local_rows), on="grid_cell_id", validate="one_to_one")

    # Frozen Mulegé focal geography.
    out["distance_to_Mulege_km"] = [
        haversine_km(lat, lon, MULEGE_LAT, MULEGE_LON)
        for lat, lon in zip(out["centroid_latitude"], out["centroid_longitude"])
    ]
    out["Mulege_focal_primary_75km"] = out["distance_to_Mulege_km"] <= 75
    for radius in SENSITIVITY_RADII_KM:
        out[f"Mulege_focal_sensitivity_{int(radius)}km"] = (
            out["distance_to_Mulege_km"] <= radius
        )
    out["central_Gulf_latitude_corridor_26p4_27p4"] = out[
        "centroid_latitude"
    ].between(26.4, 27.4, inclusive="both")

    # Frozen composite, created without C3_fraction.
    out["inverse_dominant_fraction"] = 1.0 - pd.to_numeric(
        out["dominant_fraction_of_covered_land"], errors="coerce"
    )
    components = {
        "z_effective_ecoregions": "ecoregion_effective_number",
        "z_overlap_count_ge1pct": "ecoregion_overlap_count_ge1pct_mapped",
        "z_local_unique_ecoregions": "local_unique_dominant_ecoregions_focal_plus_neighbors",
        "z_adjacent_cross_boundary_fraction": "adjacent_cross_boundary_fraction",
        "z_inverse_dominant_fraction": "inverse_dominant_fraction",
    }
    for zname, raw in components.items():
        out[zname] = safe_z(out[raw])
    zcols = list(components)
    out["ecoregion_junction_score_components_available"] = out[zcols].notna().sum(axis=1)
    out["ecoregion_junction_score_primary"] = out[zcols].mean(axis=1, skipna=True)
    out["ecoregion_junction_score_primary_eligible"] = (
        out["ecoregion_junction_score_components_available"] >= 4
    )
    out.loc[
        ~out["ecoregion_junction_score_primary_eligible"],
        "ecoregion_junction_score_primary"
    ] = np.nan

    if len(out) != n_analysis_cells or out["grid_cell_id"].nunique() != n_analysis_cells:
        raise ValueError(
            "Final table does not preserve the exact 13H1 analysis cell set: "
            f"{len(out)} rows, {out['grid_cell_id'].nunique()} unique cells, "
            f"expected {n_analysis_cells}."
        )

    analysis_edge_count = int(
        (
            edges["cell_i"].isin(analysis_set)
            & edges["cell_j"].isin(analysis_set)
        ).sum()
    )

    audit = pd.DataFrame([
        ["13H1 analysis cells", len(cells), n_analysis_cells, "PASS"],
        ["Step10B reference crosswalk cells", len(cross), 205, "PASS"],
        ["Step10B crosswalk cells retained", len(cross_analysis), n_analysis_cells, "PASS"],
        ["Step10B overlap rows retained", len(overlaps), "report", "PASS"],
        ["Step10D reference rook-neighbor edges", len(edges), 284,
         "PASS" if len(edges) == 284 else "REVIEW"],
        ["Step10D analysis-set rook-neighbor edges", analysis_edge_count, "report", "PASS"],
        ["final predictor cells", len(out), n_analysis_cells, "PASS"],
        ["junction-score eligible cells",
         int(out["ecoregion_junction_score_primary_eligible"].sum()),
         "report", "PASS"],
        ["Mulege 75-km focal cells",
         int(out["Mulege_focal_primary_75km"].sum()),
         "report", "PASS"],
    ], columns=["check", "observed", "expected", "status"])

    summary_rows = []
    for col in [
        "distance_to_Mulege_km",
        "ecoregion_overlap_count_any",
        "ecoregion_overlap_count_ge1pct_mapped",
        "ecoregion_overlap_count_ge5pct_mapped",
        "ecoregion_effective_number",
        "ecoregion_overlap_simpson_diversity",
        "local_unique_dominant_ecoregions_focal_plus_neighbors",
        "adjacent_cross_boundary_fraction",
        "local_distinct_ecoregion_interfaces",
        "inverse_dominant_fraction",
        "ecoregion_junction_score_primary",
    ]:
        s = pd.to_numeric(out[col], errors="coerce").dropna()
        summary_rows.append({
            "variable": col, "n": len(s), "mean": s.mean(),
            "sd": s.std(ddof=1), "median": s.median(),
            "min": s.min(), "max": s.max()
        })
    summary = pd.DataFrame(summary_rows)

    focal = pd.DataFrame([{
        "focal_name": "Heroica_Mulege",
        "latitude": MULEGE_LAT,
        "longitude": MULEGE_LON,
        "primary_radius_km": PRIMARY_RADIUS_KM,
        "sensitivity_radii_km": "50;100;125",
        "analysis_status": "EXPLORATORY_POST_HOC_LOCATION_FIXED_BEFORE_TESTING",
        "notes": (
            "The map feature was noticed visually; location and radii are frozen "
            "before 13H3 outcome testing."
        ),
    }])

    predictor_path = outdir / "13H2_ecotone_predictors_analysis_cells.csv"
    focal_path = outdir / "13H2_mulege_focal_definition.csv"
    summary_path = outdir / "13H2_predictor_summary.csv"
    audit_path = outdir / "13H2_input_join_audit.csv"
    readme_path = outdir / "13H2_README.txt"
    manifest_path = outdir / "13H2_run_manifest.txt"

    out.to_csv(predictor_path, index=False)
    focal.to_csv(focal_path, index=False)
    summary.to_csv(summary_path, index=False)
    audit.to_csv(audit_path, index=False)

    readme_path.write_text(
        """PHASE 13H2 — FROZEN MULEGÉ / ECOREGION-JUNCTION PREDICTORS

No C3/N0 hypothesis test was performed.

The Mulegé feature is exploratory/post hoc because it was noticed visually.
The focal center and 75-km primary radius, plus 50/100/125-km sensitivity
radii, are frozen here before testing.

The retained Step 10 ecoregion reference contains 205 occupied cells. Predictor
construction is restricted to the exact 13H1 / Phase 13C-eligible analysis set,
rather than requiring all 205 reference cells.

The primary junction score is an equal-weight mean of standardized:
1. effective number of overlapping ecoregions;
2. overlap richness at >=1% of mapped ecoregion area;
3. unique dominant ecoregions in focal cell plus rook neighbors;
4. fraction of neighboring edges crossing an ecoregion boundary;
5. inverse dominant ecoregion fraction.

Next: 13H3 tests the fixed local anomaly and junction association with spatial
and sampling controls.
""",
        encoding="utf-8",
    )

    manifest_path.write_text(
        "\n".join([
            "PHASE 13H2 RUN MANIFEST",
            f"python={sys.version.replace(chr(10), ' ')}",
            f"project_root={root}",
            f"crosswalk={cross_path}",
            f"crosswalk_sha256={sha256(cross_path)}",
            f"overlaps={overlaps_path}",
            f"overlaps_sha256={sha256(overlaps_path)}",
            f"neighbors={neighbors_path}",
            f"neighbors_sha256={sha256(neighbors_path)}",
            f"cell_input={cell_path}",
            f"cell_input_sha256={sha256(cell_path)}",
            f"n_reference_cells=205",
            f"n_analysis_cells={len(out)}",
            f"n_reference_neighbor_edges={len(edges)}",
            f"n_analysis_neighbor_edges={analysis_edge_count}",
            f"mulege_center={MULEGE_LAT},{MULEGE_LON}",
            "primary_radius_km=75",
            "sensitivity_radii_km=50,100,125",
            "outcome_used_to_construct_predictors=NO",
            "hypothesis_testing_performed=NO",
        ]) + "\n",
        encoding="utf-8",
    )

    png_path = outdir / "13H2_junction_score_map.png"
    svg_path = outdir / "13H2_junction_score_map.svg"
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 10))
        sc = ax.scatter(
            out["centroid_longitude"], out["centroid_latitude"],
            c=out["ecoregion_junction_score_primary"], s=35
        )
        ax.scatter([MULEGE_LON], [MULEGE_LAT], marker="*", s=180,
                   label="Heroica Mulegé")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title("Frozen ecoregion-junction score")
        ax.legend()
        fig.colorbar(sc, ax=ax, label="Ecoregion-junction score")
        fig.tight_layout()
        fig.savefig(png_path, dpi=300, bbox_inches="tight")
        fig.savefig(svg_path, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        print(f"WARNING: QC map not created: {e}")

    print("\nJOIN AUDIT")
    print(audit.to_string(index=False))

    print("\nFOCAL / INTERFACE COUNTS")
    for c in [
        "Mulege_focal_sensitivity_50km",
        "Mulege_focal_primary_75km",
        "Mulege_focal_sensitivity_100km",
        "Mulege_focal_sensitivity_125km",
        "overlaps_Costa_Central_del_Golfo",
        "overlaps_Sierra_de_la_Giganta",
        "overlaps_Costa_Golfo_and_Sierra_Giganta",
        "dominant_second_Costa_Golfo_Sierra_Giganta_interface",
    ]:
        print(f"  {c:58s}: {int(out[c].fillna(False).astype(bool).sum())}")

    print("\nFILES WRITTEN")
    for p in [
        predictor_path, focal_path, summary_path, audit_path,
        png_path, svg_path, readme_path, manifest_path
    ]:
        if p.exists():
            print(f"  {p}")

    print("\nSTATUS: PASS")
    print("13H2 froze Mulegé and ecoregion-junction predictors without testing C3/N0.")
    print("Next step: 13H3 local anomaly and ecoregion-junction inference.")


if __name__ == "__main__":
    main()
