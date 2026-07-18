#!/usr/bin/env python3
"""
Step 11H — Separate C3 and N0 geographic grouping analysis.

Purpose
-------
Analyze ballooning-capable (C3 = D1+D2+D3) and fixed non-ballooning (N0)
assemblages separately to ask whether latitude bands and ecoregions group
differently by genus composition.

Primary metric: Jaccard dissimilarity.
Sensitivity metric: Simpson replacement.
Standardization: equal-cell Monte Carlo resampling using the retained 25-km grid.
Clustering: average-linkage hierarchical clustering of median resampled Jaccard matrices.
Ordination: classical PCoA of median Jaccard matrices.
C3-vs-N0 structural similarity: Spearman Mantel permutation test.

This is exploratory/supplementary and does not replace paired Step 11G.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon as MplPolygon
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr


BAND_ORDER = ["23-24N", "24-26N", "26-28N", "28-30N", "30-32N"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--project-root", required=True, help="Full Baja Ballooning project root.")
    p.add_argument(
        "--output-dir",
        default=None,
        help="Defaults to <project-root>/04_analysis_USE _THIS/11H_separate_trait_grouping",
    )
    p.add_argument("--iterations", type=int, default=5000)
    p.add_argument("--mantel-permutations", type=int, default=9999)
    p.add_argument("--seed", type=int, default=20260717)
    return p.parse_args()


def clean(value) -> str:
    return str(value or "").strip()


def norm(value) -> str:
    return re.sub(r"[^a-z0-9]+", "", clean(value).casefold())


def find_file(root: Path, exact_relatives: list[str], filename: str) -> Path:
    for rel in exact_relatives:
        p = root / rel
        if p.is_file():
            return p

    candidates = []
    for p in root.glob(f"**/{filename}"):
        s = str(p).casefold()
        if any(
            token in s
            for token in [
                "/08_archive/",
                "/archive/",
                "__macosx",
                "_incorrect_",
                "/alternative",
                "/backup",
                "/old_",
            ]
        ):
            continue
        candidates.append(p)

    if not candidates:
        raise FileNotFoundError(f"Could not locate {filename} under {root}")

    # Prefer retained/current analysis paths.
    def score(p: Path):
        s = str(p).casefold()
        return (
            "04_analysis_use _this" in s,
            "c3_pipeline_rebuild" in s,
            p.stat().st_mtime,
        )

    return sorted(candidates, key=score, reverse=True)[0]


def jaccard(a: set[str], b: set[str]) -> float:
    u = len(a | b)
    return 0.0 if u == 0 else 1.0 - len(a & b) / u


def simpson(a: set[str], b: set[str]) -> float:
    shared = len(a & b)
    only_a = len(a - b)
    only_b = len(b - a)
    den = shared + min(only_a, only_b)
    return 0.0 if den == 0 else min(only_a, only_b) / den


def aggregate_set(cell_to_genera, cells):
    out = set()
    for cell in cells:
        out |= cell_to_genera.get(cell, set())
    return out


def resample_matrices(
    group_cells: dict[str, list[str]],
    labels: list[str],
    cell_to_genera: dict[str, set[str]],
    equal_n: int,
    metric,
    iterations: int,
    seed: int,
):
    rng = np.random.default_rng(seed)
    n = len(labels)
    mats = np.empty((iterations, n, n), dtype=float)

    for it in range(iterations):
        sets = []
        for label in labels:
            cells = group_cells[label]
            if len(cells) < equal_n:
                raise ValueError(
                    f"{label} has only {len(cells)} eligible cells; needs {equal_n}."
                )
            sampled = (
                cells
                if len(cells) == equal_n
                else rng.choice(cells, size=equal_n, replace=False).tolist()
            )
            sets.append(aggregate_set(cell_to_genera, sampled))

        mat = np.zeros((n, n), dtype=float)
        for i in range(n):
            for j in range(i + 1, n):
                d = metric(sets[i], sets[j])
                mat[i, j] = d
                mat[j, i] = d
        mats[it] = mat

    return mats


def save_matrix_outputs(mats, labels, prefix: Path):
    med = np.median(mats, axis=0)
    lo = np.quantile(mats, 0.025, axis=0)
    hi = np.quantile(mats, 0.975, axis=0)

    pd.DataFrame(med, index=labels, columns=labels).to_csv(
        prefix.with_name(prefix.name + "_median_matrix.csv")
    )

    rows = []
    for i, j in itertools.combinations(range(len(labels)), 2):
        rows.append(
            {
                "group_1": labels[i],
                "group_2": labels[j],
                "median": med[i, j],
                "ci_low_2.5": lo[i, j],
                "ci_high_97.5": hi[i, j],
            }
        )
    pd.DataFrame(rows).to_csv(
        prefix.with_name(prefix.name + "_pairwise_intervals.csv"), index=False
    )
    return med


def pcoa(distance_matrix):
    n = distance_matrix.shape[0]
    d2 = distance_matrix ** 2
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ d2 @ j
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    positive = eigvals > 1e-12
    vals = eigvals[positive]
    vecs = eigvecs[:, positive]
    coords = vecs * np.sqrt(vals)
    explained = vals / vals.sum() if vals.sum() else np.zeros_like(vals)
    return coords, explained


def plot_dendrogram(matrix, labels, title, outfile, horizontal=False):
    z = linkage(squareform(matrix, checks=False), method="average")
    fig, ax = plt.subplots(figsize=(9, 6 if horizontal else 5.5))
    if horizontal:
        dendrogram(z, labels=labels, orientation="right", ax=ax, leaf_font_size=9)
        ax.set_xlabel("Average-linkage Jaccard dissimilarity")
    else:
        dendrogram(z, labels=labels, ax=ax, leaf_rotation=30, leaf_font_size=10)
        ax.set_ylabel("Average-linkage Jaccard dissimilarity")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return z


def plot_pcoa(matrix, labels, title, outfile):
    coords, explained = pcoa(matrix)
    x = coords[:, 0] if coords.shape[1] >= 1 else np.zeros(len(labels))
    y = coords[:, 1] if coords.shape[1] >= 2 else np.zeros(len(labels))
    p1 = explained[0] * 100 if len(explained) >= 1 else 0
    p2 = explained[1] * 100 if len(explained) >= 2 else 0

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(x, y, s=50)
    for xi, yi, label in zip(x, y, labels):
        ax.annotate(label, (xi, yi), xytext=(5, 4), textcoords="offset points", fontsize=9)
    ax.axhline(0, linewidth=0.7)
    ax.axvline(0, linewidth=0.7)
    ax.set_xlabel(f"PCoA 1 ({p1:.1f}%)")
    ax.set_ylabel(f"PCoA 2 ({p2:.1f}%)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return p1, p2


def grouping_stability(mats, labels):
    pairs = list(itertools.combinations(range(len(labels)), 2))
    nearest = Counter()
    distinct = Counter()

    for mat in mats:
        pair_values = np.array([mat[i, j] for i, j in pairs])
        minimum = pair_values.min()
        hits = np.where(np.isclose(pair_values, minimum))[0]
        for hit in hits:
            i, j = pairs[hit]
            nearest[(labels[i], labels[j])] += 1 / len(hits)

        mean_distance = np.array(
            [np.mean(np.delete(mat[i], i)) for i in range(len(labels))]
        )
        maximum = mean_distance.max()
        hits = np.where(np.isclose(mean_distance, maximum))[0]
        for hit in hits:
            distinct[labels[hit]] += 1 / len(hits)

    total = len(mats)
    nearest_rows = [
        (a, b, count / total) for (a, b), count in nearest.most_common()
    ]
    distinct_rows = [
        (label, count / total) for label, count in distinct.most_common()
    ]
    return nearest_rows, distinct_rows


def mantel_spearman(mat1, mat2, permutations, seed):
    rng = np.random.default_rng(seed)
    tri = np.triu_indices_from(mat1, k=1)
    obs = float(spearmanr(mat1[tri], mat2[tri]).statistic)
    null = np.empty(permutations)
    n = mat1.shape[0]

    for i in range(permutations):
        perm = rng.permutation(n)
        pmat = mat2[np.ix_(perm, perm)]
        null[i] = spearmanr(mat1[tri], pmat[tri]).statistic

    p = (1 + np.sum(np.abs(null) >= abs(obs))) / (permutations + 1)
    return obs, float(p)



def geojson_polygons(geometry):
    """Return exterior rings for Polygon or MultiPolygon GeoJSON geometry."""
    gtype = geometry.get("type")
    coords = geometry.get("coordinates", [])
    if gtype == "Polygon":
        return [coords[0]] if coords else []
    if gtype == "MultiPolygon":
        return [poly[0] for poly in coords if poly]
    return []



def load_mainland_outline(mainland_outline_path: Path):
    """
    Read the validated mainland Baja outline and return exterior rings in EPSG:4326.

    The retained pipeline's preferred input is:
      10A_mainland_outline_largest_component.gpkg

    GeoPandas is imported locally so Step 11H's non-map analyses remain isolated
    from GIS dependencies until the publication map is created.
    """
    try:
        import geopandas as gpd
    except ImportError as exc:
        raise RuntimeError(
            "GeoPandas is required to draw the Baja mainland outline. "
            "Install geopandas in the pipeline environment."
        ) from exc

    gdf = gpd.read_file(mainland_outline_path)
    if gdf.empty:
        raise RuntimeError(
            f"Mainland outline file contains no features: {mainland_outline_path}"
        )

    if gdf.crs is None:
        raise RuntimeError(
            f"Mainland outline has no CRS: {mainland_outline_path}"
        )

    gdf = gdf.to_crs("EPSG:4326")

    rings = []
    for geom in gdf.geometry:
        if geom is None or geom.is_empty:
            continue
        if geom.geom_type == "Polygon":
            rings.append(list(geom.exterior.coords))
        elif geom.geom_type == "MultiPolygon":
            for poly in geom.geoms:
                rings.append(list(poly.exterior.coords))

    if not rings:
        raise RuntimeError(
            f"No Polygon/MultiPolygon geometry found in: {mainland_outline_path}"
        )

    return rings


def draw_mainland_outline(ax, outline_rings):
    """Draw a clean Baja mainland silhouette behind the analytical grid."""
    patches = [
        MplPolygon(ring, closed=True)
        for ring in outline_rings
        if len(ring) >= 3
    ]
    if not patches:
        return

    # Light land fill establishes the peninsula silhouette.
    land = PatchCollection(
        patches,
        facecolor="white",
        edgecolor="none",
        linewidths=0,
        zorder=0,
    )
    ax.add_collection(land)

    # Strong exterior coastline/border outline.
    coast = PatchCollection(
        patches,
        facecolor="none",
        edgecolor="black",
        linewidths=1.25,
        zorder=8,
    )
    ax.add_collection(coast)


def build_baja_comparison_map(
    geojson_path: Path,
    mainland_outline_path: Path,
    lookup: pd.DataFrame,
    latitude_results: dict,
    out: Path,
    iterations: int,
    equal_cells: int,
):
    """
    Build the Step 11H two-panel Baja map directly from the current analysis, using the validated mainland Baja outline.

    Left: C3 ballooning-capable assemblages.
    Right: N0 fixed non-ballooning assemblages.

    The nearest-pair bands and most-distinct band are calculated from the
    current Monte Carlo results, so the figure never depends on hard-coded
    biological results.
    """
    geo = json.loads(geojson_path.read_text(encoding="utf-8"))
    outline_rings = load_mainland_outline(mainland_outline_path)

    cell_polys = {}
    all_xy = []

    for feature in geo.get("features", []):
        props = feature.get("properties", {})
        cell_id = None
        for key in ["grid_cell_id", "cell_id", "id"]:
            if key in props and props[key] is not None:
                cell_id = str(props[key])
                break
        if cell_id is None:
            continue

        polys = geojson_polygons(feature.get("geometry", {}))
        if not polys:
            continue

        cell_polys[cell_id] = polys
        for poly in polys:
            all_xy.extend(poly)

    if not all_xy:
        raise RuntimeError(f"No usable cell polygons found in {geojson_path}")

    xs = np.array([point[0] for point in all_xy], dtype=float)
    ys = np.array([point[1] for point in all_xy], dtype=float)
    xmin, xmax = float(xs.min()), float(xs.max())
    ymin, ymax = float(ys.min()), float(ys.max())

    display = {
        "23-24N": "23–24°N",
        "24-26N": "24–26°N",
        "26-28N": "26–28°N",
        "28-30N": "28–30°N",
        "30-32N": "30–32°N",
    }
    band_spans = {
        "23-24N": (23, 24),
        "24-26N": (24, 26),
        "26-28N": (26, 28),
        "28-30N": (28, 30),
        "30-32N": (30, 32),
    }
    band_centers = {
        band: sum(span) / 2 for band, span in band_spans.items()
    }

    def draw_panel(ax, cls, title):
        # Draw the actual Baja mainland silhouette first so cells and band
        # annotations are clearly situated within the peninsula.
        draw_mainland_outline(ax, outline_rings)

        result = latitude_results[cls]
        nearest = result["nearest"][0]
        distinct = result["distinct"][0]

        nearest_bands = {nearest[0], nearest[1]}
        distinct_band = distinct[0]

        for band in BAND_ORDER:
            ids = lookup.loc[
                lookup["centroid_latitude_band"].astype(str) == band,
                "grid_cell_id",
            ].dropna().astype(str)

            patches = []
            for cell_id in ids:
                for poly in cell_polys.get(cell_id, []):
                    patches.append(MplPolygon(poly, closed=True))

            if not patches:
                continue

            if band == distinct_band:
                collection = PatchCollection(
                    patches,
                    alpha=0.35,
                    hatch="////",
                    linewidths=0.45,
                    zorder=3,
                )
            elif band in nearest_bands:
                collection = PatchCollection(
                    patches,
                    alpha=0.18,
                    linewidths=0.75,
                    zorder=3,
                )
            else:
                collection = PatchCollection(
                    patches,
                    alpha=0.07,
                    linewidths=0.35,
                    zorder=2,
                )
            ax.add_collection(collection)

        for latitude in [24, 26, 28, 30]:
            ax.axhline(latitude, linewidth=0.7, linestyle="--", zorder=5)

        for band in BAND_ORDER:
            suffix = ""
            if band == distinct_band:
                suffix = "  ← most distinct"
            elif band in nearest_bands:
                suffix = "  ← closest-pair band"

            ax.text(
                xmax + 0.12,
                band_centers[band],
                f"{display[band]}{suffix}",
                va="center",
                fontsize=8.5,
            )

        y0 = min(band_spans[nearest[0]][0], band_spans[nearest[1]][0])
        y1 = max(band_spans[nearest[0]][1], band_spans[nearest[1]][1])
        xbr = xmin - 0.15

        ax.plot([xbr, xbr], [y0, y1], linewidth=2.2)
        ax.plot([xbr, xbr + 0.08], [y0, y0], linewidth=2.2)
        ax.plot([xbr, xbr + 0.08], [y1, y1], linewidth=2.2)
        ax.text(
            xbr - 0.04,
            (y0 + y1) / 2,
            f"Closest pair\n{nearest[2] * 100:.1f}%",
            ha="right",
            va="center",
            fontsize=8.5,
        )

        ax.text(
            xmin + 0.15,
            band_centers[distinct_band],
            f"Most distinct\n{distinct[1] * 100:.1f}% of resamples",
            fontsize=9,
            va="center",
        )

        ax.set_xlim(xmin - 0.8, xmax + 2.1)
        ax.set_ylim(max(22.8, ymin - 0.1), min(32.3, ymax + 0.1))
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(title, fontsize=13)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

    fig, axes = plt.subplots(1, 2, figsize=(12, 9), sharey=True)

    draw_panel(axes[0], "C3", "Ballooning (C3)")
    draw_panel(axes[1], "N0", "Non-ballooning (N0)")

    fig.suptitle(
        "Separate geographic organization of ballooning and non-ballooning arachnid assemblages",
        fontsize=15,
        y=0.96,
    )

    c3 = latitude_results["C3"]
    n0 = latitude_results["N0"]

    fig.text(
        0.5,
        0.075,
        (
            f"C3 assemblages were most often isolated at {display[c3['distinct'][0][0]]}, "
            f"whereas N0 assemblages were most often isolated at "
            f"{display[n0['distinct'][0][0]]}. "
            "The panels summarize each trait group independently."
        ),
        ha="center",
        va="center",
        fontsize=10.5,
        wrap=True,
    )

    fig.text(
        0.5,
        0.035,
        (
            f"Closest-pair and most-distinct frequencies are from {iterations:,} "
            f"equal-cell Monte Carlo resamples ({equal_cells} occupied 25-km cells per latitude band)."
        ),
        ha="center",
        va="center",
        fontsize=9,
    )

    fig.tight_layout(rect=[0.03, 0.12, 0.98, 0.93])

    png = out / "11H_Baja_C3_vs_N0_geographic_grouping.png"
    svg = out / "11H_Baja_C3_vs_N0_geographic_grouping_EDITABLE.svg"
    pdf = out / "11H_Baja_C3_vs_N0_geographic_grouping.pdf"

    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    return png, svg, pdf

def main():
    args = parse_args()
    root = Path(args.project_root).expanduser().resolve()

    if args.output_dir:
        out = Path(args.output_dir).expanduser().resolve()
    else:
        preferred = root / "04_analysis_USE _THIS"
        base = preferred if preferred.is_dir() else root / "04_analysis"
        out = base / "11H_separate_trait_grouping"

    out.mkdir(parents=True, exist_ok=True)

    incidence_path = find_file(
        root,
        ["02_data_clean/08_grid25km_incidence/10_ballooning_final_genus_grid25km_incidence_long.csv"],
        "10_ballooning_final_genus_grid25km_incidence_long.csv",
    )
    lookup_path = find_file(
        root,
        ["02_data_clean/08_grid25km_incidence/10_common_grid25km_cell_lookup.csv"],
        "10_common_grid25km_cell_lookup.csv",
    )
    trait_path = find_file(
        root,
        [
            "04_analysis_USE _THIS/11G_trait_partitioned_equal_cell/11G_normalized_trait_lookup.csv",
            "04_analysis/11G_trait_partitioned_equal_cell/11G_normalized_trait_lookup.csv",
        ],
        "11G_normalized_trait_lookup.csv",
    )
    crosswalk_path = find_file(
        root,
        [],
        "10B_cell_ecoregion_crosswalk.csv",
    )
    summary10c_path = find_file(
        root,
        [],
        "10C_equal_cell_summary.csv",
    )

    geojson_path = find_file(
        root,
        ["02_data_clean/08_grid25km_incidence/10_common_grid25km_cells.geojson"],
        "10_common_grid25km_cells.geojson",
    )

    mainland_outline_path = find_file(
        root,
        [
            "04_analysis/C3_pipeline_rebuild/09_C3_biogeographic_concordance/"
            "10A_ecoregion_gis_audit/10A_mainland_outline_largest_component.gpkg",
            "04_analysis_USE _THIS/C3_pipeline_rebuild/09_C3_biogeographic_concordance/"
            "10A_ecoregion_gis_audit/10A_mainland_outline_largest_component.gpkg",
        ],
        "10A_mainland_outline_largest_component.gpkg",
    )

    incidence = pd.read_csv(incidence_path)
    lookup = pd.read_csv(lookup_path)
    trait = pd.read_csv(trait_path)
    crosswalk = pd.read_csv(crosswalk_path)
    summary10c = pd.read_csv(summary10c_path)

    required_inc = {"grid_cell_id", "genus"}
    required_trait = {"genus", "analysis_class"}
    required_lookup = {"grid_cell_id", "centroid_latitude_band"}
    required_cw = {
        "step10b_cell_id",
        "dominant_ecoregion",
        "primary_assignment_eligible",
    }

    for required, table, name in [
        (required_inc, incidence, "incidence"),
        (required_trait, trait, "trait"),
        (required_lookup, lookup, "lookup"),
        (required_cw, crosswalk, "crosswalk"),
    ]:
        missing = required - set(table.columns)
        if missing:
            raise RuntimeError(f"{name} missing columns: {sorted(missing)}")

    data = incidence.merge(
        trait[["genus", "analysis_class"]],
        on="genus",
        how="left",
        validate="many_to_one",
    )
    if data["analysis_class"].isna().any():
        missing = sorted(data.loc[data["analysis_class"].isna(), "genus"].unique())
        raise RuntimeError(f"Genera missing analysis_class: {missing[:30]}")

    classes = {"C3": "Ballooning (C3)", "N0": "Non-ballooning (N0)"}

    # Latitude groups and equal-cell size.
    bands = [b for b in BAND_ORDER if b in set(lookup["centroid_latitude_band"])]
    if len(bands) != 5:
        raise RuntimeError(f"Expected five retained latitude bands, found {bands}")

    band_cells = {
        band: lookup.loc[
            lookup["centroid_latitude_band"] == band, "grid_cell_id"
        ].dropna().astype(str).tolist()
        for band in bands
    }
    latitude_equal_n = min(len(v) for v in band_cells.values())

    # Prefer the formal 10 ecoregions from Step 10C.
    eco_rows = summary10c.copy()
    if "analysis_set" in eco_rows.columns:
        eco_rows = eco_rows[eco_rows["analysis_set"].astype(str) == "primary"]
    if "metric" in eco_rows.columns:
        eco_rows = eco_rows[eco_rows["metric"].astype(str) == "total_richness"]

    if "ecoregion" not in eco_rows.columns:
        raise RuntimeError("Step 10C summary does not contain 'ecoregion'.")

    ecoregions = eco_rows["ecoregion"].dropna().astype(str).drop_duplicates().tolist()

    eligible = crosswalk["primary_assignment_eligible"].astype(str).str.lower().isin(
        ["true", "t", "1", "yes"]
    )
    eco_cells = {
        eco: crosswalk.loc[
            eligible & (crosswalk["dominant_ecoregion"].astype(str) == eco),
            "step10b_cell_id",
        ].dropna().astype(str).tolist()
        for eco in ecoregions
    }
    ecoregions = [eco for eco in ecoregions if len(eco_cells[eco]) > 0]
    ecoregion_equal_n = min(len(eco_cells[eco]) for eco in ecoregions)

    # Per-class cell-to-genus maps.
    cellsets = {}
    support_rows = []
    for cls, label in classes.items():
        subset = data[data["analysis_class"] == cls]
        cellsets[cls] = (
            subset.groupby("grid_cell_id")["genus"].apply(lambda s: set(s.astype(str))).to_dict()
        )
        support_rows.append(
            {
                "class": cls,
                "label": label,
                "unique_genera": int(subset["genus"].nunique()),
                "occupied_cells": int(subset["grid_cell_id"].nunique()),
                "genus_cell_presences": int(len(subset)),
            }
        )

    pd.DataFrame(support_rows).to_csv(out / "11H_analysis_support_audit.csv", index=False)

    results = {}
    stability_rows = []
    pcoa_rows = []

    for scale, group_cells, labels, equal_n in [
        ("latitude", band_cells, bands, latitude_equal_n),
        ("ecoregion", eco_cells, ecoregions, ecoregion_equal_n),
    ]:
        for class_index, (cls, label) in enumerate(classes.items()):
            seed = args.seed + class_index + (100 if scale == "ecoregion" else 0)

            j_mats = resample_matrices(
                group_cells,
                labels,
                cellsets[cls],
                equal_n,
                jaccard,
                args.iterations,
                seed,
            )
            s_mats = resample_matrices(
                group_cells,
                labels,
                cellsets[cls],
                equal_n,
                simpson,
                args.iterations,
                seed,
            )

            j_med = save_matrix_outputs(
                j_mats, labels, out / f"11H_{scale}_{cls}_jaccard"
            )
            save_matrix_outputs(
                s_mats, labels, out / f"11H_{scale}_{cls}_simpson"
            )

            z = plot_dendrogram(
                j_med,
                labels,
                f"{label}: {scale} grouping (Jaccard)",
                out / f"11H_{scale}_{cls}_jaccard_dendrogram.png",
                horizontal=(scale == "ecoregion"),
            )
            p1, p2 = plot_pcoa(
                j_med,
                labels,
                f"{label}: {scale} composition (Jaccard PCoA)",
                out / f"11H_{scale}_{cls}_jaccard_pcoa.png",
            )
            nearest, distinct = grouping_stability(j_mats, labels)

            results[(scale, cls)] = {
                "jaccard_median": j_med,
                "nearest": nearest,
                "distinct": distinct,
                "leaf_order": [labels[i] for i in leaves_list(z)],
            }
            pcoa_rows.append(
                {
                    "scale": scale,
                    "class": cls,
                    "pcoa1_pct": p1,
                    "pcoa2_pct": p2,
                }
            )

            for rank, (a, b, freq) in enumerate(nearest, 1):
                stability_rows.append(
                    {
                        "scale": scale,
                        "class": cls,
                        "statistic": "nearest_pair_frequency",
                        "rank": rank,
                        "group_1": a,
                        "group_2": b,
                        "frequency": freq,
                    }
                )
            for rank, (group, freq) in enumerate(distinct, 1):
                stability_rows.append(
                    {
                        "scale": scale,
                        "class": cls,
                        "statistic": "most_distinct_group_frequency",
                        "rank": rank,
                        "group_1": group,
                        "group_2": "",
                        "frequency": freq,
                    }
                )

    pd.DataFrame(stability_rows).to_csv(
        out / "11H_grouping_stability_summary.csv", index=False
    )
    pd.DataFrame(pcoa_rows).to_csv(
        out / "11H_pcoa_variance_summary.csv", index=False
    )

    mantel_rows = []
    for scale in ["latitude", "ecoregion"]:
        c3 = results[(scale, "C3")]["jaccard_median"]
        n0 = results[(scale, "N0")]["jaccard_median"]
        rho, p = mantel_spearman(
            c3, n0, args.mantel_permutations, args.seed + 500
        )
        mantel_rows.append(
            {
                "scale": scale,
                "metric": "Jaccard",
                "spearman_mantel_rho": rho,
                "permutation_p_value": p,
                "permutations": args.mantel_permutations,
            }
        )
    pd.DataFrame(mantel_rows).to_csv(
        out / "11H_C3_vs_N0_distance_structure_similarity.csv", index=False
    )

    # Publication-style Baja map generated from the current Step 11H results.
    latitude_results = {
        "C3": results[("latitude", "C3")],
        "N0": results[("latitude", "N0")],
    }
    map_png, map_svg, map_pdf = build_baja_comparison_map(
        geojson_path=geojson_path,
        mainland_outline_path=mainland_outline_path,
        lookup=lookup,
        latitude_results=latitude_results,
        out=out,
        iterations=args.iterations,
        equal_cells=latitude_equal_n,
    )

    # Human-readable summary.
    support = {row["class"]: row for row in support_rows}
    lc3 = results[("latitude", "C3")]
    ln0 = results[("latitude", "N0")]
    ec3 = results[("ecoregion", "C3")]
    en0 = results[("ecoregion", "N0")]

    lines = [
        "STEP 11H — SEPARATE C3 AND N0 GEOGRAPHIC GROUPING",
        "==================================================",
        "",
        "PURPOSE",
        "Analyze ballooning-capable (C3) and fixed non-ballooning (N0) assemblages",
        "separately to test whether their geographic composition groups differently.",
        "",
        "DATA SUPPORT",
        f"C3: {support['C3']['unique_genera']} genera, "
        f"{support['C3']['occupied_cells']} occupied cells, "
        f"{support['C3']['genus_cell_presences']} genus-cell presences.",
        f"N0: {support['N0']['unique_genera']} genera, "
        f"{support['N0']['occupied_cells']} occupied cells, "
        f"{support['N0']['genus_cell_presences']} genus-cell presences.",
        "",
        "STANDARDIZATION",
        f"Latitude bands: {latitude_equal_n} occupied 25-km cells per band.",
        f"Ecoregions: {ecoregion_equal_n} primary-eligible cells per ecoregion.",
        f"Monte Carlo iterations: {args.iterations}.",
        "",
        "LATITUDE-BAND GROUPING",
        f"C3 nearest pair: {lc3['nearest'][0][0]} + {lc3['nearest'][0][1]} "
        f"({lc3['nearest'][0][2]*100:.1f}% of resamples).",
        f"N0 nearest pair: {ln0['nearest'][0][0]} + {ln0['nearest'][0][1]} "
        f"({ln0['nearest'][0][2]*100:.1f}% of resamples).",
        f"C3 most-distinct band: {lc3['distinct'][0][0]} "
        f"({lc3['distinct'][0][1]*100:.1f}% of resamples).",
        f"N0 most-distinct band: {ln0['distinct'][0][0]} "
        f"({ln0['distinct'][0][1]*100:.1f}% of resamples).",
        f"C3 clustering leaf order: {' -> '.join(lc3['leaf_order'])}",
        f"N0 clustering leaf order: {' -> '.join(ln0['leaf_order'])}",
        "",
        "ECOREGION GROUPING",
        f"C3 nearest pair: {ec3['nearest'][0][0]} + {ec3['nearest'][0][1]} "
        f"({ec3['nearest'][0][2]*100:.1f}% of resamples).",
        f"N0 nearest pair: {en0['nearest'][0][0]} + {en0['nearest'][0][1]} "
        f"({en0['nearest'][0][2]*100:.1f}% of resamples).",
        f"C3 most-distinct ecoregion: {ec3['distinct'][0][0]} "
        f"({ec3['distinct'][0][1]*100:.1f}% of resamples).",
        f"N0 most-distinct ecoregion: {en0['distinct'][0][0]} "
        f"({en0['distinct'][0][1]*100:.1f}% of resamples).",
        "",
        "C3 VS N0 DISTANCE-STRUCTURE SIMILARITY",
    ]
    for row in mantel_rows:
        lines.append(
            f"{row['scale'].capitalize()}: Mantel rho="
            f"{row['spearman_mantel_rho']:.3f}, "
            f"permutation p={row['permutation_p_value']:.4f}."
        )

    lines.extend(
        [
            "",
            "INTERPRETATION GUIDE",
            "A high nearest-pair frequency indicates a stable compositional cluster.",
            "A high most-distinct frequency indicates a region repeatedly isolated from",
            "the rest of that trait group's assemblage.",
            "Mantel rho measures whether C3 and N0 share the same overall distance pattern.",
            "",
            "STATUS",
            "Exploratory/supplementary. Review biological interpretation before inclusion",
            "in the public retained pipeline.",
        ]
    )

    # Generate a caption from the actual current-run results.
    c3_nearest = lc3["nearest"][0]
    n0_nearest = ln0["nearest"][0]
    c3_distinct = lc3["distinct"][0]
    n0_distinct = ln0["distinct"][0]
    latitude_mantel = next(
        row for row in mantel_rows if row["scale"] == "latitude"
    )

    caption = (
        "Figure X. Separate geographic organization of ballooning-capable (C3) "
        "and non-ballooning (N0) arachnid assemblages across Baja California. "
        f"Genus composition was analyzed independently for C3 and N0 assemblages "
        f"using Jaccard dissimilarity after standardizing each latitude band to "
        f"{latitude_equal_n} occupied 25-km cells over {args.iterations:,} Monte "
        "Carlo resamples. Brackets identify the pair of latitude bands most "
        "frequently recovered as compositionally closest, and hatched bands "
        "identify the region most frequently recovered as the most compositionally "
        "distinct within each trait group. "
        f"Ballooning assemblages most often grouped {c3_nearest[0]} with "
        f"{c3_nearest[1]} ({c3_nearest[2]*100:.1f}% of resamples) and most "
        f"frequently isolated {c3_distinct[0]} ({c3_distinct[1]*100:.1f}%), "
        f"whereas non-ballooning assemblages most often grouped {n0_nearest[0]} "
        f"with {n0_nearest[1]} ({n0_nearest[2]*100:.1f}%) and most frequently "
        f"isolated {n0_distinct[0]} ({n0_distinct[1]*100:.1f}%). "
        f"Overall C3 and N0 latitude-band distance structures had Mantel "
        f"rho = {latitude_mantel['spearman_mantel_rho']:.3f} and permutation "
        f"p = {latitude_mantel['permutation_p_value']:.4f}."
    )
    (out / "11H_BAJA_MAP_CAPTION.txt").write_text(caption + "\n", encoding="utf-8")

    lines.extend(
        [
            "",
            "PUBLICATION MAP",
            f"PNG: {map_png}",
            f"Editable SVG: {map_svg}",
            f"PDF: {map_pdf}",
            f"Caption: {out / '11H_BAJA_MAP_CAPTION.txt'}",
        ]
    )

    summary_path = out / "11H_RESULTS_SUMMARY.txt"
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    provenance = {
        "analysis": "Step 11H separate C3 and N0 geographic grouping",
        "project_root": str(root),
        "incidence_source": str(incidence_path),
        "trait_source": str(trait_path),
        "cell_lookup_source": str(lookup_path),
        "ecoregion_crosswalk_source": str(crosswalk_path),
        "step10c_summary_source": str(summary10c_path),
        "grid_geojson_source": str(geojson_path),
        "mainland_outline_source": str(mainland_outline_path),
        "iterations": args.iterations,
        "mantel_permutations": args.mantel_permutations,
        "seed": args.seed,
        "latitude_equal_cells": latitude_equal_n,
        "ecoregion_equal_cells": ecoregion_equal_n,
        "trait_definition": {
            "C3": "D1 + D2 + D3 ballooning",
            "N0": "fixed non-ballooning reference",
            "D4": "excluded",
        },
    }
    (out / "11H_PROVENANCE.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )

    print("\n".join(lines))
    print(f"\nOUTPUT_DIR={out}")
    print(f"SUMMARY={summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
