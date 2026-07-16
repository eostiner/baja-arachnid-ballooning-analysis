#!/usr/bin/env python3
"""
Step 10I — final map with Ballooning:Non-ballooning ratios by ecoregion
and latitude band.

This is a visualization-only synthesis of completed analyses.

Primary groups
--------------
Ballooning      = C3 = D1 + D2 + D3
Non-ballooning  = fixed N0 reference
D4              = excluded

Ecoregion ratios
----------------
Computed independently in every Step 10C equal-cell iteration:
    C3 pooled genus richness / N0 pooled genus richness
using 8 occupied 25-km cells per adequately sampled ecoregion.

Latitude-band ratios
--------------------
Read from Step 10G:
    C3 pooled genus richness / N0 pooled genus richness
using 22 occupied 25-km cells per latitude band and 5,000 iterations.

The figure does not introduce a new hypothesis test.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
import traceback
from pathlib import Path
from typing import Any

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, Rectangle, Wedge
from matplotlib import patheffects


EXPECTED_BANDS = ["23–24°N", "24–26°N", "26–28°N", "28–30°N", "30–32°N"]
FORMAL_MIN_CELLS = 8

ABBREVIATIONS = {
    "Bosques de la Sierra de la Laguna": "BSL",
    "Chaparral": "CH",
    "Costa Central del Golfo": "CGC",
    "Desierto Central": "CD",
    "Desierto de San Felipe": "SFD",
    "Desierto de Vizcaíno": "VD",
    "Matorral Costero": "CS",
    "Matorral Costero Rosetófilo": "RCS",
    "Matorrales Tropicales": "TS",
    "Planicies de Magdalena": "MP",
    "Selvas Bajas del Cabo": "CLDF",
    "Sierra de la Giganta": "SG",
    "Sierras de Juárez y San Pedro Mártir": "JSM",
}

ENGLISH_LABELS = {
    "Bosques de la Sierra de la Laguna": "Sierra de la Laguna Forests",
    "Chaparral": "Chaparral",
    "Costa Central del Golfo": "Central Gulf Coast",
    "Desierto Central": "Central Desert",
    "Desierto de San Felipe": "San Felipe Desert",
    "Desierto de Vizcaíno": "Vizcaíno Desert",
    "Matorral Costero": "Coastal Scrub",
    "Matorral Costero Rosetófilo": "Rosette Coastal Scrub",
    "Matorrales Tropicales": "Tropical Scrub",
    "Planicies de Magdalena": "Magdalena Plains",
    "Selvas Bajas del Cabo": "Cape Lowland Dry Forest",
    "Sierra de la Giganta": "Sierra de la Giganta",
    "Sierras de Juárez y San Pedro Mártir": "Juárez–San Pedro Mártir Ranges",
}

# Small offsets reduce overlap while retaining each symbol close to its ecoregion.
DONUT_OFFSETS = {
    "Matorral Costero": (-0.28, 0.24),
    "Chaparral": (0.22, 0.08),
    "Matorral Costero Rosetófilo": (-0.22, -0.03),
    "Desierto de San Felipe": (0.23, 0.07),
    "Desierto Central": (0.02, 0.10),
    "Desierto de Vizcaíno": (-0.22, 0.00),
    "Costa Central del Golfo": (0.27, 0.05),
    "Sierra de la Giganta": (0.20, 0.02),
    "Planicies de Magdalena": (-0.16, -0.04),
    "Matorrales Tropicales": (0.16, -0.05),
}

# Keep the same visual vocabulary as the preceding manuscript figures.
C_BALLOON = "#E67E22"
C_NON = "#2E9D38"
C_RATIO = "#2C7FB8"


def norm(x: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(x).lower())


def read_csv(path: Path) -> pd.DataFrame:
    last = None
    for enc in ("utf-8-sig", "utf-8", "latin1"):
        try:
            return pd.read_csv(path, encoding=enc, low_memory=False)
        except Exception as exc:
            last = exc
    raise RuntimeError(f"Could not read {path}: {last}")


def detect_col(df: pd.DataFrame, preferred: list[str], contains: str | None = None) -> str:
    lookup = {norm(c): str(c) for c in df.columns}
    for item in preferred:
        if norm(item) in lookup:
            return lookup[norm(item)]
    if contains:
        for col in df.columns:
            if contains in norm(col):
                return str(col)
    raise KeyError(f"Could not identify a required column from {preferred}; columns={list(df.columns)}")


def find_existing(candidates: list[Path], role: str) -> Path:
    for path in candidates:
        if path.exists() and path.is_file():
            return path
    raise FileNotFoundError(
        f"Could not locate {role}. Checked:\n" + "\n".join(str(path) for path in candidates)
    )


def md5_file(path: Path, block: int = 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(block)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def input_paths(root: Path) -> dict[str, Path]:
    bio = root / "04_analysis" / "C3_pipeline_rebuild" / "09_C3_biogeographic_concordance"
    return {
        "ecoregions": find_existing(
            [bio / "10A_ecoregion_gis_audit" / "10A_ecoregions_validated_mainland_only.gpkg"],
            "validated mainland ecoregion GIS",
        ),
        "breaks": find_existing(
            [
                bio / "10E_published_break_concordance" / "10E_break_registry_used.csv",
                bio / "10E_published_break_concordance" / "10E_PUBLISHED_BREAK_REGISTRY.csv",
            ],
            "published-break registry",
        ),
        "ecoregion_iterations": find_existing(
            [bio / "10C_equal_cell_ecoregion_richness" / "10C_equal_cell_iteration_results.csv"],
            "Step 10C equal-cell ecoregion iteration results",
        ),
        "ecoregion_sample_sizes": find_existing(
            [bio / "10C_equal_cell_ecoregion_richness" / "10C_ecoregion_sample_sizes.csv"],
            "Step 10C ecoregion sample sizes",
        ),
        "band_summaries": find_existing(
            [bio / "10G_band_ratio_synthesis" / "10G_band_ratio_summaries.csv"],
            "Step 10G band-ratio summaries",
        ),
        "band_tests": find_existing(
            [bio / "10G_band_ratio_synthesis" / "10G_adjacent_band_ratio_tests.csv"],
            "Step 10G adjacent-band tests",
        ),
    }


def quantile(series: pd.Series, probability: float) -> float:
    values = pd.to_numeric(series, errors="coerce").dropna().to_numpy(float)
    return float(np.quantile(values, probability)) if len(values) else np.nan


def load_data(paths: dict[str, Path]) -> dict[str, Any]:
    eco = gpd.read_file(paths["ecoregions"])
    if eco.crs is None:
        raise RuntimeError("Ecoregion layer has no coordinate reference system.")
    eco = eco.to_crs(4326)
    eco_col = detect_col(eco, ["ecoregion_label", "Nombre", "name"], "ecoregion")
    eco = eco.rename(columns={eco_col: "ecoregion"})
    eco["ecoregion"] = eco["ecoregion"].astype(str)

    breaks = read_csv(paths["breaks"])
    break_lat_col = detect_col(
        breaks, ["anchor_latitude", "latitude", "break_latitude"], "latitude"
    )
    breaks = breaks.rename(columns={break_lat_col: "anchor_latitude"})
    breaks["anchor_latitude"] = pd.to_numeric(breaks["anchor_latitude"], errors="coerce")
    breaks = breaks.dropna(subset=["anchor_latitude"]).sort_values("anchor_latitude")

    iterations = read_csv(paths["ecoregion_iterations"])
    required = {
        "analysis_set",
        "iteration",
        "ecoregion",
        "equal_cells",
        "C3_richness",
        "N0_richness",
    }
    missing = required - set(iterations.columns)
    if missing:
        raise RuntimeError(f"Step 10C iteration table is missing columns: {sorted(missing)}")
    iterations = iterations[
        iterations["analysis_set"].astype(str).str.lower().eq("primary")
    ].copy()
    for col in ["equal_cells", "C3_richness", "N0_richness"]:
        iterations[col] = pd.to_numeric(iterations[col], errors="coerce")
    iterations = iterations.dropna(subset=["ecoregion", "C3_richness", "N0_richness"])
    if iterations.empty:
        raise RuntimeError("No primary Step 10C ecoregion iterations were found.")

    sample_sizes = read_csv(paths["ecoregion_sample_sizes"])
    sample_eco_col = detect_col(
        sample_sizes, ["ecoregion", "ecoregion_label", "dominant_ecoregion"], "ecoregion"
    )
    sample_sizes = sample_sizes.rename(columns={sample_eco_col: "ecoregion"})
    n_col = detect_col(
        sample_sizes,
        ["n_primary_cells", "primary_cells", "eligible_primary_cells"],
        "primarycells",
    )
    sample_sizes = sample_sizes.rename(columns={n_col: "n_primary_cells"})
    sample_sizes["n_primary_cells"] = pd.to_numeric(
        sample_sizes["n_primary_cells"], errors="coerce"
    )

    band = read_csv(paths["band_summaries"])
    required_band = {
        "latitude_band",
        "median_log2_ratio",
        "log2_ratio_q025",
        "log2_ratio_q975",
        "median_B_to_N_ratio",
        "mean_ballooning_share",
        "equal_cells",
        "iterations",
    }
    missing_band = required_band - set(band.columns)
    if missing_band:
        raise RuntimeError(f"Step 10G band summary is missing columns: {sorted(missing_band)}")
    for col in required_band - {"latitude_band"}:
        band[col] = pd.to_numeric(band[col], errors="coerce")
    present = set(band["latitude_band"].astype(str))
    missing_expected = [name for name in EXPECTED_BANDS if name not in present]
    if missing_expected:
        raise RuntimeError(f"Missing expected latitude bands: {missing_expected}")

    tests = read_csv(paths["band_tests"])
    for col in [
        "boundary_latitude",
        "random_split_two_sided_p",
        "BH_q_four_adjacent_boundaries",
    ]:
        if col in tests.columns:
            tests[col] = pd.to_numeric(tests[col], errors="coerce")

    return {
        "ecoregions": eco,
        "breaks": breaks,
        "iterations": iterations,
        "sample_sizes": sample_sizes,
        "band": band,
        "tests": tests,
    }


def summarize_ecoregions(data: dict[str, Any]) -> pd.DataFrame:
    frame = data["iterations"].copy()
    # No zero counts occur in the frozen analysis, but the correction keeps the
    # visualization defined if a future sensitivity subset includes a zero.
    frame["zero_correction_used"] = (frame["C3_richness"] <= 0) | (frame["N0_richness"] <= 0)
    frame["B_to_N_ratio"] = np.where(
        frame["zero_correction_used"],
        (frame["C3_richness"] + 0.5) / (frame["N0_richness"] + 0.5),
        frame["C3_richness"] / frame["N0_richness"],
    )
    frame["log2_ratio"] = np.log2(frame["B_to_N_ratio"])
    frame["classified_richness"] = frame["C3_richness"] + frame["N0_richness"]
    frame["ballooning_share"] = np.where(
        frame["classified_richness"] > 0,
        frame["C3_richness"] / frame["classified_richness"],
        np.nan,
    )

    rows: list[dict[str, Any]] = []
    for name, group in frame.groupby("ecoregion", sort=True):
        rows.append(
            {
                "ecoregion": str(name),
                "iterations": int(group["iteration"].nunique()),
                "equal_cells": int(group["equal_cells"].dropna().iloc[0]),
                "mean_ballooning_richness": float(group["C3_richness"].mean()),
                "mean_non_ballooning_richness": float(group["N0_richness"].mean()),
                "mean_classified_richness": float(group["classified_richness"].mean()),
                "median_B_to_N_ratio": float(group["B_to_N_ratio"].median()),
                "ratio_q025": quantile(group["B_to_N_ratio"], 0.025),
                "ratio_q975": quantile(group["B_to_N_ratio"], 0.975),
                "median_log2_ratio": float(group["log2_ratio"].median()),
                "log2_ratio_q025": quantile(group["log2_ratio"], 0.025),
                "log2_ratio_q975": quantile(group["log2_ratio"], 0.975),
                "mean_ballooning_share": float(group["ballooning_share"].mean()),
                "zero_correction_fraction": float(group["zero_correction_used"].mean()),
            }
        )
    summary = pd.DataFrame(rows).merge(
        data["sample_sizes"][["ecoregion", "n_primary_cells"]],
        on="ecoregion",
        how="left",
    )
    summary["formal_comparison"] = summary["n_primary_cells"] >= summary["equal_cells"]
    return summary


def ratio_tick_label(value: float) -> str:
    if abs(value) < 1e-10:
        return "1:1"
    if value > 0:
        return f"{2 ** value:.2g}:1"
    return f"1:{2 ** (-value):.2g}"


def draw_donut(
    ax,
    x: float,
    y: float,
    total: float,
    ballooning: float,
    non_ballooning: float,
    max_total: float,
) -> float:
    radius = 0.145 + 0.205 * math.sqrt(max(total, 0.0) / max_total)
    values = [max(ballooning, 0.0), max(non_ballooning, 0.0)]
    denominator = sum(values)
    start = 90.0
    for value, color in zip(values, [C_BALLOON, C_NON]):
        angle = 360.0 * value / denominator if denominator > 0 else 0.0
        ax.add_patch(
            Wedge(
                (x, y),
                radius,
                start,
                start + angle,
                facecolor=color,
                edgecolor="white",
                linewidth=0.7,
                zorder=8,
            )
        )
        start += angle
    ax.add_patch(
        Circle(
            (x, y),
            radius,
            facecolor="none",
            edgecolor="0.12",
            linewidth=0.8,
            zorder=9,
        )
    )
    ax.add_patch(
        Circle(
            (x, y),
            radius * 0.43,
            facecolor="white",
            edgecolor="0.25",
            linewidth=0.45,
            zorder=10,
        )
    )
    ax.text(
        x,
        y,
        f"{total:.0f}",
        ha="center",
        va="center",
        fontsize=7.1,
        fontweight="bold",
        zorder=11,
    )
    return radius


def draw_map(ax, data: dict[str, Any], eco_summary: pd.DataFrame) -> None:
    eco = data["ecoregions"].copy()
    labels = sorted(eco["ecoregion"].unique())
    cmap = plt.get_cmap("tab20")
    polygon_colors = {name: cmap(index % 20) for index, name in enumerate(labels)}

    summary_lookup = eco_summary.set_index("ecoregion")
    for name in labels:
        polygons = eco[eco["ecoregion"] == name]
        formal = (
            name in summary_lookup.index
            and bool(summary_lookup.loc[name, "formal_comparison"])
        )
        polygons.plot(
            ax=ax,
            facecolor=polygon_colors[name] if formal else (0.91, 0.91, 0.91, 1.0),
            edgecolor="0.25",
            linewidth=0.55,
            hatch=None if formal else "///",
            zorder=1,
        )

    for latitude in data["breaks"]["anchor_latitude"]:
        latitude = float(latitude)
        ax.axhline(
            latitude,
            color="0.14",
            linestyle=(0, (5, 3)),
            linewidth=1.0,
            zorder=4,
        )
        ax.text(
            -117.20,
            latitude + 0.04,
            f"{latitude:.0f}°N",
            fontsize=7.1,
            fontweight="bold",
            ha="center",
            va="bottom",
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="0.3", lw=0.6),
            zorder=12,
        )

    formal_summary = eco_summary[eco_summary["formal_comparison"]].copy()
    max_total = float(formal_summary["mean_classified_richness"].max())

    # One representative point per dissolved ecoregion.
    dissolved = eco.dissolve(by="ecoregion")
    points = dissolved.representative_point()

    for row in formal_summary.itertuples(index=False):
        name = str(row.ecoregion)
        if name not in points.index:
            continue
        point = points.loc[name]
        dx, dy = DONUT_OFFSETS.get(name, (0.0, 0.0))
        x, y = float(point.x) + dx, float(point.y) + dy
        radius = draw_donut(
            ax,
            x,
            y,
            float(row.mean_classified_richness),
            float(row.mean_ballooning_richness),
            float(row.mean_non_ballooning_richness),
            max_total,
        )
        ratio_text = f"{float(row.median_B_to_N_ratio):.2f}:1"
        label = f"{ABBREVIATIONS.get(name, name)}  {ratio_text}"
        text_y = y - radius - 0.055
        text = ax.text(
            x,
            text_y,
            label,
            ha="center",
            va="top",
            fontsize=6.8,
            fontweight="bold",
            zorder=13,
        )
        text.set_path_effects(
            [patheffects.withStroke(linewidth=2.6, foreground="white")]
        )

    ax.set_xlim(-117.6, -108.9)
    ax.set_ylim(22.65, 33.15)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(color="0.90", linewidth=0.45, zorder=0)

    # Compact visual legend in the northeastern ocean margin.
    x0, y0 = -110.05, 32.75
    ax.add_patch(
        Rectangle(
            (x0 - 0.18, y0 - 1.60),
            1.30,
            1.50,
            facecolor="white",
            edgecolor="0.40",
            linewidth=0.7,
            alpha=0.95,
            zorder=20,
        )
    )
    ax.add_patch(
        Rectangle(
            (x0, y0 - 0.25),
            0.18,
            0.16,
            facecolor=C_BALLOON,
            edgecolor="0.2",
            zorder=21,
        )
    )
    ax.text(
        x0 + 0.24,
        y0 - 0.17,
        "Ballooning",
        fontsize=7.2,
        va="center",
        zorder=21,
    )
    ax.add_patch(
        Rectangle(
            (x0, y0 - 0.51),
            0.18,
            0.16,
            facecolor=C_NON,
            edgecolor="0.2",
            zorder=21,
        )
    )
    ax.text(
        x0 + 0.24,
        y0 - 0.43,
        "Non-ballooning",
        fontsize=7.2,
        va="center",
        zorder=21,
    )
    ax.text(
        x0,
        y0 - 0.69,
        "Donut area = mean classified\nrichness from 8 occupied cells",
        fontsize=6.65,
        va="top",
        zorder=21,
    )
    ax.text(
        x0,
        y0 - 1.05,
        "Center = C3 + N0 richness\nLabel = median B:N ratio",
        fontsize=6.65,
        va="top",
        zorder=21,
    )
    ax.text(
        x0,
        y0 - 1.39,
        "D4 excluded",
        fontsize=6.65,
        va="top",
        fontweight="bold",
        zorder=21,
    )


def draw_band_inset(ax, band: pd.DataFrame, tests: pd.DataFrame) -> None:
    order = list(reversed(EXPECTED_BANDS))
    table = (
        band.assign(latitude_band=band["latitude_band"].astype(str))
        .set_index("latitude_band")
        .loc[order]
        .reset_index()
    )
    y = np.arange(len(order))
    med = table["median_log2_ratio"].to_numpy(float)
    low = table["log2_ratio_q025"].to_numpy(float)
    high = table["log2_ratio_q975"].to_numpy(float)

    ax.axvline(0.0, color="0.25", linestyle="--", linewidth=1.0, zorder=1)
    ax.errorbar(
        med,
        y,
        xerr=[med - low, high - med],
        fmt="o",
        color=C_RATIO,
        ecolor=C_RATIO,
        linewidth=1.35,
        capsize=3.4,
        markersize=5.8,
        zorder=4,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=8.2)

    xmin = min(float(low.min()) - 0.35, -2.25)
    xmax = max(float(high.max()) + 0.70, 0.95)
    ax.set_xlim(xmin, xmax)
    ticks = np.arange(math.floor(xmin), math.ceil(xmax) + 1, 1.0)
    ax.set_xticks(ticks)
    ax.set_xticklabels([ratio_tick_label(value) for value in ticks], fontsize=7.6)
    ax.set_xlabel("Ballooning : non-ballooning ratio\n(log₂-spaced axis)", fontsize=8.2)
    ax.set_title(
        "Latitude-band balance\n22 occupied cells per band",
        fontsize=10.3,
        fontweight="bold",
        loc="left",
        pad=7,
    )
    ax.grid(axis="x", color="0.90", linewidth=0.55)

    for index, row in table.iterrows():
        ratio = float(row["median_B_to_N_ratio"])
        share = 100.0 * float(row["mean_ballooning_share"])
        position = min(float(row["log2_ratio_q975"]) + 0.10, xmax - 0.03)
        horizontal = "left"
        if position >= xmax - 0.18:
            position = xmax - 0.03
            horizontal = "right"
        ax.text(
            position,
            index,
            f"{ratio:.2f}:1  ({share:.0f}% B)",
            ha=horizontal,
            va="center",
            fontsize=7.5,
        )



def ecoregion_key() -> str:
    ordered = ["CS", "CH", "RCS", "SFD", "CD", "VD", "CGC", "SG", "MP", "TS"]
    reverse = {short: full for full, short in ABBREVIATIONS.items()}
    items = [
        f"{short} = {ENGLISH_LABELS.get(reverse[short], reverse[short])}"
        for short in ordered
        if short in reverse
    ]
    lines = []
    for start in range(0, len(items), 2):
        lines.append("     ".join(items[start : start + 2]))
    return "\n".join(lines)


def make_figure(
    out_dir: Path,
    data: dict[str, Any],
    eco_summary: pd.DataFrame,
    dpi: int,
) -> Path:
    publication = out_dir / "publication_outputs"
    publication.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(12.8, 10.8))
    ax_map = fig.add_axes([0.055, 0.115, 0.66, 0.79])
    ax_band = fig.add_axes([0.755, 0.585, 0.22, 0.27])

    draw_map(ax_map, data, eco_summary)
    draw_band_inset(ax_band, data["band"], data["tests"])

    fig.suptitle(
        "Baja arachnid dispersal-group balance",
        fontsize=17.0,
        fontweight="bold",
        y=0.965,
    )
    fig.text(
        0.5,
        0.935,
        "Equal-cell Ballooning (C3 = D1–D3) versus fixed Non-ballooning (N0); D4 excluded",
        ha="center",
        fontsize=9.4,
    )

    qvalues = pd.to_numeric(
        data["tests"].get(
            "BH_q_four_adjacent_boundaries", pd.Series(dtype=float)
        ),
        errors="coerce",
    )
    minimum_q = float(qvalues.min()) if qvalues.notna().any() else np.nan
    if np.isfinite(minimum_q) and minimum_q >= 0.05:
        test_note = "No adjacent-band shift remained significant after BH correction."
    elif np.isfinite(minimum_q):
        test_note = f"Minimum adjacent-band BH q = {minimum_q:.3f}."
    else:
        test_note = "Adjacent-band corrected tests unavailable."

    fig.text(
        0.755,
        0.515,
        "Statistical basis",
        fontsize=9.3,
        fontweight="bold",
        ha="left",
    )
    fig.text(
        0.755,
        0.492,
        "Map labels: median B:N ratio across 5,000 draws of 8 occupied\n"
        "cells per ecoregion. Inset bars: 2.5th–97.5th percentiles across\n"
        "5,000 draws of 22 occupied cells per band.\n\n"
        + test_note
        + "\nIntervals describe equal-cell resampling uncertainty; they are not\n"
        "separate significance tests against a 1:1 ratio.",
        fontsize=7.25,
        ha="left",
        va="top",
        color="0.30",
        bbox=dict(boxstyle="round,pad=0.45", fc="white", ec="0.75", lw=0.7),
    )

    fig.text(
        0.055,
        0.043,
        ecoregion_key(),
        ha="left",
        va="bottom",
        fontsize=6.8,
        color="0.26",
        linespacing=1.35,
    )
    fig.text(
        0.98,
        0.043,
        "Dashed lines mark a priori test latitudes at 24°, 26°, 28°, and 30°N.\n"
        "Sparse ecoregions (<8 eligible cells) are hatched and not assigned ratio symbols.",
        ha="right",
        va="bottom",
        fontsize=6.9,
        color="0.34",
    )

    png = publication / "10I_final_map_ecoregion_and_band_ratios.png"
    pdf = publication / "10I_final_map_ecoregion_and_band_ratios.pdf"
    svg = publication / "10I_final_map_ecoregion_and_band_ratios.svg"
    fig.savefig(png, dpi=dpi, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    return png


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root")
    parser.add_argument("--dpi", type=int, default=400)
    args = parser.parse_args()

    root = Path(args.project_root).expanduser().resolve()
    out_dir = (
        root
        / "04_analysis"
        / "C3_pipeline_rebuild"
        / "09_C3_biogeographic_concordance"
        / "10I_final_map_and_ratios"
    )
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        paths = input_paths(root)
        data = load_data(paths)
        eco_summary = summarize_ecoregions(data)
        eco_summary.to_csv(out_dir / "10I_ecoregion_ratio_summary.csv", index=False)
        data["band"].to_csv(out_dir / "10I_band_ratio_summary.csv", index=False)
        data["tests"].to_csv(out_dir / "10I_adjacent_band_ratio_tests.csv", index=False)

        manifest = pd.DataFrame(
            [
                {
                    "role": role,
                    "file": str(path),
                    "md5": md5_file(path),
                    "size_bytes": path.stat().st_size,
                }
                for role, path in paths.items()
            ]
        )
        manifest.to_csv(out_dir / "10I_input_manifest.csv", index=False)

        figure = make_figure(out_dir, data, eco_summary, args.dpi)

        qvalues = pd.to_numeric(
            data["tests"].get(
                "BH_q_four_adjacent_boundaries", pd.Series(dtype=float)
            ),
            errors="coerce",
        )
        summary = {
            "step": "10I",
            "audit_status": "PASS_VISUAL_SYNTHESIS_ONLY",
            "primary_groups": {
                "ballooning": "C3 = D1 + D2 + D3",
                "non_ballooning": "fixed N0",
                "D4": "excluded",
            },
            "ecoregion_equal_cells": int(eco_summary["equal_cells"].dropna().iloc[0]),
            "band_equal_cells": int(data["band"]["equal_cells"].dropna().iloc[0]),
            "iterations": int(data["band"]["iterations"].dropna().iloc[0]),
            "formal_ecoregions": int(eco_summary["formal_comparison"].sum()),
            "minimum_adjacent_band_BH_q": (
                float(qvalues.min()) if qvalues.notna().any() else None
            ),
            "figure": str(figure),
        }
        (out_dir / "10I_analysis_summary.json").write_text(
            json.dumps(summary, indent=2), encoding="utf-8"
        )

        caption = (
            "Figure X. Equal-cell Ballooning:Non-ballooning genus-richness balance "
            "across Baja California ecoregions and latitude bands. Donut area is "
            "proportional to mean classified richness (C3 + N0), and wedges show "
            "the relative contributions of Ballooning genera (C3 = D1–D3) and the "
            "fixed Non-ballooning reference (N0), calculated from 5,000 draws of "
            "eight occupied 25-km cells per adequately sampled ecoregion. Labels "
            "beside donuts give the median Ballooning:Non-ballooning richness ratio. "
            "The inset shows median ratios and 2.5th–97.5th percentile subset-resampling "
            "intervals from 5,000 draws of 22 occupied cells per latitude band. "
            "D4 genera were excluded. Dashed lines mark a priori test "
            "latitudes corresponding to commonly discussed biogeographic transitions; "
            "they are not inferred barriers. Subset-resampling intervals describe equal-cell "
            "resampling uncertainty and are not separate significance tests against "
            "a 1:1 ratio; adjacent-band changes were evaluated with local random-split "
            "tests and Benjamini–Hochberg correction."
        )
        (out_dir / "10I_FIGURE_CAPTION_DRAFT.txt").write_text(
            caption + "\n", encoding="utf-8"
        )
        (out_dir / "README_RESULTS_FIRST.txt").write_text(
            "STEP 10I COMPLETE\n\n"
            "Primary figure:\n"
            "  publication_outputs/Figure_3_Biogeographic_Dispersal_Balance.png\n"
            "  publication_outputs/Figure_3_Biogeographic_Dispersal_Balance.pdf\n"
            "  publication_outputs/Figure_3_Biogeographic_Dispersal_Balance.svg\n"
            "  publication_outputs/Figure_3_Biogeographic_Dispersal_Balance.jpg\n\n"
            "Supporting tables:\n"
            "  10I_ecoregion_ratio_summary.csv\n"
            "  10I_band_ratio_summary.csv\n"
            "  10I_adjacent_band_ratio_tests.csv\n"
            "  10I_FIGURE_CAPTION_DRAFT.txt\n\n"
            "This step is visualization only and does not add a new hypothesis test.\n",
            encoding="utf-8",
        )

        print("STEP 10I COMPLETE")
        print("AUDIT_STATUS=PASS_VISUAL_SYNTHESIS_ONLY")
        print(f"OUTPUT_DIR={out_dir}")
        print(f"FIGURE={figure}")
        print(f"FORMAL_ECOREGIONS={int(eco_summary['formal_comparison'].sum())}")
        return 0

    except Exception as exc:
        (out_dir / "10I_FAILURE.txt").write_text(
            f"STEP 10I FAILED: {exc}\n\n{traceback.format_exc()}",
            encoding="utf-8",
        )
        print(f"STEP 10I FAILED: {exc}")
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
