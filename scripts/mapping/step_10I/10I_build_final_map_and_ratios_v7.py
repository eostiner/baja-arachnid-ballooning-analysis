#!/usr/bin/env python3
"""Step 10I v7.1 — integrated map, latitude-band ratios, and legend.

This wrapper changes only the figure layout. All data validation, equal-cell
summaries, ratios, intervals, and adjacent-band tests remain those of the
validated Step 10I implementation.
"""
from __future__ import annotations

import importlib.util
import math
import shutil
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import patheffects
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle, Wedge

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("step10i_base", HERE / "_base_v2.py")
base = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(base)

plt.style.use("default")
mpl.rcParams.update({
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.edgecolor": "white",
    "axes.facecolor": "white",
    "axes.edgecolor": "#333333",
    "axes.labelcolor": "#222222",
    "text.color": "#222222",
    "xtick.color": "#222222",
    "ytick.color": "#222222",
    "font.size": 10.5,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

C_BALLOON = "#E67E22"
C_NON = "#2E9D38"
C_RATIO = "#2C7FB8"
C_BREAK = "#555555"

ECO_COLORS = {
    "Matorral Costero": "#B8D8E8",
    "Chaparral": "#D9E7A8",
    "Matorral Costero Rosetófilo": "#F5C58F",
    "Desierto de San Felipe": "#F0D992",
    "Desierto Central": "#F6B877",
    "Desierto de Vizcaíno": "#C9B18C",
    "Costa Central del Golfo": "#F5D0DE",
    "Sierra de la Giganta": "#C7B6D8",
    "Planicies de Magdalena": "#EEE7A8",
    "Matorrales Tropicales": "#F1A89F",
    "Bosques de la Sierra de la Laguna": "#D6C7E8",
    "Selvas Bajas del Cabo": "#F4B0A4",
    "Sierras de Juárez y San Pedro Mártir": "#D8D8D8",
}

DONUT_OFFSETS = {
    "Matorral Costero": (-0.38, 0.23),
    "Chaparral": (0.34, 0.16),
    "Matorral Costero Rosetófilo": (-0.35, -0.04),
    "Desierto de San Felipe": (0.38, 0.17),
    "Desierto Central": (-0.04, 0.14),
    "Desierto de Vizcaíno": (-0.34, 0.03),
    "Costa Central del Golfo": (0.36, 0.08),
    "Sierra de la Giganta": (0.29, 0.00),
    "Planicies de Magdalena": (-0.23, -0.04),
    "Matorrales Tropicales": (0.25, -0.04),
}

BAND_MIDPOINTS = {
    "23–24°N": 23.50,
    "24–26°N": 25.00,
    "26–28°N": 27.00,
    "28–30°N": 29.00,
    "30–32°N": 31.00,
}


def draw_donut(ax, x, y, total, ballooning, non_ballooning, max_total):
    radius = 0.18 + 0.22 * math.sqrt(max(total, 0.0) / max_total)
    values = [max(ballooning, 0.0), max(non_ballooning, 0.0)]
    denominator = sum(values)
    start = 90.0
    for value, color in zip(values, [C_BALLOON, C_NON]):
        angle = 360.0 * value / denominator if denominator > 0 else 0.0
        ax.add_patch(Wedge(
            (x, y), radius, start, start + angle,
            facecolor=color, edgecolor="white", linewidth=1.15, zorder=8,
        ))
        start += angle
    ax.add_patch(Circle(
        (x, y), radius, facecolor="none", edgecolor="#222222",
        linewidth=1.0, zorder=9,
    ))
    ax.add_patch(Circle(
        (x, y), radius * 0.44, facecolor="white", edgecolor="#555555",
        linewidth=0.55, zorder=10,
    ))
    ax.text(
        x, y, f"{total:.0f}", ha="center", va="center",
        fontsize=9.0, fontweight="bold", zorder=11,
    )
    return radius


def draw_symbol_legend(ax):
    x0, y0 = -117.45, 22.95
    width, height = 2.60, 2.35
    ax.add_patch(FancyBboxPatch(
        (x0, y0), width, height,
        boxstyle="round,pad=0.06,rounding_size=0.05",
        facecolor="white", edgecolor="#666666", linewidth=0.8,
        alpha=0.98, zorder=30,
    ))
    ax.text(
        x0 + 0.14, y0 + height - 0.17, "How to read the map",
        fontsize=9.4, fontweight="bold", va="top", zorder=31,
    )

    ax.add_patch(Rectangle(
        (x0 + 0.16, y0 + height - 0.50), 0.20, 0.17,
        facecolor=C_BALLOON, edgecolor="#333333", linewidth=0.7, zorder=31,
    ))
    ax.text(
        x0 + 0.44, y0 + height - 0.415,
        "Ballooning (C3 = D1–D3)",
        fontsize=7.8, va="center", zorder=31,
    )
    ax.add_patch(Rectangle(
        (x0 + 0.16, y0 + height - 0.78), 0.20, 0.17,
        facecolor=C_NON, edgecolor="#333333", linewidth=0.7, zorder=31,
    ))
    ax.text(
        x0 + 0.44, y0 + height - 0.695,
        "Non-ballooning (fixed N0)",
        fontsize=7.8, va="center", zorder=31,
    )

    # Fixed-size example donut.
    cx, cy, radius = x0 + 0.36, y0 + 1.03, 0.25
    ax.add_patch(Wedge(
        (cx, cy), radius, 90, 235,
        facecolor=C_BALLOON, edgecolor="white", linewidth=1.0, zorder=31,
    ))
    ax.add_patch(Wedge(
        (cx, cy), radius, 235, 450,
        facecolor=C_NON, edgecolor="white", linewidth=1.0, zorder=31,
    ))
    ax.add_patch(Circle(
        (cx, cy), radius, facecolor="none", edgecolor="#222222",
        linewidth=0.9, zorder=32,
    ))
    ax.add_patch(Circle(
        (cx, cy), radius * 0.44, facecolor="white", edgecolor="#555555",
        linewidth=0.5, zorder=33,
    ))
    ax.text(cx, cy, "40", ha="center", va="center",
            fontsize=7.6, fontweight="bold", zorder=34)

    ax.text(
        x0 + 0.72, y0 + 1.29,
        "Donut area = equal-cell C3 + N0 richness",
        fontsize=7.45, va="center", zorder=31,
    )
    ax.text(
        x0 + 0.72, y0 + 1.03,
        "Center = expected classified richness",
        fontsize=7.45, va="center", zorder=31,
    )
    ax.text(
        x0 + 0.72, y0 + 0.77,
        "Text below donut = median B:N ratio",
        fontsize=7.45, va="center", zorder=31,
    )

    # Keep the two line-symbol explanations on separate rows so the public-
    # facing wording remains legible in both raster and vector exports.
    ax.plot(
        [x0 + 0.16, x0 + 0.48], [y0 + 0.45, y0 + 0.55],
        color=C_BREAK, linestyle=(0, (5, 4)), linewidth=1.0, zorder=31,
    )
    ax.text(
        x0 + 0.56, y0 + 0.45,
        "A priori test latitude",
        fontsize=7.35, va="center", zorder=31,
    )
    ax.plot(
        [x0 + 0.16, x0 + 0.51], [y0 + 0.22, y0 + 0.32],
        color=C_RATIO, linewidth=1.7, zorder=31,
    )
    ax.scatter(
        [x0 + 0.335], [y0 + 0.32], s=26, color=C_RATIO,
        edgecolor="white", linewidth=0.6, zorder=32,
    )
    ax.text(
        x0 + 0.56, y0 + 0.22,
        "Band median + subset interval",
        fontsize=7.10, va="center", zorder=31,
    )
    ax.text(
        x0 + 0.16, y0 + 0.035,
        f"D4 excluded. Hatched ecoregions had <{base.FORMAL_MIN_CELLS} eligible cells.",
        fontsize=6.85, va="bottom", color="#444444", zorder=31,
    )


def draw_integrated_band_ratios(ax, band, tests):
    table = band.copy()
    table["latitude_band"] = table["latitude_band"].astype(str)
    table = table.set_index("latitude_band").loc[base.EXPECTED_BANDS].reset_index()

    panel_left = -108.95
    panel_right = -106.55
    scale_left = -108.18
    scale_right = -107.10
    label_x = -108.80
    value_x = -106.67
    log_min, log_max = -2.25, 0.85

    ax.add_patch(FancyBboxPatch(
        (panel_left, 22.82), panel_right - panel_left, 9.82,
        boxstyle="round,pad=0.06,rounding_size=0.05",
        facecolor="white", edgecolor="#777777", linewidth=0.8,
        alpha=0.98, zorder=25,
    ))
    ax.text(
        panel_left + 0.13, 32.43, "Latitude-band balance",
        fontsize=10.0, fontweight="bold", ha="left", va="top", zorder=26,
    )
    ax.text(
        panel_left + 0.13, 32.15,
        "Median B:N ratio and 2.5th–97.5th percentile interval",
        fontsize=7.25, ha="left", va="top", color="#444444", zorder=26,
    )
    band_equal_cells = int(pd.to_numeric(table["equal_cells"], errors="coerce").dropna().iloc[0])
    band_iterations = int(pd.to_numeric(table["iterations"], errors="coerce").dropna().iloc[0])
    ax.text(
        panel_left + 0.13, 31.92,
        f"{band_equal_cells} occupied cells per band; {band_iterations:,} draws",
        fontsize=7.25, ha="left", va="top", color="#444444", zorder=26,
    )

    def map_x(value):
        return scale_left + (value - log_min) / (log_max - log_min) * (scale_right - scale_left)

    x_equal = map_x(0.0)
    ax.plot(
        [x_equal, x_equal], [23.35, 31.55],
        color="#4A4A4A", linestyle="--", linewidth=1.0, zorder=26,
    )
    ax.text(
        x_equal, 31.73, "1:1", fontsize=6.9, ha="center",
        va="bottom", color="#333333", fontweight="bold", zorder=26,
    )
    ax.text(
        panel_left + 0.16, 31.48, "higher N0", fontsize=6.3, ha="left",
        va="bottom", color=C_NON, fontweight="bold", zorder=26,
    )
    ax.text(
        panel_right - 0.16, 31.48, "higher C3", fontsize=6.3, ha="right",
        va="bottom", color=C_BALLOON, fontweight="bold", zorder=26,
    )

    for index, row in enumerate(table.itertuples(index=False)):
        band_name = str(row.latitude_band)
        y = BAND_MIDPOINTS[band_name]
        med = float(row.median_log2_ratio)
        low = float(row.log2_ratio_q025)
        high = float(row.log2_ratio_q975)
        ratio = float(row.median_B_to_N_ratio)
        share = 100.0 * float(row.mean_ballooning_share)

        if index % 2 == 0:
            ax.add_patch(Rectangle(
                (panel_left + 0.05, y - 0.57),
                panel_right - panel_left - 0.10,
                1.14,
                facecolor="#FAFAFA", edgecolor="none", zorder=25,
            ))
        ax.plot(
            [panel_left + 0.08, panel_right - 0.08], [y - 0.57, y - 0.57],
            color="#E6E6E6", linewidth=0.65, zorder=26,
        )
        ax.text(
            label_x, y + 0.22, band_name,
            fontsize=8.0, fontweight="bold", ha="left", va="center", zorder=27,
        )
        ax.plot(
            [map_x(low), map_x(high)], [y, y],
            color=C_RATIO, linewidth=2.0, zorder=27,
        )
        ax.plot(
            [map_x(low), map_x(low)], [y - 0.11, y + 0.11],
            color=C_RATIO, linewidth=1.2, zorder=27,
        )
        ax.plot(
            [map_x(high), map_x(high)], [y - 0.11, y + 0.11],
            color=C_RATIO, linewidth=1.2, zorder=27,
        )
        ax.scatter(
            [map_x(med)], [y], s=40, color=C_RATIO,
            edgecolor="white", linewidth=0.7, zorder=28,
        )
        ax.text(
            value_x, y + 0.11, f"{ratio:.2f}:1",
            fontsize=8.0, fontweight="bold", ha="right", va="center", zorder=27,
        )
        ax.text(
            value_x, y - 0.17, f"{share:.0f}% ballooning",
            fontsize=6.9, ha="right", va="center", color="#444444", zorder=27,
        )

    qvalues = pd.to_numeric(
        tests.get("BH_q_four_adjacent_boundaries", pd.Series(dtype=float)),
        errors="coerce",
    )
    min_q = float(qvalues.min()) if qvalues.notna().any() else np.nan
    if np.isfinite(min_q) and min_q >= 0.05:
        note = "No adjacent-band ratio shift remained significant\nafter BH correction."
    elif np.isfinite(min_q):
        note = f"Minimum adjacent-band BH q = {min_q:.3f}."
    else:
        note = "Corrected adjacent-band tests unavailable."
    ax.text(
        panel_left + 0.13, 23.00, note,
        fontsize=6.6, ha="left", va="bottom", color="#4A4A4A", zorder=27,
    )


def draw_map(ax, data: dict[str, Any], eco_summary: pd.DataFrame):
    eco = data["ecoregions"].copy()
    labels = sorted(eco["ecoregion"].unique())
    lookup = eco_summary.set_index("ecoregion")

    for name in labels:
        polygons = eco[eco["ecoregion"] == name]
        formal = name in lookup.index and bool(lookup.loc[name, "formal_comparison"])
        polygons.plot(
            ax=ax,
            facecolor=ECO_COLORS.get(name, "#E5E5E5") if formal else "#EFEFEF",
            edgecolor="#525252", linewidth=0.75,
            hatch=None if formal else "///", zorder=1,
        )

    for latitude in data["breaks"]["anchor_latitude"]:
        latitude = float(latitude)
        ax.axhline(
            latitude, color=C_BREAK, linestyle=(0, (5, 4)),
            linewidth=0.95, alpha=0.80, zorder=4,
        )
        ax.text(
            -117.55, latitude + 0.035, f"{latitude:.0f}°N",
            fontsize=8.6, fontweight="bold", ha="left", va="bottom",
            color="#444444",
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.9),
            zorder=12,
        )

    formal = eco_summary[eco_summary["formal_comparison"]].copy()
    max_total = float(formal["mean_classified_richness"].max())
    dissolved = eco.dissolve(by="ecoregion")
    points = dissolved.representative_point()

    for row in formal.itertuples(index=False):
        name = str(row.ecoregion)
        if name not in points.index:
            continue
        point = points.loc[name]
        dx, dy = DONUT_OFFSETS.get(name, (0.0, 0.0))
        x, y = float(point.x) + dx, float(point.y) + dy
        radius = draw_donut(
            ax, x, y,
            float(row.mean_classified_richness),
            float(row.mean_ballooning_richness),
            float(row.mean_non_ballooning_richness),
            max_total,
        )
        ratio = float(row.median_B_to_N_ratio)
        label = f"{base.ABBREVIATIONS.get(name, name)}  {ratio:.2f}:1"
        txt = ax.text(
            x, y - radius - 0.075, label,
            ha="center", va="top", fontsize=8.3, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.84),
            zorder=13,
        )
        txt.set_path_effects([
            patheffects.withStroke(linewidth=1.5, foreground="white")
        ])

    draw_symbol_legend(ax)
    draw_integrated_band_ratios(ax, data["band"], data["tests"])

    ax.set_xlim(-117.8, -106.35)
    ax.set_ylim(22.65, 33.15)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Longitude", labelpad=4)
    ax.set_ylabel("Latitude", labelpad=5)
    ax.grid(color="#E6E6E6", linewidth=0.55, zorder=0)


def draw_ecoregion_key(fig):
    ordered = ["CS", "CH", "RCS", "SFD", "CD", "VD", "CGC", "SG", "MP", "TS"]
    reverse = {short: full for full, short in base.ABBREVIATIONS.items()}
    items = []
    for short in ordered:
        full = reverse[short]
        label = base.ENGLISH_LABELS.get(full, full)
        items.append(f"{short} = {label}")

    left = "   •   ".join(items[:5])
    right = "   •   ".join(items[5:])
    fig.text(
        0.05, 0.084, left, ha="left", va="center",
        fontsize=7.4, color="#333333",
    )
    fig.text(
        0.05, 0.058, right, ha="left", va="center",
        fontsize=7.4, color="#333333",
    )


def make_figure(out_dir: Path, data: dict[str, Any], eco_summary: pd.DataFrame, dpi: int):
    publication = out_dir / "publication_outputs"
    publication.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(16.2, 10.8), facecolor="white")
    ax = fig.add_axes([0.045, 0.145, 0.91, 0.76], facecolor="white")

    draw_map(ax, data, eco_summary)

    fig.suptitle(
        "Baja arachnid dispersal-group balance",
        fontsize=19, fontweight="bold", y=0.985, color="#1F1F1F",
    )
    fig.text(
        0.5, 0.951,
        "Equal-cell Ballooning (C3 = D1–D3) versus fixed Non-ballooning (N0); D4 excluded",
        ha="center", fontsize=10.7, color="#333333",
    )
    fig.text(
        0.05, 0.111, "Ecoregion abbreviations",
        ha="left", va="center", fontsize=9.0, fontweight="bold",
    )
    draw_ecoregion_key(fig)
    eco_equal_cells = int(pd.to_numeric(eco_summary["equal_cells"], errors="coerce").dropna().iloc[0])
    eco_iterations = int(data["iterations"]["iteration"].nunique())
    band_equal_cells = int(pd.to_numeric(data["band"]["equal_cells"], errors="coerce").dropna().iloc[0])
    band_iterations = int(pd.to_numeric(data["band"]["iterations"], errors="coerce").dropna().iloc[0])
    fig.text(
        0.95, 0.027,
        f"Map ratios: {eco_iterations:,} draws of {eco_equal_cells} occupied cells per adequately sampled ecoregion. "
        f"Band bars: 2.5th–97.5th percentile intervals from {band_iterations:,} draws of {band_equal_cells} occupied cells.",
        ha="right", va="center", fontsize=7.4, color="#555555",
    )

    png = publication / "10I_final_integrated_map_and_ratios_v7.png"
    pdf = publication / "10I_final_integrated_map_and_ratios_v7.pdf"
    svg = publication / "10I_final_integrated_map_and_ratios_v7.svg"
    jpg = publication / "10I_final_integrated_map_and_ratios_v7.jpg"

    for path in (png, pdf, svg):
        fig.savefig(
            path,
            dpi=dpi if path.suffix == ".png" else None,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="white",
            transparent=False,
        )
    fig.savefig(
        jpg, dpi=dpi, bbox_inches="tight",
        facecolor="white", edgecolor="white", transparent=False,
        pil_kwargs={"quality": 96, "subsampling": 0},
    )
    aliases = {
        png: publication / "Figure_3_Biogeographic_Dispersal_Balance.png",
        pdf: publication / "Figure_3_Biogeographic_Dispersal_Balance.pdf",
        svg: publication / "Figure_3_Biogeographic_Dispersal_Balance.svg",
        jpg: publication / "Figure_3_Biogeographic_Dispersal_Balance.jpg",
    }
    for source, target in aliases.items():
        shutil.copy2(source, target)
    plt.close(fig)
    return png


base.make_figure = make_figure

if __name__ == "__main__":
    raise SystemExit(base.main())
