#!/usr/bin/env python3
"""
Step 11N — Build a map with an inset boundary-effects panel.

Panel A:
    Baja California is the geographic background.
    Each latitude band has one compact results box containing:
      q=0  Richness
      q=1  Common-genera diversity
      q=2  Dominant-genera diversity
    Values are ballooning-capable / non-ballooning.

Panel B:
    An inset effect plot showing ballooning-minus-non-ballooning contrasts in
    the two additive Baselga components of Jaccard dissimilarity:
      Jaccard replacement
      Jaccard nestedness-resultant

Important:
    Panel B has no Hill q values. Hill q-specific boundary quantities are Hill
    beta diversity and are conceptually different from Baselga replacement and
    nestedness. They should be shown separately in the supplement.

Required Step 11K files:
    09_iNEXT_COVERAGE_STANDARDIZED_HILL.csv
    06_TRAIT_CONTRAST_SUMMARY.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd

mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["svg.hashsalt"] = "baja-ballooning-step11N-v8"

BANDS = ["23-24N", "24-26N", "26-28N", "28-30N", "30-32N"]
BAND_CENTERS = {
    "23-24N": 23.48,
    "24-26N": 25.0,
    "26-28N": 27.0,
    "28-30N": 29.0,
    "30-32N": 31.0,
}
BOUNDARIES = [
    ("23-24N", "24-26N"),
    ("24-26N", "26-28N"),
    ("26-28N", "28-30N"),
    ("28-30N", "30-32N"),
]

Q_STYLES = {
    0: ("Richness", "#009E73"),
    1: ("Common genera", "#0072B2"),
    2: ("Dominant genera", "#E69F00"),
}
REPLACEMENT_COLOR = "#D55E00"
NESTEDNESS_COLOR = "#0072B2"
MAP_FILL = "#DCEAF3"
MAP_EDGE = "#4B5563"
LATITUDE_LINE = "#9CA3AF"
PANEL_EDGE = "#6B7280"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build Step 11N map with inset Baselga boundary-effects panel."
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        required=True,
        help="Root of Baja_Ballooning_Pipeline.",
    )
    parser.add_argument(
        "--qc-dir",
        type=Path,
        default=None,
        help="Completed Step 11K directory. Auto-detected when omitted.",
    )
    parser.add_argument(
        "--map-gpkg",
        type=Path,
        default=None,
        help="Validated mainland-only GeoPackage. Auto-detected when omitted.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Defaults to <qc-dir>/figures/main_figure_v8.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=400,
        help="PNG output resolution.",
    )
    parser.add_argument(
        "--no-map",
        action="store_true",
        help="Diagnostic mode without the local GIS layer.",
    )
    return parser.parse_args()


def first_existing(paths: Iterable[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path.resolve()
    return None


def locate_qc_dir(root: Path, override: Path | None) -> Path:
    if override is not None:
        if not override.is_dir():
            raise FileNotFoundError(f"QC directory does not exist: {override}")
        return override.resolve()

    candidates = [
        root / "04_analysis" / "11K_publication_nestedness_replacement_QC",
        root / "04_analysis_USE _THIS" / "11K_publication_nestedness_replacement_QC",
    ]
    found = first_existing(candidates)
    if found is not None:
        return found

    matches = [
        p for p in root.rglob("11K_publication_nestedness_replacement_QC")
        if p.is_dir()
    ]
    if len(matches) == 1:
        return matches[0].resolve()
    if len(matches) > 1:
        raise RuntimeError(
            "Multiple Step 11K directories found. Supply --qc-dir explicitly:\n"
            + "\n".join(str(p) for p in matches)
        )
    raise FileNotFoundError(
        f"Could not locate 11K_publication_nestedness_replacement_QC under {root}"
    )


def locate_map_gpkg(root: Path, override: Path | None) -> Path:
    if override is not None:
        if not override.is_file():
            raise FileNotFoundError(f"Map GeoPackage does not exist: {override}")
        return override.resolve()

    rel = (
        Path("C3_pipeline_rebuild")
        / "09_C3_biogeographic_concordance"
        / "10A_ecoregion_gis_audit"
        / "10A_ecoregions_validated_mainland_only.gpkg"
    )
    candidates = [
        root / "04_analysis_USE _THIS" / rel,
        root / "04_analysis" / rel,
    ]
    found = first_existing(candidates)
    if found is not None:
        return found

    matches = list(root.rglob("10A_ecoregions_validated_mainland_only.gpkg"))
    if len(matches) == 1:
        return matches[0].resolve()
    if len(matches) > 1:
        preferred = [
            p for p in matches if "04_analysis_USE _THIS" in str(p)
        ]
        if len(preferred) == 1:
            return preferred[0].resolve()
        raise RuntimeError(
            "Multiple validated mainland GeoPackages found. "
            "Supply --map-gpkg explicitly:\n"
            + "\n".join(str(p) for p in matches)
        )
    raise FileNotFoundError(
        "Could not locate 10A_ecoregions_validated_mainland_only.gpkg."
    )


def require_file(directory: Path, name: str) -> Path:
    path = directory / name
    if not path.is_file():
        raise FileNotFoundError(f"Required file missing: {path}")
    return path


def load_tables(qc_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    hill = pd.read_csv(
        require_file(qc_dir, "09_iNEXT_COVERAGE_STANDARDIZED_HILL.csv")
    )
    contrast = pd.read_csv(
        require_file(qc_dir, "06_TRAIT_CONTRAST_SUMMARY.csv")
    )

    band_hill = hill[
        hill["scope"].eq("band")
        & hill["assemblage"].isin(BANDS)
        & hill["trait_class"].isin(["ballooning", "non_ballooning"])
        & hill["Order.q"].isin([0, 1, 2])
    ].copy()

    expected = len(BANDS) * 2 * 3
    if len(band_hill) != expected:
        raise ValueError(
            f"Expected {expected} band-level Hill records; found {len(band_hill)}."
        )

    target_coverage = band_hill["target_coverage"].dropna().unique()
    if len(target_coverage) != 1:
        raise ValueError(
            "Band-level Hill estimates do not share one target coverage."
        )

    adjacent = contrast[contrast["adjacent"].eq(True)]
    for band_1, band_2 in BOUNDARIES:
        rows = adjacent[
            adjacent["band_1"].eq(band_1)
            & adjacent["band_2"].eq(band_2)
            & adjacent["metric"].isin(
                ["jaccard_turnover", "jaccard_nestedness"]
            )
        ]
        if set(rows["metric"]) != {
            "jaccard_turnover",
            "jaccard_nestedness",
        }:
            raise ValueError(
                f"Missing Baselga effects for {band_1} versus {band_2}."
            )
        if not (rows["n"] == 5000).all():
            raise ValueError(
                f"Expected 5,000 iterations for {band_1} versus {band_2}."
            )

    return band_hill, contrast


def band_label(band: str) -> str:
    lo, hi = band.replace("N", "").split("-")
    return f"{lo}–{hi}°N"


def get_q(
    hill: pd.DataFrame,
    band: str,
    trait: str,
    order_q: int,
) -> float:
    rows = hill[
        hill["scope"].eq("band")
        & hill["assemblage"].eq(band)
        & hill["trait_class"].eq(trait)
        & hill["Order.q"].eq(order_q)
    ]
    if len(rows) != 1:
        raise ValueError(
            f"Expected one record for {band}/{trait}/q={order_q}; found {len(rows)}."
        )
    return float(rows.iloc[0]["qD"])


def get_effect(
    contrast: pd.DataFrame,
    band_1: str,
    band_2: str,
    metric: str,
) -> pd.Series:
    rows = contrast[
        contrast["adjacent"].eq(True)
        & contrast["band_1"].eq(band_1)
        & contrast["band_2"].eq(band_2)
        & contrast["metric"].eq(metric)
    ]
    if len(rows) != 1:
        raise ValueError(
            f"Expected one effect for {band_1}/{band_2}/{metric}; found {len(rows)}."
        )
    return rows.iloc[0]


def draw_baja(ax: plt.Axes, gpkg: Path) -> None:
    try:
        import geopandas as gpd
    except ImportError as exc:
        raise RuntimeError(
            "geopandas is required for the map. Install repository dependencies "
            "or use --no-map for diagnostics."
        ) from exc

    gdf = gpd.read_file(gpkg)
    if gdf.empty:
        raise ValueError(f"Validated Baja GeoPackage is empty: {gpkg}")
    if gdf.crs is None:
        raise ValueError(f"Validated Baja GeoPackage has no CRS: {gpkg}")

    gdf = gdf.to_crs("EPSG:4326")
    geometry = (
        gdf.geometry.union_all()
        if hasattr(gdf.geometry, "union_all")
        else gdf.geometry.unary_union
    )
    outline = gpd.GeoSeries([geometry], crs="EPSG:4326")
    outline.plot(
        ax=ax,
        facecolor=MAP_FILL,
        edgecolor="none",
        linewidth=0,
        alpha=0.88,
        zorder=0,
    )
    outline.boundary.plot(
        ax=ax,
        color=MAP_EDGE,
        linewidth=1.2,
        zorder=1,
    )


def draw_band_box(
    ax: plt.Axes,
    hill: pd.DataFrame,
    band: str,
    x0: float,
    center_y: float,
) -> None:
    width = 2.55
    height = 0.92
    y0 = center_y - height / 2

    panel = FancyBboxPatch(
        (x0, y0),
        width,
        height,
        boxstyle="round,pad=0.035,rounding_size=0.055",
        facecolor="white",
        edgecolor=PANEL_EDGE,
        linewidth=0.85,
        alpha=0.94,
        zorder=7,
    )
    ax.add_patch(panel)

    ax.text(
        x0 + width / 2,
        y0 + height - 0.12,
        f"{band_label(band)}    B / N",
        ha="center",
        va="center",
        fontsize=9.1,
        fontweight="bold",
        color="#111827",
        zorder=8,
    )

    row_y = [
        y0 + 0.55,
        y0 + 0.34,
        y0 + 0.13,
    ]
    for order_q, y in zip([0, 1, 2], row_y):
        descriptor, color = Q_STYLES[order_q]
        b_value = get_q(hill, band, "ballooning", order_q)
        n_value = get_q(hill, band, "non_ballooning", order_q)

        ax.text(
            x0 + 0.12,
            y,
            f"q={order_q}  {descriptor}",
            ha="left",
            va="center",
            fontsize=7.7,
            fontweight="bold",
            color=color,
            zorder=8,
        )
        ax.text(
            x0 + width - 0.12,
            y,
            f"{b_value:.1f} / {n_value:.1f}",
            ha="right",
            va="center",
            fontsize=8.2,
            fontweight="bold",
            color="#111827",
            zorder=8,
        )


def draw_boundary_inset(
    parent_ax: plt.Axes,
    contrast: pd.DataFrame,
) -> plt.Axes:
    inset = parent_ax.inset_axes([0.725, 0.08, 0.245, 0.82])
    inset.set_facecolor((1, 1, 1, 0.96))

    for spine in inset.spines.values():
        spine.set_color(PANEL_EDGE)
        spine.set_linewidth(0.9)

    inset.axvline(
        0,
        color="#6B7280",
        linewidth=0.9,
        linestyle="--",
        zorder=1,
    )

    y_positions = [1, 2, 3, 4]
    replacement_offset = 0.10
    nestedness_offset = -0.10

    for (band_1, band_2), y in zip(BOUNDARIES, y_positions):
        replacement = get_effect(
            contrast,
            band_1,
            band_2,
            "jaccard_turnover",
        )
        nestedness = get_effect(
            contrast,
            band_1,
            band_2,
            "jaccard_nestedness",
        )

        for row, y_value, color, marker in [
            (
                replacement,
                y + replacement_offset,
                REPLACEMENT_COLOR,
                "o",
            ),
            (
                nestedness,
                y + nestedness_offset,
                NESTEDNESS_COLOR,
                "s",
            ),
        ]:
            median = float(row["median"])
            lower = float(row["p025"])
            upper = float(row["p975"])
            resolved = lower > 0 or upper < 0

            inset.errorbar(
                median,
                y_value,
                xerr=np.array(
                    [[median - lower], [upper - median]]
                ),
                fmt=marker,
                markersize=5.5,
                markerfacecolor=color if marker == "o" else "white",
                markeredgecolor=color,
                markeredgewidth=1.1,
                color=color,
                linewidth=1.35,
                capsize=3,
                zorder=4,
            )
            if resolved:
                inset.text(
                    median,
                    y_value + 0.13,
                    "*",
                    ha="center",
                    va="bottom",
                    fontsize=10,
                    fontweight="bold",
                    color=color,
                    zorder=5,
                )

    inset.set_yticks(y_positions)
    inset.set_yticklabels(
        [
            "23–24 ↔ 24–26",
            "24–26 ↔ 26–28",
            "26–28 ↔ 28–30",
            "28–30 ↔ 30–32",
        ],
        fontsize=7.2,
    )
    inset.set_ylim(0.55, 4.55)
    inset.set_xlim(-0.70, 0.70)
    inset.set_xticks([-0.6, -0.3, 0.0, 0.3, 0.6])
    inset.tick_params(axis="x", labelsize=7)
    inset.grid(
        axis="x",
        color="#E5E7EB",
        linewidth=0.55,
        zorder=0,
    )

    inset.text(
        0.02,
        0.98,
        "B",
        transform=inset.transAxes,
        ha="left",
        va="top",
        fontsize=16,
        fontweight="bold",
    )
    inset.text(
        0.15,
        0.98,
        "Boundary effects",
        transform=inset.transAxes,
        ha="left",
        va="top",
        fontsize=9.3,
        fontweight="bold",
    )
    inset.text(
        0.15,
        0.925,
        "Baselga Δ = B − N\n(no Hill q values)",
        transform=inset.transAxes,
        ha="left",
        va="top",
        fontsize=7.0,
        color="#4B5563",
        linespacing=1.05,
    )

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color=REPLACEMENT_COLOR,
            markerfacecolor=REPLACEMENT_COLOR,
            markersize=5,
            linewidth=1.3,
            label="Replacement",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color=NESTEDNESS_COLOR,
            markerfacecolor="white",
            markeredgecolor=NESTEDNESS_COLOR,
            markersize=5,
            linewidth=1.3,
            label="Nestedness",
        ),
    ]
    inset.legend(
        handles=handles,
        loc="lower right",
        frameon=False,
        fontsize=6.8,
    )
    inset.set_xlabel(
        "Trait contrast",
        fontsize=7.3,
    )
    return inset


def build_figure(
    hill: pd.DataFrame,
    contrast: pd.DataFrame,
    gpkg: Path | None,
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(11.2, 10.2))

    if gpkg is not None:
        draw_baja(ax, gpkg)

    for latitude in [24, 26, 28, 30, 32]:
        ax.axhline(
            latitude,
            color=LATITUDE_LINE,
            linewidth=0.85,
            linestyle="--",
            zorder=2,
        )

    # A consistent left-side column keeps the q values readable and leaves the
    # right side available for the boundary-effects inset.
    box_x = -118.38
    for band in BANDS:
        draw_band_box(
            ax=ax,
            hill=hill,
            band=band,
            x0=box_x,
            center_y=BAND_CENTERS[band],
        )

    draw_boundary_inset(
        parent_ax=ax,
        contrast=contrast,
    )

    ax.text(
        0.015,
        0.985,
        "A",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=18,
        fontweight="bold",
    )

    ax.set_xlim(-118.75, -107.0)
    ax.set_ylim(22.75, 32.60)
    ax.set_xticks([-118, -116, -114, -112, -110, -108])
    ax.set_xticklabels(
        ["118°W", "116°W", "114°W", "112°W", "110°W", "108°W"]
    )
    ax.set_yticks([23, 24, 26, 28, 30, 32])
    ax.set_yticklabels(["23°N", "24°N", "26°N", "28°N", "30°N", "32°N"])
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(
        "Latitudinal diversity and the structure of compositional turnover\n"
        "in Baja California arachnids",
        fontsize=17.5,
        pad=14,
    )

    trait_key = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="#111827",
            markeredgecolor="#111827",
            markersize=5.5,
            label="B = ballooning-capable",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            linestyle="none",
            markerfacecolor="white",
            markeredgecolor="#111827",
            markersize=5.5,
            label="N = non-ballooning",
        ),
    ]
    ax.legend(
        handles=trait_key,
        loc="lower left",
        bbox_to_anchor=(0.01, 0.01),
        frameon=False,
        fontsize=7.8,
    )

    target_coverage = float(hill["target_coverage"].dropna().unique()[0])
    fig.text(
        0.06,
        0.025,
        (
            "Panel A reports coverage-standardized effective genus diversity "
            f"at common coverage C = {target_coverage:.4f}. q=0 is genus richness; "
            "q=1 increasingly weights common genera; q=2 places the greatest "
            "weight on dominant or frequently detected genera. Values are "
            "ballooning-capable / non-ballooning. Panel B is not q-specific: "
            "it reports ballooning-minus-non-ballooning contrasts in the additive "
            "Baselga replacement and nestedness-resultant components of Jaccard "
            "dissimilarity. Points are medians and bars are 2.5th–97.5th "
            "percentile intervals from 5,000 paired equal-cell iterations; "
            "* marks intervals excluding zero."
        ),
        ha="left",
        va="bottom",
        fontsize=8.5,
        wrap=True,
    )

    fig.tight_layout(rect=[0.035, 0.075, 0.995, 0.98])
    return fig


def write_caption(output_dir: Path) -> Path:
    caption = """# Figure caption draft

**Figure 1. Latitudinal diversity and the structure of compositional turnover
in Baja California arachnids.** **(A)** Coverage-standardized effective genus
diversity for ballooning-capable and non-ballooning assemblages within five
latitude bands. q = 0 represents genus richness, q = 1 gives greater weight to
common genera, and q = 2 gives the greatest weight to dominant or frequently
detected genera. Values are ballooning-capable / non-ballooning at a common
sample-coverage threshold. **(B)** Ballooning-capable minus non-ballooning
contrasts in the Jaccard replacement and nestedness-resultant components across
adjacent latitude-band boundaries. Panel B is not q-specific: these are
Baselga components of pairwise compositional dissimilarity rather than Hill
diversity orders. Points are medians and horizontal intervals are the
2.5th–97.5th percentiles from 5,000 paired equal-cell iterations. Asterisks
identify intervals that exclude zero.

Hill beta diversity at q = 0, 1, and 2, observed coverage, total Jaccard
dissimilarity, Simpson replacement, gamma diversity, and the full iNEXT curves
are reported in the supplementary material.
"""
    path = output_dir / "FIGURE_1_CAPTION_DRAFT.md"
    path.write_text(caption, encoding="utf-8")
    return path


def main() -> int:
    args = parse_args()
    root = args.project_root.expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Project root does not exist: {root}")

    qc_dir = locate_qc_dir(root, args.qc_dir)
    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else qc_dir / "figures" / "main_figure_v8"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    gpkg = None if args.no_map else locate_map_gpkg(root, args.map_gpkg)
    hill, contrast = load_tables(qc_dir)

    fig = build_figure(
        hill=hill,
        contrast=contrast,
        gpkg=gpkg,
    )

    stem = "Figure_1_map_with_boundary_effect_inset_v8"
    outputs = [
        output_dir / f"{stem}.png",
        output_dir / f"{stem}.pdf",
        output_dir / f"{stem}.svg",
    ]
    fig.savefig(outputs[0], dpi=args.dpi, bbox_inches="tight")
    fig.savefig(outputs[1], bbox_inches="tight")
    fig.savefig(outputs[2], bbox_inches="tight")
    plt.close(fig)

    caption = write_caption(output_dir)

    print("STEP 11N COMPLETE")
    print(f"PROJECT_ROOT={root}")
    print(f"QC_DIR={qc_dir}")
    print(f"MAP_GPKG={gpkg if gpkg is not None else 'not drawn'}")
    for output in outputs:
        print(f"OUTPUT={output}")
    print(f"OUTPUT={caption}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"STEP 11N FAILED: {exc}", file=sys.stderr)
        raise
