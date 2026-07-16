#!/usr/bin/env python3
"""
STEP 11F — Equal-cell incidence-based rarefaction
Baja Ballooning Publication

Purpose
-------
Create directly comparable rarefaction curves for all five Baja latitude bands
using the same spatial standardization as Step 11B:

  * occupied, approximately equal-area 25-km grid cells are sampling units;
  * each band is restricted to the common maximum of 22 cells (or the minimum
    available cell count if the frozen dataset changes);
  * for each iteration, cells are sampled without replacement and randomly
    ordered;
  * cumulative genus richness is calculated from 1 through 22 cells;
  * the procedure is repeated 2,000 times;
  * the mean expected richness and 95% equal-cell resampling envelope are plotted.

Analyses
--------
  1. primary
  2. taxonomy_strict
  3. explicit_low_confidence_exclusion

The same cell draw and ordering are used for all three analyses within each
band and iteration, making sensitivity comparisons paired.

No third-party packages are required for the analysis. Matplotlib is optional:
when available, PNG, PDF, and SVG figures are written; otherwise a standalone
SVG is still created.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import random
import shutil
import statistics
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence


SCRIPT_VERSION = "11F_mean_expected_publication_v2_2026-07-16"
DEFAULT_ITERATIONS = 2000
DEFAULT_SEED = 20260713

BAND_ORDER = [
    "23-24N",
    "24-26N",
    "26-28N",
    "28-30N",
    "30-32N",
]

BAND_LABELS = {
    "23-24N": "23–24°N",
    "24-26N": "24–26°N",
    "26-28N": "26–28°N",
    "28-30N": "28–30°N",
    "30-32N": "30–32°N",
}

# Matches the established Step 11 visual identity.
BAND_COLORS = {
    "23-24N": "#5B2A86",
    "24-26N": "#B23A8A",
    "26-28N": "#E07A36",
    "28-30N": "#9A8F2B",
    "30-32N": "#169C97",
}

ANALYSIS_ORDER = [
    "primary",
    "taxonomy_strict",
    "explicit_low_confidence_exclusion",
]

ANALYSIS_LABELS = {
    "primary": "Primary biodiversity dataset",
    "taxonomy_strict": "Taxonomy-strict sensitivity analysis",
    "explicit_low_confidence_exclusion": (
        "Sensitivity analysis excluding explicitly LOW-confidence genera"
    ),
}


class Tee:
    def __init__(self, *streams: Any) -> None:
        self.streams = streams

    def write(self, text: str) -> int:
        for stream in self.streams:
            stream.write(text)
            stream.flush()
        return len(text)

    def flush(self) -> None:
        for stream in self.streams:
            stream.flush()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def first_existing(
    paths: Sequence[Path], label: str, required: bool = True
) -> Path | None:
    for path in paths:
        expanded = path.expanduser()
        if expanded.is_file():
            return expanded.resolve()
    if required:
        attempted = "\n".join(str(path.expanduser()) for path in paths)
        raise FileNotFoundError(f"Could not find {label}. Tried:\n{attempted}")
    return None


def find_field(
    fields: Sequence[str],
    candidates: Sequence[str],
    label: str,
    required: bool = True,
) -> str | None:
    lower_to_original = {field.strip().lower(): field for field in fields}
    for candidate in candidates:
        if candidate.lower() in lower_to_original:
            return lower_to_original[candidate.lower()]
    if required:
        raise ValueError(
            f"Could not identify {label}. Tried: {', '.join(candidates)}. "
            f"Available fields: {', '.join(fields)}"
        )
    return None


def normalize_band(value: str) -> str:
    text = value.strip().upper()
    text = text.replace("–", "-").replace("—", "-").replace("°", "")
    text = "".join(text.split())
    if text.endswith("N"):
        text = text[:-1]
    return f"{text}N"


def read_csv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {path}")
        rows = [dict(row) for row in reader]
        return list(reader.fieldnames), rows


def write_csv(
    path: Path, rows: Iterable[dict[str, Any]], fieldnames: Sequence[str]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def read_incidence_matrix(path: Path) -> dict[str, Any]:
    """Read a genus × grid-cell 0/1 matrix and convert cell columns to bitmasks."""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Empty incidence matrix: {path}") from exc

        if len(header) < 2:
            raise ValueError(f"Incidence matrix has fewer than two columns: {path}")

        genus_field = header[0]
        cells = [cell.strip() for cell in header[1:]]
        if any(not cell for cell in cells):
            raise ValueError(f"Blank grid-cell column in {path}")
        if len(set(cells)) != len(cells):
            raise ValueError(f"Duplicate grid-cell columns in {path}")

        genera: list[str] = []
        rows_by_genus: dict[str, list[int]] = {}

        for row_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(
                    f"Row {row_number} in {path} has {len(row)} fields; "
                    f"expected {len(header)}"
                )

            genus = row[0].strip()
            if not genus:
                raise ValueError(f"Blank genus at row {row_number} in {path}")
            key = genus.casefold()
            if key in rows_by_genus:
                raise ValueError(f"Duplicate genus {genus!r} in {path}")

            try:
                values = [int(value) for value in row[1:]]
            except ValueError as exc:
                raise ValueError(
                    f"Non-integer matrix value at row {row_number} in {path}"
                ) from exc

            if any(value not in (0, 1) for value in values):
                raise ValueError(
                    f"Matrix values other than 0 or 1 at row {row_number} in {path}"
                )

            genera.append(genus)
            rows_by_genus[key] = values

    return {
        "path": path,
        "genus_field": genus_field,
        "genera": genera,
        "cells": cells,
        "rows_by_genus": rows_by_genus,
    }


def align_matrix(
    matrix: dict[str, Any],
    genera: Sequence[str],
    cells: Sequence[str],
) -> list[int]:
    """Align a matrix to the primary genus/cell universe and return cell masks."""
    matrix_cells = matrix["cells"]
    cell_position = {cell: index for index, cell in enumerate(matrix_cells)}

    missing_cells = [cell for cell in cells if cell not in cell_position]
    if missing_cells:
        raise ValueError(
            f"Matrix {matrix['path']} lacks required cells: {missing_cells[:10]}"
        )

    cell_masks = [0] * len(cells)

    for genus_index, genus in enumerate(genera):
        key = genus.casefold()
        if key not in matrix["rows_by_genus"]:
            raise ValueError(f"Matrix {matrix['path']} lacks genus {genus!r}")

        source_values = matrix["rows_by_genus"][key]
        bit = 1 << genus_index

        for target_cell_index, cell in enumerate(cells):
            if source_values[cell_position[cell]] == 1:
                cell_masks[target_cell_index] |= bit

    return cell_masks


def load_cell_bands(path: Path, required_cells: Sequence[str]) -> dict[str, str]:
    fields, rows = read_csv_rows(path)
    cell_field = find_field(
        fields,
        ["grid_cell_id", "cell_id", "grid25km_cell_id"],
        "grid-cell ID field",
    )
    band_field = find_field(
        fields,
        ["centroid_latitude_band", "latitude_band", "lat_band"],
        "latitude-band field",
    )

    mapping: dict[str, str] = {}
    for row in rows:
        cell = row.get(cell_field, "").strip()
        band = normalize_band(row.get(band_field, ""))
        if cell:
            mapping[cell] = band

    missing = [cell for cell in required_cells if cell not in mapping]
    if missing:
        raise ValueError(
            f"Cell lookup lacks {len(missing)} required cells. "
            f"Examples: {missing[:10]}"
        )

    invalid = sorted(
        {mapping[cell] for cell in required_cells if mapping[cell] not in BAND_ORDER}
    )
    if invalid:
        raise ValueError(f"Cells assigned outside the five latitude bands: {invalid}")

    return {cell: mapping[cell] for cell in required_cells}


def load_low_confidence_mask(path: Path, genera: Sequence[str]) -> dict[str, Any]:
    """Return a bitmask retaining all genera except those explicitly coded LOW."""
    fields, rows = read_csv_rows(path)
    genus_field = find_field(
        fields,
        ["genus", "analysis_genus"],
        "trait genus field",
    )
    confidence_field = find_field(
        fields,
        [
            "final_confidence",
            "trait_final_confidence",
            "trait_confidence",
            "trait_ballooning_confidence",
            "current_trait_confidence",
        ],
        "trait confidence field",
        required=False,
    )

    # Confidence is optional in the authoritative publication trait table.
    # When absent, retain every genus in the low-confidence sensitivity and
    # label confidence as UNSPECIFIED. This makes the sensitivity identical
    # to the primary analysis rather than aborting a total-richness analysis.
    lookup: dict[str, str] = {}
    for row in rows:
        genus = row.get(genus_field, "").strip()
        if not genus:
            continue
        key = genus.casefold()
        if key in lookup:
            raise ValueError(f"Duplicate genus {genus!r} in trait table {path}")
        confidence = (
            row.get(confidence_field, "").strip().upper()
            if confidence_field is not None
            else "UNSPECIFIED"
        )
        lookup[key] = confidence or "UNSPECIFIED"

    missing = [genus for genus in genera if genus.casefold() not in lookup]
    if missing:
        raise ValueError(
            "Genera missing from trait table: " + ", ".join(missing[:20])
        )

    all_mask = (1 << len(genera)) - 1
    low_mask = 0
    confidence_counts: dict[str, int] = defaultdict(int)

    for genus_index, genus in enumerate(genera):
        confidence = lookup[genus.casefold()]
        confidence_counts[confidence] += 1
        if confidence == "LOW":
            low_mask |= 1 << genus_index

    return {
        "keep_mask": all_mask & ~low_mask,
        "low_mask": low_mask,
        "confidence_counts": dict(sorted(confidence_counts.items())),
        "confidence_field": confidence_field,
    }


def percentile(values: Sequence[int | float], probability: float) -> float:
    ordered = sorted(float(value) for value in values)
    if not ordered:
        return math.nan
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def summarize_values(values: Sequence[int | float]) -> dict[str, float | int]:
    clean = [float(value) for value in values if math.isfinite(float(value))]
    if not clean:
        return {
            "n_iterations": 0,
            "mean": math.nan,
            "median": math.nan,
            "standard_deviation": math.nan,
            "p025": math.nan,
            "p975": math.nan,
            "minimum": math.nan,
            "maximum": math.nan,
        }
    return {
        "n_iterations": len(clean),
        "mean": statistics.fmean(clean),
        "median": statistics.median(clean),
        "standard_deviation": statistics.stdev(clean) if len(clean) > 1 else 0.0,
        "p025": percentile(clean, 0.025),
        "p975": percentile(clean, 0.975),
        "minimum": min(clean),
        "maximum": max(clean),
    }


def archive_existing(output_dir: Path, archive_root: Path) -> Path | None:
    if not output_dir.exists() or not any(output_dir.iterdir()):
        return None
    timestamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    archive_dir = archive_root / f"11F_equal_cell_rarefaction_{timestamp}"
    archive_dir.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(output_dir, archive_dir)
    shutil.rmtree(output_dir)
    return archive_dir


def svg_escape(text: Any) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def write_fallback_svg(
    path: Path,
    summary_rows: Sequence[dict[str, Any]],
    analysis: str,
    sample_size: int,
    iterations: int,
) -> None:
    """Write a standalone SVG without requiring Matplotlib."""
    rows = [row for row in summary_rows if row["analysis"] == analysis]
    width, height = 1100, 760
    left, right, top, bottom = 105, 40, 100, 105
    plot_width = width - left - right
    plot_height = height - top - bottom

    y_max = max(float(row["p975"]) for row in rows)
    y_max = max(10.0, math.ceil(y_max / 10.0) * 10.0)

    def x_pos(k: int) -> float:
        return left + plot_width * (k - 1) / max(1, sample_size - 1)

    def y_pos(value: float) -> float:
        return top + plot_height * (y_max - value) / y_max

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#111}.grid{stroke:#dddddd;stroke-width:1}.axis{stroke:#111;stroke-width:1.5}</style>',
        f'<text x="{width/2}" y="38" text-anchor="middle" font-size="25" font-weight="bold">Equal-cell incidence-based rarefaction of arachnid genera</text>',
        f'<text x="{width/2}" y="70" text-anchor="middle" font-size="16">{svg_escape(ANALYSIS_LABELS[analysis])}; {iterations:,} resamples; 1–{sample_size} occupied 25-km cells per band</text>',
    ]

    for tick in range(0, 6):
        value = y_max * tick / 5
        y = y_pos(value)
        parts.append(f'<line class="grid" x1="{left}" x2="{width-right}" y1="{y}" y2="{y}"/>')
        parts.append(f'<text x="{left-12}" y="{y+5}" text-anchor="end" font-size="15">{value:.0f}</text>')

    for k in sorted(set([1, 5, 10, 15, 20, sample_size])):
        if 1 <= k <= sample_size:
            x = x_pos(k)
            parts.append(f'<line class="grid" x1="{x}" x2="{x}" y1="{top}" y2="{height-bottom}"/>')
            parts.append(f'<text x="{x}" y="{height-bottom+30}" text-anchor="middle" font-size="15">{k}</text>')

    parts.extend([
        f'<line class="axis" x1="{left}" x2="{left}" y1="{top}" y2="{height-bottom}"/>',
        f'<line class="axis" x1="{left}" x2="{width-right}" y1="{height-bottom}" y2="{height-bottom}"/>',
        f'<text x="{width/2}" y="{height-35}" text-anchor="middle" font-size="18">Number of occupied 25-km grid cells</text>',
        f'<text transform="translate(28,{top+plot_height/2}) rotate(-90)" text-anchor="middle" font-size="18">Mean expected genus richness</text>',
    ])

    for band_index, band in enumerate(BAND_ORDER):
        band_rows = sorted(
            [row for row in rows if row["latitude_band"] == band],
            key=lambda row: int(row["n_cells"]),
        )
        upper_points = " ".join(
            f"{x_pos(int(row['n_cells'])):.2f},{y_pos(float(row['p975'])):.2f}"
            for row in band_rows
        )
        lower_points = " ".join(
            f"{x_pos(int(row['n_cells'])):.2f},{y_pos(float(row['p025'])):.2f}"
            for row in reversed(band_rows)
        )
        mean_points = " ".join(
            f"{x_pos(int(row['n_cells'])):.2f},{y_pos(float(row['mean'])):.2f}"
            for row in band_rows
        )
        color = BAND_COLORS[band]
        parts.append(
            f'<polygon points="{upper_points} {lower_points}" fill="{color}" fill-opacity="0.16" stroke="none"/>'
        )
        parts.append(
            f'<polyline points="{mean_points}" fill="none" stroke="{color}" stroke-width="3.5"/>'
        )

        legend_y = top + 24 + band_index * 28
        legend_x = width - right - 220
        parts.append(f'<line x1="{legend_x}" x2="{legend_x+34}" y1="{legend_y}" y2="{legend_y}" stroke="{color}" stroke-width="4"/>')
        parts.append(f'<text x="{legend_x+44}" y="{legend_y+5}" font-size="15">{svg_escape(BAND_LABELS[band])}</text>')

    parts.append(f'<text x="{left}" y="{height-12}" font-size="13">Lines show mean expected richness; ribbons show 2.5th–97.5th percentile equal-cell resampling envelopes.</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def plot_with_matplotlib(
    figure_dir: Path,
    summary_rows: Sequence[dict[str, Any]],
    analysis: str,
    sample_size: int,
    iterations: int,
) -> list[dict[str, Any]]:
    """Create SVG, PNG, and PDF figures when Matplotlib is available."""
    status: list[dict[str, Any]] = []
    basename = f"Figure_1_equal_cell_rarefaction_{analysis}"

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        svg_path = figure_dir / f"{basename}.svg"
        write_fallback_svg(
            svg_path, summary_rows, analysis, sample_size, iterations
        )
        status.append(
            {
                "analysis": analysis,
                "format": "svg",
                "success": True,
                "path": str(svg_path),
                "note": f"Fallback SVG used because Matplotlib was unavailable: {exc}",
            }
        )
        return status

    rows = [row for row in summary_rows if row["analysis"] == analysis]
    fig, ax = plt.subplots(figsize=(10.5, 7.1))

    for band in BAND_ORDER:
        band_rows = sorted(
            [row for row in rows if row["latitude_band"] == band],
            key=lambda row: int(row["n_cells"]),
        )
        x = [int(row["n_cells"]) for row in band_rows]
        mean_expected = [float(row["mean"]) for row in band_rows]
        lower = [float(row["p025"]) for row in band_rows]
        upper = [float(row["p975"]) for row in band_rows]
        color = BAND_COLORS[band]
        ax.fill_between(x, lower, upper, color=color, alpha=0.16, linewidth=0)
        ax.plot(x, mean_expected, color=color, linewidth=2.6, label=BAND_LABELS[band])

    ax.set_title(
        "Equal-cell incidence-based rarefaction of arachnid genera",
        fontsize=16,
        pad=18,
    )
    ax.text(
        0.5,
        1.015,
        (
            f"{ANALYSIS_LABELS[analysis]}; {iterations:,} resamples; "
            f"1–{sample_size} occupied 25-km cells per band"
        ),
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=10.5,
    )
    ax.set_xlabel("Number of occupied 25-km grid cells", fontsize=12)
    ax.set_ylabel("Mean expected genus richness", fontsize=12)
    ax.set_xlim(1, sample_size)
    ax.set_xticks(sorted(set([1, 5, 10, 15, 20, sample_size])))
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.25, linewidth=0.8)
    ax.legend(title="Latitude band", frameon=False, loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.text(
        0.01,
        0.01,
        "Lines show mean expected richness; ribbons show 2.5th–97.5th percentile equal-cell resampling envelopes.",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.035, 1, 0.96))

    for extension, kwargs in [
        ("png", {"dpi": 400}),
        ("pdf", {}),
        ("svg", {}),
    ]:
        path = figure_dir / f"{basename}.{extension}"
        try:
            fig.savefig(path, bbox_inches="tight", **kwargs)
            status.append(
                {
                    "analysis": analysis,
                    "format": extension,
                    "success": True,
                    "path": str(path),
                    "note": "",
                }
            )
        except Exception as exc:
            status.append(
                {
                    "analysis": analysis,
                    "format": extension,
                    "success": False,
                    "path": str(path),
                    "note": str(exc),
                }
            )

    plt.close(fig)
    return status


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Equal-cell incidence-based rarefaction for Baja arachnid genera."
        )
    )
    parser.add_argument(
        "project_root",
        nargs="?",
        default="~/Desktop/Baja_Ballooning_Pipeline",
        help="Baja Ballooning Pipeline project root.",
    )
    parser.add_argument(
        "iterations",
        nargs="?",
        type=int,
        default=DEFAULT_ITERATIONS,
        help=f"Number of Monte Carlo iterations; default {DEFAULT_ITERATIONS}.",
    )
    parser.add_argument(
        "seed",
        nargs="?",
        type=int,
        default=DEFAULT_SEED,
        help=f"Random-number seed; default {DEFAULT_SEED}.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.iterations < 100:
        raise ValueError("Use at least 100 Monte Carlo iterations.")

    project_root = Path(args.project_root).expanduser().resolve()
    if not project_root.is_dir():
        raise FileNotFoundError(f"Project root not found: {project_root}")

    analysis_ready = project_root / "ANALYSIS_READY_INPUTS"
    grid_fallback = project_root / "02_data_clean" / "08_grid25km_incidence"
    step11_dir = project_root / "04_analysis" / "11_latitude_band_diversity_turnover"
    trait_fallback = project_root / "02_data_clean" / "07_final_trait_merge"

    output_dir = project_root / "04_analysis" / "11F_equal_cell_rarefaction"
    figure_dir = output_dir / "figures"
    archive_root = project_root / "08_archive"

    archived = archive_existing(output_dir, archive_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    log_path = output_dir / "11F_analysis_log.txt"
    with log_path.open("w", encoding="utf-8") as log_handle:
        original_stdout = sys.stdout
        original_stderr = sys.stderr
        sys.stdout = Tee(original_stdout, log_handle)
        sys.stderr = Tee(original_stderr, log_handle)

        try:
            print("STEP 11F STARTED")
            print(f"Version: {SCRIPT_VERSION}")
            print(f"Project: {project_root}")
            print(f"Iterations: {args.iterations}")
            print(f"Random seed: {args.seed}")
            if archived is not None:
                print(f"Archived prior Step 11F output: {archived}")

            primary_path = first_existing(
                [
                    analysis_ready
                    / "02_incidence_matrices_25km"
                    / "10_biodiversity_final_genus_by_grid25km_incidence.csv",
                    grid_fallback
                    / "10_biodiversity_final_genus_by_grid25km_incidence.csv",
                ],
                "primary biodiversity incidence matrix",
            )
            strict_path = first_existing(
                [
                    analysis_ready
                    / "02_incidence_matrices_25km"
                    / "10_biodiversity_taxonomy_strict_genus_by_grid25km_incidence.csv",
                    grid_fallback
                    / "10_biodiversity_taxonomy_strict_genus_by_grid25km_incidence.csv",
                ],
                "taxonomy-strict incidence matrix",
            )
            cell_lookup_path = first_existing(
                [
                    analysis_ready
                    / "04_spatial_reference"
                    / "10_common_grid25km_cell_lookup.csv",
                    step11_dir / "11_grid25km_cell_to_latitude_band.csv",
                    grid_fallback / "10_common_grid25km_cell_lookup.csv",
                ],
                "25-km cell lookup",
            )
            trait_path = first_existing(
                [
                    analysis_ready
                    / "03_trait_tables"
                    / "07_reviewed_genus_trait_lookup_final.csv",
                    step11_dir
                    / "11_biodiversity_final_genus_by_latitude_band_incidence.csv",
                    trait_fallback / "07_reviewed_genus_trait_lookup_final.csv",
                ],
                "final genus trait table",
            )

            print(f"Primary matrix: {primary_path}")
            print(f"Taxonomy-strict matrix: {strict_path}")
            print(f"Cell lookup: {cell_lookup_path}")
            print(f"Trait table: {trait_path}")

            primary_raw = read_incidence_matrix(primary_path)
            strict_raw = read_incidence_matrix(strict_path)
            genera = primary_raw["genera"]
            cells = primary_raw["cells"]

            if any(genus.casefold() == "fesa" for genus in genera):
                raise ValueError("Fesa remains in the primary genus universe.")
            if any(genus.casefold() == "rhipicephalus" for genus in genera):
                print(
                    "NOTICE: genus Rhipicephalus occurs in the incidence matrix. "
                    "The frozen upstream dataset should already exclude "
                    "Rhipicephalus sanguineus from primary biodiversity analyses."
                )

            primary_masks = align_matrix(primary_raw, genera, cells)
            strict_masks = align_matrix(strict_raw, genera, cells)
            cell_bands = load_cell_bands(cell_lookup_path, cells)
            trait_info = load_low_confidence_mask(trait_path, genera)

            band_indices: dict[str, list[int]] = {
                band: [
                    index
                    for index, cell in enumerate(cells)
                    if cell_bands[cell] == band
                ]
                for band in BAND_ORDER
            }
            available_cells = {
                band: len(indices) for band, indices in band_indices.items()
            }
            if any(count == 0 for count in available_cells.values()):
                raise ValueError(
                    f"At least one latitude band has no occupied cells: {available_cells}"
                )

            sample_size = min(available_cells.values())
            if sample_size != 22:
                print(
                    f"NOTICE: common equal-cell limit is {sample_size}, not the "
                    "previously documented 22. The script will use the minimum "
                    "available count from the current frozen inputs."
                )

            print(f"Genera in primary universe: {len(genera)}")
            print(f"Occupied grid cells: {len(cells)}")
            print(f"Available cells by band: {available_cells}")
            print(f"Equal-cell rarefaction limit: {sample_size}")
            print(
                "Explicitly LOW-confidence genera excluded in sensitivity: "
                f"{trait_info['low_mask'].bit_count()}"
            )

            values: dict[tuple[str, str, int], list[int]] = defaultdict(list)
            raw_path = output_dir / "11F_equal_cell_rarefaction_iterations.csv.gz"
            rng = random.Random(args.seed)

            with gzip.open(raw_path, "wt", encoding="utf-8", newline="") as raw_handle:
                raw_fields = [
                    "iteration",
                    "analysis",
                    "latitude_band",
                    "n_cells",
                    "genus_richness",
                ]
                writer = csv.DictWriter(raw_handle, fieldnames=raw_fields)
                writer.writeheader()

                for iteration in range(1, args.iterations + 1):
                    for band in BAND_ORDER:
                        ordered_indices = rng.sample(
                            band_indices[band], sample_size
                        )

                        cumulative_primary = 0
                        cumulative_strict = 0
                        for k, cell_index in enumerate(ordered_indices, start=1):
                            cumulative_primary |= primary_masks[cell_index]
                            cumulative_strict |= strict_masks[cell_index]

                            richness_by_analysis = {
                                "primary": cumulative_primary.bit_count(),
                                "taxonomy_strict": cumulative_strict.bit_count(),
                                "explicit_low_confidence_exclusion": (
                                    cumulative_primary
                                    & trait_info["keep_mask"]
                                ).bit_count(),
                            }

                            for analysis in ANALYSIS_ORDER:
                                richness = richness_by_analysis[analysis]
                                values[(analysis, band, k)].append(richness)
                                writer.writerow(
                                    {
                                        "iteration": iteration,
                                        "analysis": analysis,
                                        "latitude_band": band,
                                        "n_cells": k,
                                        "genus_richness": richness,
                                    }
                                )

                    if iteration % max(1, args.iterations // 10) == 0:
                        print(
                            f"Completed {iteration:,} / {args.iterations:,} iterations"
                        )

            summary_rows: list[dict[str, Any]] = []
            for analysis in ANALYSIS_ORDER:
                for band in BAND_ORDER:
                    for k in range(1, sample_size + 1):
                        summary = summarize_values(values[(analysis, band, k)])
                        summary_rows.append(
                            {
                                "analysis": analysis,
                                "latitude_band": band,
                                "latitude_band_label": BAND_LABELS[band],
                                "n_cells": k,
                                **summary,
                            }
                        )

            summary_path = output_dir / "11F_equal_cell_rarefaction_summary.csv"
            summary_fields = [
                "analysis",
                "latitude_band",
                "latitude_band_label",
                "n_cells",
                "n_iterations",
                "mean",
                "median",
                "standard_deviation",
                "p025",
                "p975",
                "minimum",
                "maximum",
            ]
            write_csv(summary_path, summary_rows, summary_fields)

            endpoint_rows = [
                row for row in summary_rows if int(row["n_cells"]) == sample_size
            ]
            endpoint_path = output_dir / "11F_equal_cell_richness_at_common_limit.csv"
            write_csv(endpoint_path, endpoint_rows, summary_fields)

            design_rows = [
                {
                    "latitude_band": band,
                    "latitude_band_label": BAND_LABELS[band],
                    "available_occupied_cells": available_cells[band],
                    "common_rarefaction_limit_cells": sample_size,
                    "nominal_cell_area_km2": 625,
                    "nominal_sampled_grid_area_at_limit_km2": sample_size * 625,
                }
                for band in BAND_ORDER
            ]
            design_path = output_dir / "11F_sampling_design_by_band.csv"
            write_csv(
                design_path,
                design_rows,
                [
                    "latitude_band",
                    "latitude_band_label",
                    "available_occupied_cells",
                    "common_rarefaction_limit_cells",
                    "nominal_cell_area_km2",
                    "nominal_sampled_grid_area_at_limit_km2",
                ],
            )

            export_status: list[dict[str, Any]] = []
            for analysis in ANALYSIS_ORDER:
                export_status.extend(
                    plot_with_matplotlib(
                        figure_dir,
                        summary_rows,
                        analysis,
                        sample_size,
                        args.iterations,
                    )
                )

            export_status_path = output_dir / "11F_figure_export_status.csv"
            write_csv(
                export_status_path,
                export_status,
                ["analysis", "format", "success", "path", "note"],
            )

            caption = (
                "Figure 1. Equal-cell incidence-based rarefaction of arachnid "
                "genus richness across five latitude bands of the Baja California "
                "Peninsula. Occupied, approximately equal-area 25-km grid cells "
                "were treated as incidence sampling units. For each latitude band, "
                f"cells were randomly sampled without replacement and ordered from "
                f"1 to {sample_size} cells across {args.iterations:,} Monte Carlo "
                "iterations. Lines show mean expected accumulated genus richness and "
                "shaded ribbons show the 2.5th–97.5th percentile subset-resampling "
                "envelope. The common cell limit standardizes the number and nominal "
                "area of occupied sampling units among bands; it does not equalize "
                "within-cell collection intensity, detectability, temporal coverage, "
                "or the total land area of each latitude band."
            )
            (output_dir / "Figure_1_caption.txt").write_text(
                caption + "\n", encoding="utf-8"
            )

            input_paths = [
                primary_path,
                strict_path,
                cell_lookup_path,
                trait_path,
            ]
            input_manifest = [
                {
                    "path": str(path),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for path in input_paths
            ]
            write_csv(
                output_dir / "11F_input_manifest.csv",
                input_manifest,
                ["path", "bytes", "sha256"],
            )

            provenance = {
                "step": "11F_equal_cell_rarefaction",
                "script_version": SCRIPT_VERSION,
                "created_utc": utc_now(),
                "project_root": str(project_root),
                "iterations": args.iterations,
                "random_seed": args.seed,
                "sampling_design": {
                    "sampling_unit": "occupied approximately equal-area 25-km grid cell",
                    "without_replacement": True,
                    "common_limit_cells": sample_size,
                    "available_cells_by_band": available_cells,
                    "nominal_cell_area_km2": 625,
                    "nominal_sampled_grid_area_at_limit_km2": sample_size * 625,
                    "same_cell_order_used_across_sensitivity_analyses": True,
                },
                "analyses": {
                    "primary": "Primary biodiversity genus-by-cell incidence matrix",
                    "taxonomy_strict": "Taxonomy-strict incidence matrix",
                    "explicit_low_confidence_exclusion": (
                        "Primary incidence matrix after excluding genera explicitly "
                        "coded LOW confidence; UNSPECIFIED legacy assignments retained"
                    ),
                },
                "interval_definition": (
                    "2.5th and 97.5th percentiles among equal-cell subsets of the "
                    "observed occupied-cell universe"
                ),
                "input_manifest": input_manifest,
                "trait_confidence_counts": trait_info["confidence_counts"],
                "validation_passed": True,
            }
            (output_dir / "11F_provenance.json").write_text(
                json.dumps(provenance, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )

            running_script = Path(__file__).resolve()
            shutil.copy2(
                running_script,
                output_dir / "11F_equal_cell_rarefaction.py",
            )

            output_files = sorted(
                path
                for path in output_dir.rglob("*")
                if path.is_file()
                and path.name != "11F_output_file_manifest.csv"
            )
            output_manifest = [
                {
                    "relative_path": str(path.relative_to(output_dir)),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for path in output_files
            ]
            write_csv(
                output_dir / "11F_output_file_manifest.csv",
                output_manifest,
                ["relative_path", "bytes", "sha256"],
            )

            failed_exports = [
                row for row in export_status if not bool(row["success"])
            ]
            if failed_exports:
                print(
                    "NOTICE: one or more optional figure exports failed. "
                    "See 11F_figure_export_status.csv."
                )

            print()
            print("=" * 78)
            print("STEP 11F COMPLETED SUCCESSFULLY")
            print("=" * 78)
            print(f"Iterations: {args.iterations:,}")
            print(f"Equal cells per latitude band: {sample_size}")
            print(
                "Primary manuscript figure: "
                f"{figure_dir / 'Figure_1_equal_cell_rarefaction_primary.png'}"
            )
            print(f"Outputs: {output_dir}")
            return 0

        finally:
            sys.stdout = original_stdout
            sys.stderr = original_stderr


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nStep 11F interrupted by user.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:
        print(f"\nSTEP 11F FAILED: {exc}", file=sys.stderr)
        raise
