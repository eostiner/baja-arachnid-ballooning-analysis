#!/usr/bin/env python3
"""
STEP 11B — Equal-cell Monte Carlo resampling
Baja Ballooning Publication

Draws the same number of occupied 25-km cells without replacement from each
latitude band, then recalculates richness, ballooning composition, and
pairwise beta diversity.

Core analyses:
  1. primary biodiversity matrix
  2. taxonomy-strict sensitivity matrix
  3. explicit LOW-confidence trait exclusion sensitivity

No third-party Python packages are required.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import re
import shutil
import statistics
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence


SCRIPT_VERSION = "11B_C3_publication_v2_2026-07-16"
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

ANALYSIS_ORDER = [
    "primary",
    "taxonomy_strict",
    "explicit_low_confidence_exclusion",
]

PAIR_ORDER = [
    ("23-24N", "24-26N"),
    ("23-24N", "26-28N"),
    ("23-24N", "28-30N"),
    ("23-24N", "30-32N"),
    ("24-26N", "26-28N"),
    ("24-26N", "28-30N"),
    ("24-26N", "30-32N"),
    ("26-28N", "28-30N"),
    ("26-28N", "30-32N"),
    ("28-30N", "30-32N"),
]

ADJACENT_PAIRS = {
    ("23-24N", "24-26N"),
    ("24-26N", "26-28N"),
    ("26-28N", "28-30N"),
    ("28-30N", "30-32N"),
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


def first_existing(paths: Sequence[Path], label: str, required: bool = True) -> Path | None:
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
    lower_to_original = {field.lower(): field for field in fields}
    for candidate in candidates:
        if candidate.lower() in lower_to_original:
            return lower_to_original[candidate.lower()]
    if required:
        raise ValueError(
            f"Could not identify {label}. Tried: {', '.join(candidates)}. "
            f"Available fields: {', '.join(fields)}"
        )
    return None


def read_csv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {path}")
        rows = [dict(row) for row in reader]
        return list(reader.fieldnames), rows


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def read_incidence_matrix(path: Path) -> dict[str, Any]:
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
        ["grid_cell_id", "cell_id"],
        "grid-cell ID field",
    )
    band_field = find_field(
        fields,
        ["centroid_latitude_band", "latitude_band"],
        "latitude-band field",
    )

    mapping: dict[str, str] = {}
    for row in rows:
        cell = row[cell_field].strip()
        band = row[band_field].strip()
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


EVIDENCE_TOKEN_RE = re.compile(
    r"(?<![A-Z0-9])(D1|D2|D3|D4|N0|C3)(?![A-Z0-9])", re.IGNORECASE
)


def parse_evidence_class(value: Any) -> str | None:
    text = str(value or "").strip().upper()
    if not text:
        return None
    tokens = {token.upper() for token in EVIDENCE_TOKEN_RE.findall(text)}
    if tokens == {"N0"}:
        return "N0"
    if tokens == {"C3"}:
        return "C3"
    if len(tokens) == 1:
        token = next(iter(tokens))
        if token in {"D1", "D2", "D3", "D4"}:
            return token
    normalized = re.sub(r"[^a-z0-9]+", "", text.casefold())
    if normalized in {"nonballooning", "fixednonballooning", "referencenonballooning", "noballooning"}:
        return "N0"
    if normalized in {"c3", "primaryc3", "d1d2d3", "d1tod3"}:
        return "C3"
    if normalized in {"d4excluded", "excludedd4"}:
        return "D4"
    return None


def choose_evidence_field(fields: Sequence[str], rows: Sequence[dict[str, str]]) -> str:
    preferred = {
        "evidence_class", "final_evidence_class", "final_evidence_category",
        "evidence_category", "d_level", "dlevel", "trait_class",
        "primary_class", "analysis_class", "ballooning_evidence_tier",
        "ballooning_evidence_category", "final_designation", "designation",
    }
    candidates: list[tuple[float, str]] = []
    for field in fields:
        parsed = [parse_evidence_class(row.get(field, "")) for row in rows]
        nonblank = [value for value in parsed if value is not None]
        if not nonblank:
            continue
        fraction = len(nonblank) / max(1, len(rows))
        classes = set(nonblank)
        name = re.sub(r"[^a-z0-9]+", "", field.casefold())
        bonus = 100.0 if field.casefold() in preferred else 0.0
        if any(term in name for term in ("evidence", "tier", "class", "designation", "decision")):
            bonus += 20.0
        if len(classes) >= 2 and fraction >= 0.25:
            candidates.append((bonus + 100.0 * fraction + 5.0 * len(classes), field))
    if not candidates:
        raise ValueError(
            "Trait table lacks an explicit D1/D2/D3/D4/N0 or C3/N0 evidence field. "
            "Legacy binary columns are not accepted because they cannot distinguish D4 from N0."
        )
    candidates.sort(reverse=True)
    return candidates[0][1]


def load_traits(path: Path, genera: Sequence[str]) -> dict[str, Any]:
    fields, rows = read_csv_rows(path)
    genus_field = find_field(fields, ["genus", "analysis_genus"], "trait genus field")
    evidence_field = choose_evidence_field(fields, rows)
    confidence_field = find_field(
        fields,
        ["final_confidence", "trait_final_confidence", "trait_confidence", "trait_ballooning_confidence"],
        "trait confidence field",
        required=False,
    )
    lookup: dict[str, dict[str, Any]] = {}
    for row in rows:
        genus = row.get(genus_field, "").strip()
        if not genus:
            continue
        key = genus.casefold()
        if key in lookup:
            raise ValueError(f"Duplicate genus {genus!r} in trait table {path}")
        lookup[key] = {
            "evidence": parse_evidence_class(row.get(evidence_field, "")),
            "confidence": (row.get(confidence_field, "").strip().upper() if confidence_field else "UNSPECIFIED") or "UNSPECIFIED",
        }
    c3_mask = n0_mask = d4_mask = low_confidence_mask = 0
    evidence_counts: dict[str, int] = defaultdict(int)
    confidence_counts: dict[str, int] = defaultdict(int)
    missing: list[str] = []
    unresolved: list[str] = []
    for genus_index, genus in enumerate(genera):
        trait = lookup.get(genus.casefold())
        if trait is None:
            missing.append(genus); continue
        evidence = trait["evidence"]
        if evidence is None:
            unresolved.append(genus); continue
        bit = 1 << genus_index
        if evidence in {"D1", "D2", "D3", "C3"}:
            c3_mask |= bit
        elif evidence == "N0":
            n0_mask |= bit
        else:
            d4_mask |= bit
        if trait["confidence"] == "LOW":
            low_confidence_mask |= bit
        evidence_counts[evidence] += 1
        confidence_counts[trait["confidence"]] += 1
    if missing:
        raise ValueError("Genera missing from trait table: " + ", ".join(missing[:20]))
    if unresolved:
        raise ValueError("Genera lacking explicit evidence classes: " + ", ".join(unresolved[:20]))
    all_mask = (1 << len(genera)) - 1
    classified_mask = c3_mask | n0_mask
    if (classified_mask | d4_mask) != all_mask or c3_mask & n0_mask:
        raise ValueError("C3, N0, and D4 trait masks are invalid.")
    return {
        "path": path,
        "evidence_field": evidence_field,
        "balloon_mask": c3_mask,
        "n0_mask": n0_mask,
        "d4_mask": d4_mask,
        "classified_mask": classified_mask,
        "low_confidence_mask": low_confidence_mask,
        "confidence_keep_mask": all_mask & ~low_confidence_mask,
        "evidence_counts": dict(sorted(evidence_counts.items())),
        "confidence_counts": dict(sorted(confidence_counts.items())),
    }

def union_masks(cell_masks: Sequence[int], selected_indices: Sequence[int]) -> int:
    result = 0
    for index in selected_indices:
        result |= cell_masks[index]
    return result


def band_metrics(
    genus_mask: int,
    selected_indices: Sequence[int],
    cell_masks: Sequence[int],
    balloon_mask: int,
    n0_mask: int,
    d4_mask: int,
    classified_mask: int,
    keep_mask: int | None = None,
) -> dict[str, float | int]:
    if keep_mask is not None:
        genus_mask &= keep_mask

    richness = genus_mask.bit_count()
    classified_richness = (genus_mask & classified_mask).bit_count()
    ballooning_richness = (genus_mask & balloon_mask).bit_count()
    non_ballooning_richness = (genus_mask & n0_mask).bit_count()
    excluded_d4_richness = (genus_mask & d4_mask).bit_count()

    genus_cell_incidences = 0
    classified_genus_cell_incidences = 0
    ballooning_genus_cell_incidences = 0
    non_ballooning_genus_cell_incidences = 0
    excluded_d4_genus_cell_incidences = 0

    for cell_index in selected_indices:
        cell_mask = cell_masks[cell_index]
        if keep_mask is not None:
            cell_mask &= keep_mask
        genus_cell_incidences += cell_mask.bit_count()
        classified_genus_cell_incidences += (cell_mask & classified_mask).bit_count()
        ballooning_genus_cell_incidences += (cell_mask & balloon_mask).bit_count()
        non_ballooning_genus_cell_incidences += (cell_mask & n0_mask).bit_count()
        excluded_d4_genus_cell_incidences += (cell_mask & d4_mask).bit_count()

    return {
        "genus_richness": richness,
        "classified_genus_richness": classified_richness,
        "ballooning_genus_richness": ballooning_richness,
        "non_ballooning_genus_richness": non_ballooning_richness,
        "excluded_D4_genus_richness": excluded_d4_richness,
        "ballooning_genus_proportion": (
            ballooning_richness / classified_richness if classified_richness else math.nan
        ),
        "genus_cell_incidences": genus_cell_incidences,
        "classified_genus_cell_incidences": classified_genus_cell_incidences,
        "ballooning_genus_cell_incidences": ballooning_genus_cell_incidences,
        "non_ballooning_genus_cell_incidences": non_ballooning_genus_cell_incidences,
        "excluded_D4_genus_cell_incidences": excluded_d4_genus_cell_incidences,
        "ballooning_incidence_proportion": (
            ballooning_genus_cell_incidences / classified_genus_cell_incidences
            if classified_genus_cell_incidences else math.nan
        ),
    }

def beta_metrics(mask_a: int, mask_b: int) -> dict[str, float | int]:
    shared = (mask_a & mask_b).bit_count()
    unique_a = (mask_a & ~mask_b).bit_count()
    unique_b = (mask_b & ~mask_a).bit_count()

    jaccard_denominator = shared + unique_a + unique_b
    sorensen_denominator = 2 * shared + unique_a + unique_b
    minimum_unique = min(unique_a, unique_b)
    simpson_denominator = shared + minimum_unique

    jaccard = (
        (unique_a + unique_b) / jaccard_denominator
        if jaccard_denominator
        else math.nan
    )
    sorensen = (
        (unique_a + unique_b) / sorensen_denominator
        if sorensen_denominator
        else math.nan
    )
    simpson = (
        minimum_unique / simpson_denominator
        if simpson_denominator
        else math.nan
    )
    nestedness = (
        max(0.0, sorensen - simpson)
        if not math.isnan(sorensen) and not math.isnan(simpson)
        else math.nan
    )

    return {
        "shared_genera": shared,
        "unique_to_band_1": unique_a,
        "unique_to_band_2": unique_b,
        "jaccard_dissimilarity": jaccard,
        "sorensen_dissimilarity": sorensen,
        "simpson_turnover": simpson,
        "sorensen_nestedness_resultant": nestedness,
    }


def percentile(values: Sequence[float], probability: float) -> float:
    ordered = sorted(values)
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


def summarize_values(values: Sequence[float]) -> dict[str, float | int]:
    clean = [
        float(value)
        for value in values
        if value is not None and not math.isnan(float(value))
    ]

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
        "standard_deviation": (
            statistics.stdev(clean) if len(clean) > 1 else 0.0
        ),
        "p025": percentile(clean, 0.025),
        "p975": percentile(clean, 0.975),
        "minimum": min(clean),
        "maximum": max(clean),
    }


def summarize_long(
    rows: Sequence[dict[str, Any]],
    group_fields: Sequence[str],
    metric_fields: Sequence[str],
) -> list[dict[str, Any]]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        key = tuple(row[field] for field in group_fields)
        grouped[key].append(row)

    output: list[dict[str, Any]] = []

    for key in sorted(grouped, key=lambda item: tuple(str(value) for value in item)):
        group_rows = grouped[key]
        base = dict(zip(group_fields, key))

        for metric in metric_fields:
            summary = summarize_values(
                [float(row[metric]) for row in group_rows]
            )
            output.append(
                {
                    **base,
                    "metric": metric,
                    **summary,
                }
            )

    return output


def paired_difference_rows(
    rows: Sequence[dict[str, Any]],
    comparison_analysis: str,
    group_fields: Sequence[str],
    metric_fields: Sequence[str],
) -> list[dict[str, Any]]:
    primary_lookup: dict[tuple[Any, ...], dict[str, Any]] = {}
    comparison_lookup: dict[tuple[Any, ...], dict[str, Any]] = {}

    key_fields = ["iteration", *group_fields]

    for row in rows:
        key = tuple(row[field] for field in key_fields)
        if row["analysis"] == "primary":
            primary_lookup[key] = row
        elif row["analysis"] == comparison_analysis:
            comparison_lookup[key] = row

    output: list[dict[str, Any]] = []

    for key in sorted(set(primary_lookup) & set(comparison_lookup)):
        primary = primary_lookup[key]
        comparison = comparison_lookup[key]
        base = dict(zip(key_fields, key))

        for metric in metric_fields:
            primary_value = float(primary[metric])
            comparison_value = float(comparison[metric])
            if math.isnan(primary_value) or math.isnan(comparison_value):
                difference = math.nan
            else:
                difference = comparison_value - primary_value

            output.append(
                {
                    **base,
                    "comparison": f"{comparison_analysis}_minus_primary",
                    "metric": metric,
                    "primary_value": primary_value,
                    "comparison_value": comparison_value,
                    "difference": difference,
                }
            )

    return output


def summarize_differences(
    difference_rows: Sequence[dict[str, Any]],
    group_fields: Sequence[str],
) -> list[dict[str, Any]]:
    grouped: dict[tuple[Any, ...], list[float]] = defaultdict(list)

    for row in difference_rows:
        key = tuple(row[field] for field in group_fields)
        grouped[key].append(float(row["difference"]))

    output: list[dict[str, Any]] = []

    for key in sorted(grouped, key=lambda item: tuple(str(value) for value in item)):
        summary = summarize_values(grouped[key])
        output.append(
            {
                **dict(zip(group_fields, key)),
                **summary,
            }
        )

    return output


def svg_escape(text: Any) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def interval_plot_svg(
    path: Path,
    title: str,
    y_label: str,
    categories: Sequence[str],
    series: Sequence[dict[str, Any]],
    y_min: float | None = None,
    y_max: float | None = None,
) -> None:
    width = 1000
    height = 650
    left = 105
    right = 35
    top = 85
    bottom = 110
    plot_width = width - left - right
    plot_height = height - top - bottom

    all_values: list[float] = []
    for item in series:
        all_values.extend([float(item["p025"]), float(item["p975"])])

    actual_min = min(all_values) if all_values else 0.0
    actual_max = max(all_values) if all_values else 1.0

    if y_min is None:
        padding = (actual_max - actual_min) * 0.08 or 1.0
        y_min = actual_min - padding
    if y_max is None:
        padding = (actual_max - actual_min) * 0.08 or 1.0
        y_max = actual_max + padding
    if y_max <= y_min:
        y_max = y_min + 1.0

    def x_position(category_index: int, series_index: int, series_count: int) -> float:
        base = left + plot_width * (category_index + 0.5) / len(categories)
        if series_count == 1:
            return base
        spread = min(55.0, plot_width / max(1, len(categories)) * 0.28)
        return base + (series_index - (series_count - 1) / 2) * spread

    def y_position(value: float) -> float:
        return top + plot_height * (y_max - value) / (y_max - y_min)

    by_category: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for item in series:
        by_category[str(item["category"])].append(item)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>'
        'text{font-family:Arial,Helvetica,sans-serif;fill:#111}'
        '.axis{stroke:#111;stroke-width:1.5}'
        '.grid{stroke:#ddd;stroke-width:1}'
        '.interval{stroke:#333;stroke-width:3}'
        '.point{fill:#111}'
        '.median{fill:white;stroke:#111;stroke-width:2}'
        '</style>',
        f'<text x="{width/2}" y="38" text-anchor="middle" '
        f'font-size="24" font-weight="bold">{svg_escape(title)}</text>',
    ]

    tick_count = 5
    for tick in range(tick_count + 1):
        value = y_min + (y_max - y_min) * tick / tick_count
        y = y_position(value)
        parts.append(
            f'<line class="grid" x1="{left}" x2="{width-right}" y1="{y}" y2="{y}"/>'
        )
        parts.append(
            f'<text x="{left-12}" y="{y+5}" text-anchor="end" font-size="15">'
            f'{value:.2f}</text>'
        )

    parts.extend(
        [
            f'<line class="axis" x1="{left}" x2="{left}" y1="{top}" y2="{height-bottom}"/>',
            f'<line class="axis" x1="{left}" x2="{width-right}" '
            f'y1="{height-bottom}" y2="{height-bottom}"/>',
            f'<text transform="translate(28,{top+plot_height/2}) rotate(-90)" '
            f'text-anchor="middle" font-size="18">{svg_escape(y_label)}</text>',
        ]
    )

    legend_labels: list[str] = []
    for category_index, category in enumerate(categories):
        category_series = by_category.get(category, [])
        for series_index, item in enumerate(category_series):
            x = x_position(category_index, series_index, len(category_series))
            y_low = y_position(float(item["p025"]))
            y_high = y_position(float(item["p975"]))
            y_mean = y_position(float(item["mean"]))
            y_median = y_position(float(item["median"]))

            parts.append(
                f'<line class="interval" x1="{x}" x2="{x}" y1="{y_low}" y2="{y_high}"/>'
            )
            parts.append(
                f'<line class="interval" x1="{x-7}" x2="{x+7}" '
                f'y1="{y_low}" y2="{y_low}"/>'
            )
            parts.append(
                f'<line class="interval" x1="{x-7}" x2="{x+7}" '
                f'y1="{y_high}" y2="{y_high}"/>'
            )
            parts.append(
                f'<circle class="point" cx="{x}" cy="{y_mean}" r="6"/>'
            )
            parts.append(
                f'<circle class="median" cx="{x}" cy="{y_median}" r="3"/>'
            )

            label = str(item.get("series", ""))
            if label and label not in legend_labels:
                legend_labels.append(label)

        x_label = left + plot_width * (category_index + 0.5) / len(categories)
        parts.append(
            f'<text x="{x_label}" y="{height-bottom+34}" text-anchor="middle" '
            f'font-size="16">{svg_escape(category)}</text>'
        )

    parts.append(
        f'<text x="{width/2}" y="{height-30}" text-anchor="middle" font-size="14">'
        'Point = mean; open center = median; line = 95% equal-cell subset interval'
        '</text>'
    )

    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def heatmap_svg(
    path: Path,
    title: str,
    matrix: dict[tuple[str, str], float],
) -> None:
    width = 800
    height = 760
    left = 150
    top = 100
    cell_size = 105

    def shade(value: float) -> str:
        clipped = max(0.0, min(1.0, value))
        level = int(round(255 * (1 - clipped)))
        return f"rgb({level},{level},{level})"

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#111}</style>',
        f'<text x="{width/2}" y="42" text-anchor="middle" '
        f'font-size="24" font-weight="bold">{svg_escape(title)}</text>',
    ]

    for row_index, band_a in enumerate(BAND_ORDER):
        y = top + row_index * cell_size
        parts.append(
            f'<text x="{left-15}" y="{y+cell_size/2+6}" text-anchor="end" '
            f'font-size="16">{svg_escape(BAND_LABELS[band_a])}</text>'
        )

        for column_index, band_b in enumerate(BAND_ORDER):
            x = left + column_index * cell_size
            if band_a == band_b:
                value = 0.0
            else:
                key = (band_a, band_b)
                reverse_key = (band_b, band_a)
                value = matrix.get(key, matrix.get(reverse_key, math.nan))

            fill = "#f2f2f2" if math.isnan(value) else shade(value)
            text_fill = "white" if not math.isnan(value) and value > 0.55 else "black"
            label = "" if math.isnan(value) else f"{value:.2f}"

            parts.append(
                f'<rect x="{x}" y="{y}" width="{cell_size}" height="{cell_size}" '
                f'fill="{fill}" stroke="white" stroke-width="2"/>'
            )
            parts.append(
                f'<text x="{x+cell_size/2}" y="{y+cell_size/2+6}" '
                f'text-anchor="middle" font-size="18" fill="{text_fill}">'
                f'{label}</text>'
            )

    for column_index, band in enumerate(BAND_ORDER):
        x = left + column_index * cell_size + cell_size / 2
        y = top - 12
        parts.append(
            f'<text x="{x}" y="{y}" text-anchor="end" font-size="16" '
            f'transform="rotate(-45 {x} {y})">{svg_escape(BAND_LABELS[band])}</text>'
        )

    parts.append(
        f'<text x="{width/2}" y="{height-35}" text-anchor="middle" font-size="14">'
        'Cells show median equal-cell Jaccard dissimilarity'
        '</text>'
    )
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def archive_existing(output_dir: Path, archive_root: Path) -> Path | None:
    if not output_dir.exists() or not any(output_dir.iterdir()):
        return None

    timestamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    archive_dir = archive_root / f"11B_equal_cell_resampling_{timestamp}"
    archive_dir.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(output_dir, archive_dir)
    shutil.rmtree(output_dir)
    return archive_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Equal-cell Monte Carlo resampling of Baja arachnid richness, "
            "ballooning composition, and beta diversity."
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
        help=f"Number of resampling iterations; default {DEFAULT_ITERATIONS}.",
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
    step11_dir = (
        project_root
        / "04_analysis"
        / "11_latitude_band_diversity_turnover"
    )
    trait_fallback = project_root / "02_data_clean" / "07_final_trait_merge"

    output_dir = (
        project_root
        / "04_analysis"
        / "11B_equal_cell_resampling"
    )
    figure_dir = output_dir / "figures"
    archive_root = project_root / "08_archive"

    archived = archive_existing(output_dir, archive_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    log_path = output_dir / "11B_analysis_log.txt"
    with log_path.open("w", encoding="utf-8") as log_handle:
        original_stdout = sys.stdout
        original_stderr = sys.stderr
        sys.stdout = Tee(original_stdout, log_handle)
        sys.stderr = Tee(original_stderr, log_handle)

        try:
            print("STEP 11B STARTED")
            print(f"Version: {SCRIPT_VERSION}")
            print(f"Project: {project_root}")
            print(f"Iterations: {args.iterations}")
            print(f"Random seed: {args.seed}")
            if archived is not None:
                print(f"Archived prior Step 11B output: {archived}")

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
                    / "07_reviewed_genus_trait_lookup_normalized.csv",
                    analysis_ready
                    / "03_trait_tables"
                    / "07_reviewed_genus_trait_lookup_final.csv",
                    step11_dir
                    / "11_biodiversity_final_genus_by_latitude_band_incidence.csv",
                    trait_fallback
                    / "07_reviewed_genus_trait_lookup_final.csv",
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

            if len(genera) != 267:
                print(
                    f"NOTICE: primary matrix contains {len(genera)} genera, "
                    "not the previously documented 267."
                )
            if len(cells) != 205:
                print(
                    f"NOTICE: primary matrix contains {len(cells)} cells, "
                    "not the previously documented 205."
                )
            if any(genus.casefold() == "fesa" for genus in genera):
                raise ValueError("Fesa remains in the primary genus universe.")

            primary_cell_masks = align_matrix(primary_raw, genera, cells)
            strict_cell_masks = align_matrix(strict_raw, genera, cells)

            cell_bands = load_cell_bands(cell_lookup_path, cells)
            traits = load_traits(trait_path, genera)

            band_indices: dict[str, list[int]] = {
                band: [
                    index
                    for index, cell in enumerate(cells)
                    if cell_bands[cell] == band
                ]
                for band in BAND_ORDER
            }

            cell_counts = {
                band: len(indices)
                for band, indices in band_indices.items()
            }
            if any(count == 0 for count in cell_counts.values()):
                raise ValueError(f"At least one latitude band has no cells: {cell_counts}")

            sample_size = min(cell_counts.values())
            if sample_size != 22:
                print(
                    f"NOTICE: equal-cell sample size is {sample_size}; "
                    "the previously documented minimum was 22."
                )

            print("Cells by latitude band:")
            for band in BAND_ORDER:
                print(f"  {band}: {cell_counts[band]}")
            print(f"Equal-cell sample size: {sample_size}")
            print("Trait evidence field: " + traits["evidence_field"])
            print("Trait evidence counts: " + ", ".join(f"{k}={v}" for k, v in traits["evidence_counts"].items()))
            print(
                "Trait confidence counts: "
                + ", ".join(
                    f"{key}={value}"
                    for key, value in traits["confidence_counts"].items()
                )
            )

            rng = random.Random(args.seed)

            selection_counts = [0] * len(cells)
            band_iteration_rows: list[dict[str, Any]] = []
            beta_iteration_rows: list[dict[str, Any]] = []

            for iteration in range(1, args.iterations + 1):
                selected_by_band: dict[str, list[int]] = {}

                for band in BAND_ORDER:
                    selected = rng.sample(band_indices[band], sample_size)
                    selected_by_band[band] = selected
                    for cell_index in selected:
                        selection_counts[cell_index] += 1

                primary_masks = {
                    band: union_masks(primary_cell_masks, selected_by_band[band])
                    for band in BAND_ORDER
                }
                strict_masks = {
                    band: union_masks(strict_cell_masks, selected_by_band[band])
                    for band in BAND_ORDER
                }
                low_conf_masks = {
                    band: primary_masks[band] & traits["confidence_keep_mask"]
                    for band in BAND_ORDER
                }

                masks_by_analysis = {
                    "primary": primary_masks,
                    "taxonomy_strict": strict_masks,
                    "explicit_low_confidence_exclusion": low_conf_masks,
                }

                for analysis in ANALYSIS_ORDER:
                    analysis_masks = masks_by_analysis[analysis]

                    for band in BAND_ORDER:
                        keep_mask = (
                            traits["confidence_keep_mask"]
                            if analysis == "explicit_low_confidence_exclusion"
                            else None
                        )
                        source_cell_masks = (
                            strict_cell_masks
                            if analysis == "taxonomy_strict"
                            else primary_cell_masks
                        )
                        metrics = band_metrics(
                            genus_mask=analysis_masks[band],
                            selected_indices=selected_by_band[band],
                            cell_masks=source_cell_masks,
                            balloon_mask=traits["balloon_mask"],
                            n0_mask=traits["n0_mask"],
                            d4_mask=traits["d4_mask"],
                            classified_mask=traits["classified_mask"],
                            keep_mask=keep_mask,
                        )

                        band_iteration_rows.append(
                            {
                                "iteration": iteration,
                                "analysis": analysis,
                                "latitude_band": band,
                                "latitude_band_label": BAND_LABELS[band],
                                "sampled_cells": sample_size,
                                **metrics,
                            }
                        )

                    for band_a, band_b in PAIR_ORDER:
                        metrics = beta_metrics(
                            analysis_masks[band_a],
                            analysis_masks[band_b],
                        )
                        beta_iteration_rows.append(
                            {
                                "iteration": iteration,
                                "analysis": analysis,
                                "band_1": band_a,
                                "band_2": band_b,
                                "band_1_label": BAND_LABELS[band_a],
                                "band_2_label": BAND_LABELS[band_b],
                                "adjacent_bands": (
                                    (band_a, band_b) in ADJACENT_PAIRS
                                ),
                                "sampled_cells_per_band": sample_size,
                                **metrics,
                            }
                        )

                if iteration % max(1, args.iterations // 10) == 0:
                    print(
                        f"Completed {iteration:,} of {args.iterations:,} iterations"
                    )

            band_fields = [
                "iteration",
                "analysis",
                "latitude_band",
                "latitude_band_label",
                "sampled_cells",
                "genus_richness",
                "classified_genus_richness",
                "ballooning_genus_richness",
                "non_ballooning_genus_richness",
                "excluded_D4_genus_richness",
                "ballooning_genus_proportion",
                "genus_cell_incidences",
                "classified_genus_cell_incidences",
                "ballooning_genus_cell_incidences",
                "non_ballooning_genus_cell_incidences",
                "excluded_D4_genus_cell_incidences",
                "ballooning_incidence_proportion",
            ]
            beta_fields = [
                "iteration",
                "analysis",
                "band_1",
                "band_2",
                "band_1_label",
                "band_2_label",
                "adjacent_bands",
                "sampled_cells_per_band",
                "shared_genera",
                "unique_to_band_1",
                "unique_to_band_2",
                "jaccard_dissimilarity",
                "sorensen_dissimilarity",
                "simpson_turnover",
                "sorensen_nestedness_resultant",
            ]

            write_csv(
                output_dir / "11B_equal_cell_band_metrics_iterations.csv",
                band_iteration_rows,
                band_fields,
            )
            write_csv(
                output_dir / "11B_equal_cell_pairwise_beta_iterations.csv",
                beta_iteration_rows,
                beta_fields,
            )

            band_metric_fields = [
                "genus_richness",
                "classified_genus_richness",
                "ballooning_genus_richness",
                "non_ballooning_genus_richness",
                "excluded_D4_genus_richness",
                "ballooning_genus_proportion",
                "genus_cell_incidences",
                "classified_genus_cell_incidences",
                "ballooning_genus_cell_incidences",
                "non_ballooning_genus_cell_incidences",
                "excluded_D4_genus_cell_incidences",
                "ballooning_incidence_proportion",
            ]
            beta_metric_fields = [
                "shared_genera",
                "unique_to_band_1",
                "unique_to_band_2",
                "jaccard_dissimilarity",
                "sorensen_dissimilarity",
                "simpson_turnover",
                "sorensen_nestedness_resultant",
            ]

            band_summary = summarize_long(
                band_iteration_rows,
                ["analysis", "latitude_band", "latitude_band_label"],
                band_metric_fields,
            )
            beta_summary = summarize_long(
                beta_iteration_rows,
                [
                    "analysis",
                    "band_1",
                    "band_2",
                    "band_1_label",
                    "band_2_label",
                    "adjacent_bands",
                ],
                beta_metric_fields,
            )

            summary_fields = [
                "analysis",
                "latitude_band",
                "latitude_band_label",
                "metric",
                "n_iterations",
                "mean",
                "median",
                "standard_deviation",
                "p025",
                "p975",
                "minimum",
                "maximum",
            ]
            beta_summary_fields = [
                "analysis",
                "band_1",
                "band_2",
                "band_1_label",
                "band_2_label",
                "adjacent_bands",
                "metric",
                "n_iterations",
                "mean",
                "median",
                "standard_deviation",
                "p025",
                "p975",
                "minimum",
                "maximum",
            ]

            write_csv(
                output_dir / "11B_equal_cell_band_metrics_summary.csv",
                band_summary,
                summary_fields,
            )
            write_csv(
                output_dir / "11B_equal_cell_pairwise_beta_summary.csv",
                beta_summary,
                beta_summary_fields,
            )

            adjacent_summary = [
                row for row in beta_summary if str(row["adjacent_bands"]) == "True"
            ]
            write_csv(
                output_dir / "11B_equal_cell_adjacent_band_turnover_summary.csv",
                adjacent_summary,
                beta_summary_fields,
            )

            band_difference_rows: list[dict[str, Any]] = []
            beta_difference_rows: list[dict[str, Any]] = []

            for comparison in [
                "taxonomy_strict",
                "explicit_low_confidence_exclusion",
            ]:
                band_difference_rows.extend(
                    paired_difference_rows(
                        band_iteration_rows,
                        comparison,
                        ["latitude_band", "latitude_band_label"],
                        band_metric_fields,
                    )
                )
                beta_difference_rows.extend(
                    paired_difference_rows(
                        beta_iteration_rows,
                        comparison,
                        [
                            "band_1",
                            "band_2",
                            "band_1_label",
                            "band_2_label",
                            "adjacent_bands",
                        ],
                        beta_metric_fields,
                    )
                )

            band_difference_summary = summarize_differences(
                band_difference_rows,
                [
                    "comparison",
                    "latitude_band",
                    "latitude_band_label",
                    "metric",
                ],
            )
            beta_difference_summary = summarize_differences(
                beta_difference_rows,
                [
                    "comparison",
                    "band_1",
                    "band_2",
                    "band_1_label",
                    "band_2_label",
                    "adjacent_bands",
                    "metric",
                ],
            )

            band_difference_fields = [
                "comparison",
                "latitude_band",
                "latitude_band_label",
                "metric",
                "n_iterations",
                "mean",
                "median",
                "standard_deviation",
                "p025",
                "p975",
                "minimum",
                "maximum",
            ]
            beta_difference_fields = [
                "comparison",
                "band_1",
                "band_2",
                "band_1_label",
                "band_2_label",
                "adjacent_bands",
                "metric",
                "n_iterations",
                "mean",
                "median",
                "standard_deviation",
                "p025",
                "p975",
                "minimum",
                "maximum",
            ]

            write_csv(
                output_dir
                / "11B_primary_vs_sensitivity_band_difference_summary.csv",
                band_difference_summary,
                band_difference_fields,
            )
            write_csv(
                output_dir
                / "11B_primary_vs_sensitivity_beta_difference_summary.csv",
                beta_difference_summary,
                beta_difference_fields,
            )

            selection_rows = []
            for cell_index, cell in enumerate(cells):
                band = cell_bands[cell]
                expected_probability = sample_size / cell_counts[band]
                selection_rows.append(
                    {
                        "grid_cell_id": cell,
                        "latitude_band": band,
                        "latitude_band_label": BAND_LABELS[band],
                        "times_selected": selection_counts[cell_index],
                        "selection_proportion": (
                            selection_counts[cell_index] / args.iterations
                        ),
                        "expected_selection_proportion": expected_probability,
                        "difference_from_expected": (
                            selection_counts[cell_index] / args.iterations
                            - expected_probability
                        ),
                    }
                )

            write_csv(
                output_dir / "11B_cell_selection_frequencies.csv",
                selection_rows,
                [
                    "grid_cell_id",
                    "latitude_band",
                    "latitude_band_label",
                    "times_selected",
                    "selection_proportion",
                    "expected_selection_proportion",
                    "difference_from_expected",
                ],
            )

            full_band_rows: list[dict[str, Any]] = []
            full_beta_rows: list[dict[str, Any]] = []

            full_selected_by_band = band_indices
            full_primary_masks = {
                band: union_masks(primary_cell_masks, full_selected_by_band[band])
                for band in BAND_ORDER
            }
            full_strict_masks = {
                band: union_masks(strict_cell_masks, full_selected_by_band[band])
                for band in BAND_ORDER
            }
            full_low_masks = {
                band: full_primary_masks[band] & traits["confidence_keep_mask"]
                for band in BAND_ORDER
            }

            full_masks_by_analysis = {
                "primary": full_primary_masks,
                "taxonomy_strict": full_strict_masks,
                "explicit_low_confidence_exclusion": full_low_masks,
            }

            for analysis in ANALYSIS_ORDER:
                analysis_masks = full_masks_by_analysis[analysis]
                for band in BAND_ORDER:
                    keep_mask = (
                        traits["confidence_keep_mask"]
                        if analysis == "explicit_low_confidence_exclusion"
                        else None
                    )
                    source_cell_masks = (
                        strict_cell_masks
                        if analysis == "taxonomy_strict"
                        else primary_cell_masks
                    )
                    full_band_rows.append(
                        {
                            "analysis": analysis,
                            "latitude_band": band,
                            "latitude_band_label": BAND_LABELS[band],
                            "available_cells": cell_counts[band],
                            **band_metrics(
                                genus_mask=analysis_masks[band],
                                selected_indices=full_selected_by_band[band],
                                cell_masks=source_cell_masks,
                                balloon_mask=traits["balloon_mask"],
                                n0_mask=traits["n0_mask"],
                                d4_mask=traits["d4_mask"],
                                classified_mask=traits["classified_mask"],
                                keep_mask=keep_mask,
                            ),
                        }
                    )

                for band_a, band_b in PAIR_ORDER:
                    full_beta_rows.append(
                        {
                            "analysis": analysis,
                            "band_1": band_a,
                            "band_2": band_b,
                            "band_1_label": BAND_LABELS[band_a],
                            "band_2_label": BAND_LABELS[band_b],
                            "adjacent_bands": (
                                (band_a, band_b) in ADJACENT_PAIRS
                            ),
                            **beta_metrics(
                                analysis_masks[band_a],
                                analysis_masks[band_b],
                            ),
                        }
                    )

            write_csv(
                output_dir / "11B_full_cell_reference_band_metrics.csv",
                full_band_rows,
                [
                    "analysis",
                    "latitude_band",
                    "latitude_band_label",
                    "available_cells",
                    *band_metric_fields,
                ],
            )
            write_csv(
                output_dir / "11B_full_cell_reference_pairwise_beta.csv",
                full_beta_rows,
                [
                    "analysis",
                    "band_1",
                    "band_2",
                    "band_1_label",
                    "band_2_label",
                    "adjacent_bands",
                    *beta_metric_fields,
                ],
            )

            validation_rows = [
                {
                    "check": "primary_matrix_genera_unique",
                    "passed": len({genus.casefold() for genus in genera}) == len(genera),
                    "detail": f"{len(genera)} genera",
                },
                {
                    "check": "primary_matrix_cells_unique",
                    "passed": len(set(cells)) == len(cells),
                    "detail": f"{len(cells)} cells",
                },
                {
                    "check": "fesa_absent",
                    "passed": not any(
                        genus.casefold() == "fesa" for genus in genera
                    ),
                    "detail": "",
                },
                {
                    "check": "five_latitude_bands_present",
                    "passed": set(cell_counts) == set(BAND_ORDER),
                    "detail": json.dumps(cell_counts, sort_keys=True),
                },
                {
                    "check": "all_bands_have_at_least_equal_cell_sample",
                    "passed": all(
                        count >= sample_size for count in cell_counts.values()
                    ),
                    "detail": f"equal-cell sample={sample_size}",
                },
                {
                    "check": "all_traits_resolved",
                    "passed": (
                        traits["balloon_mask"]
                        | ((1 << len(genera)) - 1 ^ traits["balloon_mask"])
                    ).bit_count()
                    == len(genera),
                    "detail": json.dumps(
                        traits["confidence_counts"], sort_keys=True
                    ),
                },
                {
                    "check": "expected_band_iteration_rows",
                    "passed": len(band_iteration_rows)
                    == args.iterations * len(BAND_ORDER) * len(ANALYSIS_ORDER),
                    "detail": str(len(band_iteration_rows)),
                },
                {
                    "check": "expected_beta_iteration_rows",
                    "passed": len(beta_iteration_rows)
                    == args.iterations * len(PAIR_ORDER) * len(ANALYSIS_ORDER),
                    "detail": str(len(beta_iteration_rows)),
                },
                {
                    "check": "each_band_draw_has_equal_cell_count",
                    "passed": all(
                        int(row["sampled_cells"]) == sample_size
                        for row in band_iteration_rows
                    ),
                    "detail": str(sample_size),
                },
                {
                    "check": "smallest_band_cells_always_selected",
                    "passed": all(
                        selection_counts[index] == args.iterations
                        for band in BAND_ORDER
                        if cell_counts[band] == sample_size
                        for index in band_indices[band]
                    ),
                    "detail": "",
                },
            ]

            if not all(bool(row["passed"]) for row in validation_rows):
                write_csv(
                    output_dir / "11B_validation.csv",
                    validation_rows,
                    ["check", "passed", "detail"],
                )
                failed = [
                    row["check"] for row in validation_rows if not row["passed"]
                ]
                raise RuntimeError(
                    "Step 11B validation failed: " + ", ".join(failed)
                )

            write_csv(
                output_dir / "11B_validation.csv",
                validation_rows,
                ["check", "passed", "detail"],
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
                output_dir / "11B_input_file_manifest.csv",
                input_manifest,
                ["path", "bytes", "sha256"],
            )

            # -------------------------- figures -----------------------------

            def select_summary(
                metric: str,
                analyses: set[str],
            ) -> list[dict[str, Any]]:
                selected = [
                    row
                    for row in band_summary
                    if row["metric"] == metric
                    and row["analysis"] in analyses
                ]
                selected.sort(
                    key=lambda row: (
                        BAND_ORDER.index(row["latitude_band"]),
                        ANALYSIS_ORDER.index(row["analysis"]),
                    )
                )
                return selected

            richness_series = []
            for row in select_summary("genus_richness", {"primary"}):
                richness_series.append(
                    {
                        "category": BAND_LABELS[row["latitude_band"]],
                        "series": "Primary",
                        **row,
                    }
                )

            interval_plot_svg(
                figure_dir / "11B_equal_cell_genus_richness.svg",
                "Equal-cell standardized arachnid genus richness",
                "Genus richness in 22 occupied cells",
                [BAND_LABELS[band] for band in BAND_ORDER],
                richness_series,
                y_min=0,
            )

            balloon_series = []
            for row in select_summary(
                "ballooning_genus_proportion",
                {"primary", "explicit_low_confidence_exclusion"},
            ):
                balloon_series.append(
                    {
                        "category": BAND_LABELS[row["latitude_band"]],
                        "series": (
                            "Primary"
                            if row["analysis"] == "primary"
                            else "Explicit LOW-confidence exclusion"
                        ),
                        **row,
                    }
                )

            interval_plot_svg(
                figure_dir / "11B_equal_cell_ballooning_proportion.svg",
                "Equal-cell ballooning composition",
                "Proportion of observed genera classified as ballooning",
                [BAND_LABELS[band] for band in BAND_ORDER],
                balloon_series,
                y_min=0,
                y_max=1,
            )

            adjacent_turnover_series = []
            for row in beta_summary:
                if (
                    row["analysis"] == "primary"
                    and row["metric"] == "simpson_turnover"
                    and str(row["adjacent_bands"]) == "True"
                ):
                    adjacent_turnover_series.append(
                        {
                            "category": (
                                f"{BAND_LABELS[row['band_1']]} to "
                                f"{BAND_LABELS[row['band_2']]}"
                            ),
                            "series": "Primary",
                            **row,
                        }
                    )

            adjacent_categories = [
                (
                    f"{BAND_LABELS[band_a]} to "
                    f"{BAND_LABELS[band_b]}"
                )
                for band_a, band_b in PAIR_ORDER
                if (band_a, band_b) in ADJACENT_PAIRS
            ]

            interval_plot_svg(
                figure_dir / "11B_equal_cell_adjacent_simpson_turnover.svg",
                "Equal-cell adjacent-band faunal turnover",
                "Simpson turnover",
                adjacent_categories,
                adjacent_turnover_series,
                y_min=0,
                y_max=1,
            )

            jaccard_matrix: dict[tuple[str, str], float] = {}
            for row in beta_summary:
                if (
                    row["analysis"] == "primary"
                    and row["metric"] == "jaccard_dissimilarity"
                ):
                    jaccard_matrix[(row["band_1"], row["band_2"])] = float(
                        row["median"]
                    )

            heatmap_svg(
                figure_dir / "11B_equal_cell_jaccard_heatmap.svg",
                "Equal-cell pairwise Jaccard dissimilarity",
                jaccard_matrix,
            )

            provenance = {
                "created_utc": utc_now(),
                "script_version": SCRIPT_VERSION,
                "script_path": str(Path(__file__).resolve()),
                "project_root": str(project_root),
                "output_dir": str(output_dir),
                "iterations": args.iterations,
                "random_seed": args.seed,
                "sampling_design": {
                    "sampling_unit": "occupied 25-km equal-area grid cell",
                    "sampling_method": (
                        "Monte Carlo subsampling without replacement within "
                        "each latitude band"
                    ),
                    "same_cell_draws_used_for_all_analyses": True,
                    "equal_cells_per_band": sample_size,
                    "available_cells_by_band": cell_counts,
                    "latitude_band_order": BAND_ORDER,
                },
                "analyses": {
                    "primary": (
                        "Primary 267-genus biodiversity incidence matrix"
                    ),
                    "taxonomy_strict": (
                        "Taxonomy-strict matrix evaluated with the same "
                        "cell draw in every iteration"
                    ),
                    "explicit_low_confidence_exclusion": (
                        "Primary matrix after excluding genera explicitly "
                        "classified LOW confidence; UNSPECIFIED legacy "
                        "assignments are retained"
                    ),
                },
                "primary_trait_definition": "C3 = D1 + D2 + D3 versus fixed N0; D4 excluded",
                "trait_evidence_field": traits["evidence_field"],
                "trait_evidence_counts": traits["evidence_counts"],
                "trait_confidence_counts": traits["confidence_counts"],
                "beta_diversity_formulas": {
                    "jaccard_dissimilarity": "(b+c)/(a+b+c)",
                    "sorensen_dissimilarity": "(b+c)/(2a+b+c)",
                    "simpson_turnover": "min(b,c)/(a+min(b,c))",
                    "sorensen_nestedness_resultant": (
                        "sorensen_dissimilarity - simpson_turnover"
                    ),
                },
                "interval_definition": (
                    "2.5th and 97.5th percentiles of the Monte Carlo "
                    "resampling distribution"
                ),
                "input_manifest": input_manifest,
                "validation_passed": True,
            }
            (output_dir / "11B_provenance.json").write_text(
                json.dumps(provenance, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )

            running_script = Path(__file__).resolve()
            shutil.copy2(
                running_script,
                output_dir / "11B_equal_cell_resampling.py",
            )

            output_files = sorted(
                path
                for path in output_dir.rglob("*")
                if path.is_file()
                and path.name != "11B_output_file_manifest.csv"
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
                output_dir / "11B_output_file_manifest.csv",
                output_manifest,
                ["relative_path", "bytes", "sha256"],
            )

            print()
            print("=" * 78)
            print("STEP 11B COMPLETED SUCCESSFULLY")
            print("=" * 78)
            print(f"Iterations: {args.iterations:,}")
            print(f"Equal cells per latitude band: {sample_size}")
            print(f"Band-metric iteration rows: {len(band_iteration_rows):,}")
            print(f"Pairwise-beta iteration rows: {len(beta_iteration_rows):,}")
            print("Analyses: primary, taxonomy_strict, explicit LOW-confidence exclusion")
            print(f"Outputs: {output_dir}")

            return 0

        finally:
            sys.stdout = original_stdout
            sys.stderr = original_stderr


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nStep 11B interrupted by user.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:
        print(f"\nSTEP 11B FAILED: {exc}", file=sys.stderr)
        raise
