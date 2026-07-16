#!/usr/bin/env python3
"""
STEP 11G — Equal-cell, trait-partitioned richness and turnover
Baja Ballooning Publication

Purpose
-------
Tests whether the latitudinal taxonomic pattern differs between genera coded as
ballooning and non-ballooning. In every Monte Carlo iteration, the script draws
the same number of occupied 25-km cells without replacement from each latitude
band. It then calculates:

  1. standardized genus richness separately for ballooning and non-ballooning
     genera in every latitude band;
  2. trait-specific pairwise beta diversity (Jaccard, Sørensen, Simpson turnover,
     and Sørensen nestedness-resultant) among latitude bands;
  3. paired richness contrasts among bands within each trait class; and
  4. paired contrasts between ballooning and non-ballooning turnover for each
     band transition.

The same cell draw is used for both trait classes within an iteration, so class
and band contrasts are paired by construction.

Inputs are discovered from the standard Baja_Ballooning_Pipeline folders. No
third-party Python packages are required.

Default run
-----------
python3 11G_equal_cell_trait_partitioned_richness_turnover.py

Optional arguments
------------------
python3 11G_equal_cell_trait_partitioned_richness_turnover.py \
  ~/Desktop/Baja_Ballooning_Pipeline 2000 20260713
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import random
import re
import shutil
import statistics
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence


SCRIPT_VERSION = "11G_C3_publication_v2_2026-07-16"
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

TRAIT_CLASSES = ["ballooning", "non_ballooning"]
TRAIT_LABELS = {
    "ballooning": "Ballooning",
    "non_ballooning": "Non-ballooning",
}

ANALYSIS_ORDER = [
    "primary",
    "taxonomy_strict",
    "explicit_low_confidence_exclusion",
]

ANALYSIS_LABELS = {
    "primary": "Primary",
    "taxonomy_strict": "Taxonomy-strict sensitivity",
    "explicit_low_confidence_exclusion": "Explicit LOW-confidence exclusion",
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
    paths: Sequence[Path],
    label: str,
    required: bool = True,
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
        return list(reader.fieldnames), [dict(row) for row in reader]


def write_csv(
    path: Path,
    rows: Iterable[dict[str, Any]],
    fieldnames: Sequence[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
    temporary.replace(path)


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
    cell_position = {cell: index for index, cell in enumerate(matrix["cells"])}
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
    cell_field = find_field(fields, ["grid_cell_id", "cell_id"], "grid-cell ID")
    band_field = find_field(
        fields,
        ["centroid_latitude_band", "latitude_band"],
        "latitude band",
    )

    mapping: dict[str, str] = {}
    for row in rows:
        cell = row.get(cell_field, "").strip()
        band = row.get(band_field, "").strip()
        if cell:
            mapping[cell] = band

    missing = [cell for cell in required_cells if cell not in mapping]
    if missing:
        raise ValueError(
            f"Cell lookup lacks {len(missing)} matrix cells. Examples: {missing[:10]}"
        )

    invalid = sorted(
        {mapping[cell] for cell in required_cells if mapping[cell] not in BAND_ORDER}
    )
    if invalid:
        raise ValueError(f"Cells assigned outside the five bands: {invalid}")

    return {cell: mapping[cell] for cell in required_cells}


EVIDENCE_TOKEN_RE = re.compile(
    r"(?<![A-Z0-9])(D1|D2|D3|D4|N0|C3)(?![A-Z0-9])", re.IGNORECASE
)


def parse_evidence_class(value: Any) -> str | None:
    """Parse an explicit D1/D2/D3/D4/N0/C3 trait designation.

    C3 is accepted as an already-derived primary class. D4 is always retained
    as an excluded class and is never recoded as non-ballooning.
    """
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
    if normalized in {
        "nonballooning", "fixednonballooning", "referencenonballooning",
        "noballooning", "nonballooningreference",
    }:
        return "N0"
    if normalized in {"c3", "primaryc3", "d1d2d3", "d1tod3"}:
        return "C3"
    if normalized in {"d4excluded", "excludedd4"}:
        return "D4"
    return None


def choose_evidence_field(fields: Sequence[str], rows: Sequence[dict[str, str]]) -> str:
    preferred = [
        "evidence_class", "final_evidence_class", "final_evidence_category",
        "evidence_category", "d_level", "dlevel", "trait_class",
        "primary_class", "analysis_class", "ballooning_evidence_tier",
        "ballooning_evidence_category", "final_designation", "designation",
    ]
    lower = {field.casefold(): field for field in fields}
    candidates: list[tuple[float, str]] = []
    for field in fields:
        parsed = [parse_evidence_class(row.get(field, "")) for row in rows]
        nonblank = [value for value in parsed if value is not None]
        if not nonblank:
            continue
        fraction = len(nonblank) / max(1, len(rows))
        classes = set(nonblank)
        name = re.sub(r"[^a-z0-9]+", "", field.casefold())
        bonus = 0.0
        if field.casefold() in {item.casefold() for item in preferred}:
            bonus += 100.0
        if any(term in name for term in ("evidence", "tier", "class", "designation", "decision")):
            bonus += 20.0
        if len(classes) >= 2 and fraction >= 0.25:
            candidates.append((bonus + 100.0 * fraction + 5.0 * len(classes), field))
    if not candidates:
        raise ValueError(
            "Trait table lacks an explicit D1/D2/D3/D4/N0 or C3/N0 evidence field. "
            "Legacy binary ballooning columns are intentionally not accepted in the "
            "publication workflow because they cannot distinguish D4 from N0."
        )
    candidates.sort(reverse=True)
    return candidates[0][1]


def load_traits(path: Path, genera: Sequence[str]) -> dict[str, Any]:
    fields, rows = read_csv_rows(path)
    genus_field = find_field(fields, ["genus", "analysis_genus"], "trait genus")
    evidence_field = choose_evidence_field(fields, rows)
    confidence_field = find_field(
        fields,
        [
            "final_confidence", "trait_final_confidence", "trait_confidence",
            "trait_ballooning_confidence",
        ],
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
        evidence = parse_evidence_class(row.get(evidence_field, ""))
        confidence = (
            row.get(confidence_field, "").strip().upper()
            if confidence_field is not None else "UNSPECIFIED"
        ) or "UNSPECIFIED"
        lookup[key] = {"genus": genus, "evidence": evidence, "confidence": confidence}

    missing: list[str] = []
    unresolved: list[str] = []
    c3_mask = 0
    n0_mask = 0
    d4_mask = 0
    low_confidence_mask = 0
    evidence_counts: dict[str, int] = defaultdict(int)
    confidence_counts: dict[str, int] = defaultdict(int)
    normalized_rows: list[dict[str, Any]] = []

    for genus_index, genus in enumerate(genera):
        key = genus.casefold()
        if key not in lookup:
            missing.append(genus)
            continue
        trait = lookup[key]
        evidence = trait["evidence"]
        if evidence is None:
            unresolved.append(genus)
            continue
        analysis_class = (
            "C3" if evidence in {"D1", "D2", "D3", "C3"}
            else "N0" if evidence == "N0"
            else "D4_excluded"
        )
        bit = 1 << genus_index
        if analysis_class == "C3":
            c3_mask |= bit
        elif analysis_class == "N0":
            n0_mask |= bit
        else:
            d4_mask |= bit
        if trait["confidence"] == "LOW":
            low_confidence_mask |= bit
        evidence_counts[evidence] += 1
        confidence_counts[trait["confidence"]] += 1
        normalized_rows.append({
            "genus": genus,
            "evidence_class": evidence,
            "analysis_class": analysis_class,
            "trait_class": (
                "ballooning" if analysis_class == "C3"
                else "non_ballooning" if analysis_class == "N0"
                else "excluded_D4"
            ),
            "final_confidence": trait["confidence"],
        })

    if missing:
        raise ValueError("Genera missing from trait table: " + ", ".join(missing[:20]))
    if unresolved:
        raise ValueError(
            "Genera lacking an explicit D1/D2/D3/D4/N0 or C3/N0 designation: "
            + ", ".join(unresolved[:20])
        )

    all_mask = (1 << len(genera)) - 1
    classified_mask = c3_mask | n0_mask
    confidence_keep_mask = all_mask & ~low_confidence_mask
    if c3_mask & n0_mask or c3_mask & d4_mask or n0_mask & d4_mask:
        raise ValueError("Trait classes are not disjoint.")
    if (classified_mask | d4_mask) != all_mask:
        raise ValueError("Trait classes do not cover every incidence-matrix genus.")

    return {
        "path": path,
        "evidence_field": evidence_field,
        "ballooning_mask": c3_mask,
        "non_ballooning_mask": n0_mask,
        "d4_excluded_mask": d4_mask,
        "classified_mask": classified_mask,
        "low_confidence_mask": low_confidence_mask,
        "confidence_keep_mask": confidence_keep_mask,
        "evidence_counts": dict(sorted(evidence_counts.items())),
        "confidence_counts": dict(sorted(confidence_counts.items())),
        "normalized_rows": normalized_rows,
    }

def union_masks(cell_masks: Sequence[int], selected_indices: Sequence[int]) -> int:
    result = 0
    for index in selected_indices:
        result |= cell_masks[index]
    return result


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
        "standard_deviation": statistics.stdev(clean) if len(clean) > 1 else 0.0,
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
    grouped: dict[tuple[Any, ...], dict[str, list[float]]] = defaultdict(
        lambda: defaultdict(list)
    )
    group_examples: dict[tuple[Any, ...], dict[str, Any]] = {}

    for row in rows:
        key = tuple(row[field] for field in group_fields)
        group_examples[key] = {field: row[field] for field in group_fields}
        for metric in metric_fields:
            grouped[key][metric].append(float(row[metric]))

    output: list[dict[str, Any]] = []
    for key, metrics in grouped.items():
        base = group_examples[key]
        for metric in metric_fields:
            output.append(
                {
                    **base,
                    "metric": metric,
                    **summarize_values(metrics[metric]),
                }
            )
    return output


def archive_existing(output_dir: Path, archive_root: Path) -> Path | None:
    if not output_dir.exists() or not any(output_dir.iterdir()):
        return None
    timestamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    archive_dir = archive_root / f"11G_trait_partitioned_{timestamp}"
    archive_dir.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(output_dir, archive_dir)
    shutil.rmtree(output_dir)
    return archive_dir


def format_number(value: float, digits: int = 3) -> str:
    if value is None or math.isnan(float(value)):
        return "NA"
    return f"{float(value):.{digits}f}"


def svg_text(
    x: float,
    y: float,
    text: str,
    size: int = 18,
    anchor: str = "middle",
    weight: str = "normal",
    rotate: float | None = None,
) -> str:
    transform = f' transform="rotate({rotate} {x} {y})"' if rotate is not None else ""
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" '
        f'font-family="Arial, Helvetica, sans-serif" font-size="{size}" '
        f'font-weight="{weight}" fill="#111"{transform}>'
        f'{html.escape(text)}</text>'
    )


def make_combined_primary_svg(
    path: Path,
    richness_summary: Sequence[dict[str, Any]],
    beta_summary: Sequence[dict[str, Any]],
    sample_size: int,
    iterations: int,
) -> None:
    width, height = 1600, 1120
    margin_left, margin_right = 145, 65
    panel_width = (width - margin_left - margin_right - 115) / 2
    panel_top, panel_bottom = 160, 890
    panel_height = panel_bottom - panel_top
    panel_a_left = margin_left
    panel_b_left = margin_left + panel_width + 115

    trait_styles = {
        "ballooning": {"color": "#0072B2", "offset": -13, "dash": ""},
        "non_ballooning": {"color": "#D55E00", "offset": 13, "dash": "7,5"},
    }

    rich_rows = [
        row
        for row in richness_summary
        if row["analysis"] == "primary"
        and row["metric"] == "trait_genus_richness"
    ]
    beta_rows = [
        row
        for row in beta_summary
        if row["analysis"] == "primary"
        and row["metric"] == "jaccard_dissimilarity"
        and str(row["adjacent_bands"]) == "True"
    ]

    if not rich_rows or not beta_rows:
        raise RuntimeError("Primary summary rows needed for the combined SVG are missing.")

    richness_lookup = {
        (row["latitude_band"], row["trait_class"]): row for row in rich_rows
    }
    beta_lookup = {
        (row["band_1"], row["band_2"], row["trait_class"]): row
        for row in beta_rows
    }

    richness_ymax = max(float(row["p975"]) for row in rich_rows)
    richness_ymax = max(10.0, math.ceil(richness_ymax / 10.0) * 10.0)

    def x_for(index: int, left: float, categories: int) -> float:
        if categories == 1:
            return left + panel_width / 2
        return left + index * panel_width / (categories - 1)

    def y_rich(value: float) -> float:
        return panel_bottom - (value / richness_ymax) * panel_height

    def y_beta(value: float) -> float:
        return panel_bottom - value * panel_height

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        svg_text(
            width / 2,
            58,
            "Equal-cell, trait-partitioned arachnid taxonomy across latitude",
            size=30,
            weight="bold",
        ),
        svg_text(
            width / 2,
            96,
            f"Means and 2.5th–97.5th percentile intervals from {iterations:,} paired resampling iterations; {sample_size} occupied 25-km cells per band",
            size=17,
        ),
    ]

    # Panel A axes and labels.
    parts.extend(
        [
            svg_text(panel_a_left - 85, panel_top - 35, "A", size=27, weight="bold"),
            svg_text(
                panel_a_left + panel_width / 2,
                panel_top - 35,
                "Trait-partitioned standardized genus richness",
                size=22,
                weight="bold",
            ),
            f'<line x1="{panel_a_left}" y1="{panel_bottom}" x2="{panel_a_left + panel_width}" y2="{panel_bottom}" stroke="#111" stroke-width="2"/>',
            f'<line x1="{panel_a_left}" y1="{panel_top}" x2="{panel_a_left}" y2="{panel_bottom}" stroke="#111" stroke-width="2"/>',
            svg_text(
                panel_a_left - 95,
                panel_top + panel_height / 2,
                f"Genus richness in {sample_size} occupied cells",
                size=19,
                rotate=-90,
            ),
        ]
    )

    rich_tick_step = 10 if richness_ymax <= 100 else 20
    for tick in range(0, int(richness_ymax) + 1, rich_tick_step):
        y = y_rich(float(tick))
        parts.append(
            f'<line x1="{panel_a_left}" y1="{y:.1f}" x2="{panel_a_left + panel_width}" y2="{y:.1f}" stroke="#ddd" stroke-width="1"/>'
        )
        parts.append(svg_text(panel_a_left - 15, y + 6, str(tick), size=16, anchor="end"))

    for i, band in enumerate(BAND_ORDER):
        x = x_for(i, panel_a_left, len(BAND_ORDER))
        parts.append(svg_text(x, panel_bottom + 38, BAND_LABELS[band], size=16))

    for trait_class in TRAIT_CLASSES:
        style = trait_styles[trait_class]
        points: list[tuple[float, float]] = []
        for i, band in enumerate(BAND_ORDER):
            row = richness_lookup[(band, trait_class)]
            x = x_for(i, panel_a_left, len(BAND_ORDER)) + style["offset"]
            y = y_rich(float(row["mean"]))
            y_low = y_rich(float(row["p025"]))
            y_high = y_rich(float(row["p975"]))
            points.append((x, y))
            parts.extend(
                [
                    f'<line x1="{x:.1f}" y1="{y_low:.1f}" x2="{x:.1f}" y2="{y_high:.1f}" stroke="{style["color"]}" stroke-width="3"/>',
                    f'<line x1="{x-7:.1f}" y1="{y_low:.1f}" x2="{x+7:.1f}" y2="{y_low:.1f}" stroke="{style["color"]}" stroke-width="3"/>',
                    f'<line x1="{x-7:.1f}" y1="{y_high:.1f}" x2="{x+7:.1f}" y2="{y_high:.1f}" stroke="{style["color"]}" stroke-width="3"/>',
                    f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" fill="white" stroke="{style["color"]}" stroke-width="4"/>',
                ]
            )
        point_string = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
        dash = f' stroke-dasharray="{style["dash"]}"' if style["dash"] else ""
        parts.append(
            f'<polyline points="{point_string}" fill="none" stroke="{style["color"]}" stroke-width="3"{dash}/>'
        )

    # Panel B axes and labels.
    adjacent_pairs = [pair for pair in PAIR_ORDER if pair in ADJACENT_PAIRS]
    parts.extend(
        [
            svg_text(panel_b_left - 85, panel_top - 35, "B", size=27, weight="bold"),
            svg_text(
                panel_b_left + panel_width / 2,
                panel_top - 35,
                "Trait-partitioned adjacent-band turnover",
                size=22,
                weight="bold",
            ),
            f'<line x1="{panel_b_left}" y1="{panel_bottom}" x2="{panel_b_left + panel_width}" y2="{panel_bottom}" stroke="#111" stroke-width="2"/>',
            f'<line x1="{panel_b_left}" y1="{panel_top}" x2="{panel_b_left}" y2="{panel_bottom}" stroke="#111" stroke-width="2"/>',
            svg_text(
                panel_b_left - 95,
                panel_top + panel_height / 2,
                "Jaccard dissimilarity",
                size=19,
                rotate=-90,
            ),
        ]
    )

    for tick_index in range(0, 11, 2):
        tick = tick_index / 10
        y = y_beta(tick)
        parts.append(
            f'<line x1="{panel_b_left}" y1="{y:.1f}" x2="{panel_b_left + panel_width}" y2="{y:.1f}" stroke="#ddd" stroke-width="1"/>'
        )
        parts.append(svg_text(panel_b_left - 15, y + 6, f"{tick:.1f}", size=16, anchor="end"))

    for i, pair in enumerate(adjacent_pairs):
        x = x_for(i, panel_b_left, len(adjacent_pairs))
        label = f"{BAND_LABELS[pair[0]]}\nto {BAND_LABELS[pair[1]]}"
        first, second = label.split("\n")
        parts.append(svg_text(x, panel_bottom + 34, first, size=15))
        parts.append(svg_text(x, panel_bottom + 56, second, size=15))

    for trait_class in TRAIT_CLASSES:
        style = trait_styles[trait_class]
        points = []
        for i, pair in enumerate(adjacent_pairs):
            row = beta_lookup[(pair[0], pair[1], trait_class)]
            x = x_for(i, panel_b_left, len(adjacent_pairs)) + style["offset"]
            y = y_beta(float(row["mean"]))
            y_low = y_beta(float(row["p025"]))
            y_high = y_beta(float(row["p975"]))
            points.append((x, y))
            parts.extend(
                [
                    f'<line x1="{x:.1f}" y1="{y_low:.1f}" x2="{x:.1f}" y2="{y_high:.1f}" stroke="{style["color"]}" stroke-width="3"/>',
                    f'<line x1="{x-7:.1f}" y1="{y_low:.1f}" x2="{x+7:.1f}" y2="{y_low:.1f}" stroke="{style["color"]}" stroke-width="3"/>',
                    f'<line x1="{x-7:.1f}" y1="{y_high:.1f}" x2="{x+7:.1f}" y2="{y_high:.1f}" stroke="{style["color"]}" stroke-width="3"/>',
                    f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" fill="white" stroke="{style["color"]}" stroke-width="4"/>',
                ]
            )
        point_string = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
        dash = f' stroke-dasharray="{style["dash"]}"' if style["dash"] else ""
        parts.append(
            f'<polyline points="{point_string}" fill="none" stroke="{style["color"]}" stroke-width="3"{dash}/>'
        )

    # Legend and explanatory footer.
    legend_y = 982
    for index, trait_class in enumerate(TRAIT_CLASSES):
        style = trait_styles[trait_class]
        legend_x = width / 2 - 245 + index * 350
        dash = f' stroke-dasharray="{style["dash"]}"' if style["dash"] else ""
        parts.append(
            f'<line x1="{legend_x:.1f}" y1="{legend_y}" x2="{legend_x + 60:.1f}" y2="{legend_y}" stroke="{style["color"]}" stroke-width="4"{dash}/>'
        )
        parts.append(
            f'<circle cx="{legend_x + 30:.1f}" cy="{legend_y}" r="7" fill="white" stroke="{style["color"]}" stroke-width="4"/>'
        )
        parts.append(
            svg_text(
                legend_x + 78,
                legend_y + 7,
                TRAIT_LABELS[trait_class],
                size=18,
                anchor="start",
            )
        )

    parts.append(
        svg_text(
            width / 2,
            1055,
            "Panel A asks whether either dispersal class accumulates taxonomically within particular bands. Panel B asks whether one class contributes more strongly to adjacent-band genus replacement.",
            size=16,
        )
    )
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def make_main_values_table(
    richness_summary: Sequence[dict[str, Any]],
    beta_summary: Sequence[dict[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in richness_summary:
        if row["analysis"] == "primary" and row["metric"] == "trait_genus_richness":
            rows.append(
                {
                    "panel": "A_richness",
                    "analysis": row["analysis"],
                    "category": row["latitude_band_label"],
                    "trait_class": row["trait_class"],
                    "metric": row["metric"],
                    "mean": row["mean"],
                    "median": row["median"],
                    "p025": row["p025"],
                    "p975": row["p975"],
                }
            )
    for row in beta_summary:
        if (
            row["analysis"] == "primary"
            and row["metric"] == "jaccard_dissimilarity"
            and str(row["adjacent_bands"]) == "True"
        ):
            rows.append(
                {
                    "panel": "B_turnover",
                    "analysis": row["analysis"],
                    "category": f"{row['band_1_label']} to {row['band_2_label']}",
                    "trait_class": row["trait_class"],
                    "metric": row["metric"],
                    "mean": row["mean"],
                    "median": row["median"],
                    "p025": row["p025"],
                    "p975": row["p975"],
                }
            )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Equal-cell, trait-partitioned genus richness and beta-diversity "
            "analysis for the Baja Ballooning project."
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
        help=f"Monte Carlo iterations; default {DEFAULT_ITERATIONS}.",
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
    trait_fallback = project_root / "02_data_clean" / "07_final_trait_merge"
    output_dir = project_root / "04_analysis" / "11G_trait_partitioned_equal_cell"
    figure_dir = output_dir / "figures"
    archive_root = project_root / "08_archive"

    archived = archive_existing(output_dir, archive_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    log_path = output_dir / "11G_analysis_log.txt"
    with log_path.open("w", encoding="utf-8") as log_handle:
        original_stdout = sys.stdout
        original_stderr = sys.stderr
        sys.stdout = Tee(original_stdout, log_handle)
        sys.stderr = Tee(original_stderr, log_handle)

        try:
            print("STEP 11G STARTED")
            print(f"Version: {SCRIPT_VERSION}")
            print(f"Started UTC: {utc_now()}")
            print(f"Project: {project_root}")
            print(f"Iterations: {args.iterations:,}")
            print(f"Random seed: {args.seed}")
            if archived is not None:
                print(f"Archived prior Step 11G output: {archived}")

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
                "taxonomy-strict biodiversity incidence matrix",
                required=False,
            )
            cell_lookup_path = first_existing(
                [
                    analysis_ready
                    / "04_spatial_reference"
                    / "10_common_grid25km_cell_lookup.csv",
                    grid_fallback / "10_common_grid25km_cell_lookup.csv",
                ],
                "25-km grid-cell lookup",
            )
            trait_path = first_existing(
                [
                    analysis_ready
                    / "03_trait_tables"
                    / "07_reviewed_genus_trait_lookup_normalized.csv",
                    analysis_ready
                    / "03_trait_tables"
                    / "07_reviewed_genus_trait_lookup_final.csv",
                    trait_fallback / "07_reviewed_genus_trait_lookup_final.csv",
                ],
                "final reviewed genus trait lookup",
            )

            print(f"Primary matrix: {primary_path}")
            print(f"Taxonomy-strict matrix: {strict_path or 'not available'}")
            print(f"Cell lookup: {cell_lookup_path}")
            print(f"Trait lookup: {trait_path}")

            primary = read_incidence_matrix(primary_path)
            genera = primary["genera"]
            cells = primary["cells"]
            if "fesa" in {genus.casefold() for genus in genera}:
                raise RuntimeError("Fesa remains in the primary incidence matrix.")

            primary_cell_masks = align_matrix(primary, genera, cells)
            matrix_masks: dict[str, list[int]] = {"primary": primary_cell_masks}

            if strict_path is not None:
                strict = read_incidence_matrix(strict_path)
                matrix_masks["taxonomy_strict"] = align_matrix(strict, genera, cells)

            cell_bands = load_cell_bands(cell_lookup_path, cells)
            band_cell_indices: dict[str, list[int]] = {
                band: [
                    index
                    for index, cell in enumerate(cells)
                    if cell_bands[cell] == band
                ]
                for band in BAND_ORDER
            }
            cell_counts = {band: len(indices) for band, indices in band_cell_indices.items()}
            if any(count == 0 for count in cell_counts.values()):
                raise RuntimeError(f"At least one latitude band has no occupied cells: {cell_counts}")
            sample_size = min(cell_counts.values())
            print(f"Occupied 25-km cells by band: {cell_counts}")
            print(f"Equal cells sampled per band: {sample_size}")

            traits = load_traits(trait_path, genera)
            all_mask = (1 << len(genera)) - 1
            trait_masks = {
                "ballooning": traits["ballooning_mask"],
                "non_ballooning": traits["non_ballooning_mask"],
            }
            keep_masks = {
                "primary": all_mask,
                "taxonomy_strict": all_mask,
                "explicit_low_confidence_exclusion": traits[
                    "confidence_keep_mask"
                ],
            }

            # Low-confidence sensitivity uses the primary occurrence matrix.
            matrix_masks["explicit_low_confidence_exclusion"] = primary_cell_masks
            analyses = [
                analysis
                for analysis in ANALYSIS_ORDER
                if analysis in matrix_masks
            ]

            regional_pools: dict[tuple[str, str], int] = {}
            for analysis in analyses:
                regional_mask = union_masks(
                    matrix_masks[analysis], range(len(matrix_masks[analysis]))
                ) & keep_masks[analysis]
                for trait_class in TRAIT_CLASSES:
                    regional_pools[(analysis, trait_class)] = (
                        regional_mask & trait_masks[trait_class]
                    ).bit_count()

            print(f"Analyses: {', '.join(analyses)}")
            print(
                "Primary regional trait pools: "
                f"ballooning={regional_pools[('primary', 'ballooning')]}, "
                f"non-ballooning={regional_pools[('primary', 'non_ballooning')]}"
            )
            print(f"Trait evidence counts: {traits['evidence_counts']}")
            print(f"Trait confidence counts: {traits['confidence_counts']}")

            validation_rows = [
                {
                    "check": "primary_matrix_has_no_Fesa",
                    "passed": True,
                    "detail": "Fesa absent",
                },
                {
                    "check": "all_five_bands_have_cells",
                    "passed": all(count > 0 for count in cell_counts.values()),
                    "detail": json.dumps(cell_counts, sort_keys=True),
                },
                {
                    "check": "equal_sample_size_at_least_1",
                    "passed": sample_size >= 1,
                    "detail": str(sample_size),
                },
                {
                    "check": "trait_classes_cover_all_genera",
                    "passed": (traits["classified_mask"] | traits["d4_excluded_mask"]) == all_mask,
                    "detail": (f"classified={traits['classified_mask'].bit_count()}; "
                               f"D4_excluded={traits['d4_excluded_mask'].bit_count()}; total={len(genera)}"),
                },
                {
                    "check": "primary_trait_classes_are_disjoint",
                    "passed": (traits["ballooning_mask"] & traits["non_ballooning_mask"]) == 0,
                    "detail": "No genus assigned to both C3 and N0",
                },
                {
                    "check": "primary_trait_pools_nonempty",
                    "passed": all(
                        regional_pools[("primary", trait_class)] > 0
                        for trait_class in TRAIT_CLASSES
                    ),
                    "detail": json.dumps(
                        {
                            trait_class: regional_pools[("primary", trait_class)]
                            for trait_class in TRAIT_CLASSES
                        },
                        sort_keys=True,
                    ),
                },
            ]
            if not all(bool(row["passed"]) for row in validation_rows):
                write_csv(
                    output_dir / "11G_validation.csv",
                    validation_rows,
                    ["check", "passed", "detail"],
                )
                failed = [row["check"] for row in validation_rows if not row["passed"]]
                raise RuntimeError("Step 11G validation failed: " + ", ".join(failed))

            write_csv(
                output_dir / "11G_validation.csv",
                validation_rows,
                ["check", "passed", "detail"],
            )
            write_csv(
                output_dir / "11G_normalized_trait_lookup.csv",
                traits["normalized_rows"],
                ["genus", "evidence_class", "analysis_class", "trait_class", "final_confidence"],
            )

            input_paths = [primary_path, cell_lookup_path, trait_path]
            if strict_path is not None:
                input_paths.append(strict_path)
            input_manifest = [
                {
                    "path": str(path),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for path in input_paths
            ]
            write_csv(
                output_dir / "11G_input_file_manifest.csv",
                input_manifest,
                ["path", "bytes", "sha256"],
            )

            rng = random.Random(args.seed)
            richness_rows: list[dict[str, Any]] = []
            beta_rows: list[dict[str, Any]] = []

            for iteration in range(1, args.iterations + 1):
                selected_by_band = {
                    band: rng.sample(band_cell_indices[band], sample_size)
                    for band in BAND_ORDER
                }

                for analysis in analyses:
                    cell_masks = matrix_masks[analysis]
                    keep_mask = keep_masks[analysis]
                    band_masks = {
                        band: union_masks(cell_masks, selected_by_band[band]) & keep_mask
                        for band in BAND_ORDER
                    }

                    for band in BAND_ORDER:
                        total_mask = band_masks[band]
                        total_richness = total_mask.bit_count()
                        classified_richness = (total_mask & traits["classified_mask"]).bit_count()
                        for trait_class in TRAIT_CLASSES:
                            trait_mask = total_mask & trait_masks[trait_class]
                            trait_richness = trait_mask.bit_count()
                            regional_pool = regional_pools[(analysis, trait_class)]
                            richness_rows.append(
                                {
                                    "iteration": iteration,
                                    "analysis": analysis,
                                    "latitude_band": band,
                                    "latitude_band_label": BAND_LABELS[band],
                                    "sampled_cells_per_band": sample_size,
                                    "trait_class": trait_class,
                                    "trait_class_label": TRAIT_LABELS[trait_class],
                                    "trait_genus_richness": trait_richness,
                                    "total_genus_richness": total_richness,
                                    "classified_genus_richness": classified_richness,
                                    "trait_proportion": (
                                        trait_richness / classified_richness
                                        if classified_richness
                                        else math.nan
                                    ),
                                    "regional_trait_pool_richness": regional_pool,
                                    "sampled_fraction_of_regional_trait_pool": (
                                        trait_richness / regional_pool
                                        if regional_pool
                                        else math.nan
                                    ),
                                }
                            )

                    for band_1, band_2 in PAIR_ORDER:
                        for trait_class in TRAIT_CLASSES:
                            mask_1 = band_masks[band_1] & trait_masks[trait_class]
                            mask_2 = band_masks[band_2] & trait_masks[trait_class]
                            metrics = beta_metrics(mask_1, mask_2)
                            beta_rows.append(
                                {
                                    "iteration": iteration,
                                    "analysis": analysis,
                                    "trait_class": trait_class,
                                    "trait_class_label": TRAIT_LABELS[trait_class],
                                    "band_1": band_1,
                                    "band_2": band_2,
                                    "band_1_label": BAND_LABELS[band_1],
                                    "band_2_label": BAND_LABELS[band_2],
                                    "adjacent_bands": (band_1, band_2) in ADJACENT_PAIRS,
                                    "sampled_cells_per_band": sample_size,
                                    **metrics,
                                }
                            )

                if iteration % max(1, args.iterations // 10) == 0:
                    print(f"Completed {iteration:,}/{args.iterations:,} iterations")

            richness_fields = [
                "iteration",
                "analysis",
                "latitude_band",
                "latitude_band_label",
                "sampled_cells_per_band",
                "trait_class",
                "trait_class_label",
                "trait_genus_richness",
                "total_genus_richness",
                "classified_genus_richness",
                "trait_proportion",
                "regional_trait_pool_richness",
                "sampled_fraction_of_regional_trait_pool",
            ]
            beta_fields = [
                "iteration",
                "analysis",
                "trait_class",
                "trait_class_label",
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
                output_dir / "11G_trait_partitioned_richness_iterations.csv",
                richness_rows,
                richness_fields,
            )
            write_csv(
                output_dir / "11G_trait_partitioned_beta_iterations.csv",
                beta_rows,
                beta_fields,
            )

            richness_metric_fields = [
                "trait_genus_richness",
                "total_genus_richness",
                "classified_genus_richness",
                "trait_proportion",
                "sampled_fraction_of_regional_trait_pool",
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

            richness_summary = summarize_long(
                richness_rows,
                [
                    "analysis",
                    "latitude_band",
                    "latitude_band_label",
                    "trait_class",
                    "trait_class_label",
                    "regional_trait_pool_richness",
                ],
                richness_metric_fields,
            )
            beta_summary = summarize_long(
                beta_rows,
                [
                    "analysis",
                    "trait_class",
                    "trait_class_label",
                    "band_1",
                    "band_2",
                    "band_1_label",
                    "band_2_label",
                    "adjacent_bands",
                ],
                beta_metric_fields,
            )

            summary_stat_fields = [
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
            richness_summary_fields = [
                "analysis",
                "latitude_band",
                "latitude_band_label",
                "trait_class",
                "trait_class_label",
                "regional_trait_pool_richness",
                *summary_stat_fields,
            ]
            beta_summary_fields = [
                "analysis",
                "trait_class",
                "trait_class_label",
                "band_1",
                "band_2",
                "band_1_label",
                "band_2_label",
                "adjacent_bands",
                *summary_stat_fields,
            ]

            write_csv(
                output_dir / "11G_trait_partitioned_richness_summary.csv",
                richness_summary,
                richness_summary_fields,
            )
            write_csv(
                output_dir / "11G_trait_partitioned_beta_summary.csv",
                beta_summary,
                beta_summary_fields,
            )
            write_csv(
                output_dir / "11G_trait_partitioned_adjacent_beta_summary.csv",
                [row for row in beta_summary if str(row["adjacent_bands"]) == "True"],
                beta_summary_fields,
            )

            # Paired richness contrasts among latitude bands within each class.
            richness_by_key = {
                (
                    row["iteration"],
                    row["analysis"],
                    row["trait_class"],
                    row["latitude_band"],
                ): float(row["trait_genus_richness"])
                for row in richness_rows
            }
            richness_contrast_rows: list[dict[str, Any]] = []
            for analysis in analyses:
                for trait_class in TRAIT_CLASSES:
                    for band_1, band_2 in PAIR_ORDER:
                        values = [
                            richness_by_key[(iteration, analysis, trait_class, band_2)]
                            - richness_by_key[(iteration, analysis, trait_class, band_1)]
                            for iteration in range(1, args.iterations + 1)
                        ]
                        stats = summarize_values(values)
                        richness_contrast_rows.append(
                            {
                                "analysis": analysis,
                                "trait_class": trait_class,
                                "trait_class_label": TRAIT_LABELS[trait_class],
                                "band_1": band_1,
                                "band_2": band_2,
                                "band_1_label": BAND_LABELS[band_1],
                                "band_2_label": BAND_LABELS[band_2],
                                "adjacent_bands": (band_1, band_2) in ADJACENT_PAIRS,
                                "contrast_definition": "band_2 minus band_1 trait genus richness",
                                "probability_difference_gt_zero": sum(
                                    value > 0 for value in values
                                )
                                / len(values),
                                "probability_difference_lt_zero": sum(
                                    value < 0 for value in values
                                )
                                / len(values),
                                **stats,
                            }
                        )

            richness_contrast_fields = [
                "analysis",
                "trait_class",
                "trait_class_label",
                "band_1",
                "band_2",
                "band_1_label",
                "band_2_label",
                "adjacent_bands",
                "contrast_definition",
                "probability_difference_gt_zero",
                "probability_difference_lt_zero",
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
                output_dir / "11G_trait_richness_band_contrasts_summary.csv",
                richness_contrast_rows,
                richness_contrast_fields,
            )
            write_csv(
                output_dir / "11G_trait_richness_adjacent_band_contrasts_summary.csv",
                [
                    row
                    for row in richness_contrast_rows
                    if str(row["adjacent_bands"]) == "True"
                ],
                richness_contrast_fields,
            )

            # Paired turnover contrasts between trait classes for each band pair.
            beta_by_key = {
                (
                    row["iteration"],
                    row["analysis"],
                    row["band_1"],
                    row["band_2"],
                    row["trait_class"],
                ): row
                for row in beta_rows
            }
            turnover_contrast_rows: list[dict[str, Any]] = []
            for analysis in analyses:
                for band_1, band_2 in PAIR_ORDER:
                    for metric in ["jaccard_dissimilarity", "simpson_turnover"]:
                        values = [
                            float(
                                beta_by_key[
                                    (
                                        iteration,
                                        analysis,
                                        band_1,
                                        band_2,
                                        "ballooning",
                                    )
                                ][metric]
                            )
                            - float(
                                beta_by_key[
                                    (
                                        iteration,
                                        analysis,
                                        band_1,
                                        band_2,
                                        "non_ballooning",
                                    )
                                ][metric]
                            )
                            for iteration in range(1, args.iterations + 1)
                        ]
                        clean_values = [value for value in values if not math.isnan(value)]
                        stats = summarize_values(clean_values)
                        turnover_contrast_rows.append(
                            {
                                "analysis": analysis,
                                "band_1": band_1,
                                "band_2": band_2,
                                "band_1_label": BAND_LABELS[band_1],
                                "band_2_label": BAND_LABELS[band_2],
                                "adjacent_bands": (band_1, band_2) in ADJACENT_PAIRS,
                                "metric": metric,
                                "contrast_definition": "ballooning minus non-ballooning",
                                "probability_difference_gt_zero": (
                                    sum(value > 0 for value in clean_values)
                                    / len(clean_values)
                                    if clean_values
                                    else math.nan
                                ),
                                "probability_difference_lt_zero": (
                                    sum(value < 0 for value in clean_values)
                                    / len(clean_values)
                                    if clean_values
                                    else math.nan
                                ),
                                **stats,
                            }
                        )

            turnover_contrast_fields = [
                "analysis",
                "band_1",
                "band_2",
                "band_1_label",
                "band_2_label",
                "adjacent_bands",
                "metric",
                "contrast_definition",
                "probability_difference_gt_zero",
                "probability_difference_lt_zero",
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
                output_dir / "11G_trait_turnover_class_contrasts_summary.csv",
                turnover_contrast_rows,
                turnover_contrast_fields,
            )
            write_csv(
                output_dir / "11G_trait_turnover_adjacent_class_contrasts_summary.csv",
                [
                    row
                    for row in turnover_contrast_rows
                    if str(row["adjacent_bands"]) == "True"
                ],
                turnover_contrast_fields,
            )

            main_values = make_main_values_table(richness_summary, beta_summary)
            write_csv(
                output_dir / "11G_main_figure_values.csv",
                main_values,
                [
                    "panel",
                    "analysis",
                    "category",
                    "trait_class",
                    "metric",
                    "mean",
                    "median",
                    "p025",
                    "p975",
                ],
            )

            make_combined_primary_svg(
                figure_dir / "Figure_3_trait_partitioned_richness_turnover.svg",
                richness_summary,
                beta_summary,
                sample_size,
                args.iterations,
            )

            readme = f"""STEP 11G — Equal-cell, trait-partitioned richness and turnover

Question
--------
Does dispersal strategy leave a taxonomic signature across latitude, such that
ballooning and non-ballooning genera show different standardized richness and
turnover patterns among bands?

Design
------
- {args.iterations:,} Monte Carlo iterations.
- {sample_size} occupied 25-km cells sampled without replacement from each band.
- The same cell draw is used for ballooning and non-ballooning genera within an
  iteration, making trait-class comparisons paired.
- Primary results use the final biodiversity matrix and final reviewed binary
  trait table.
- Sensitivities: taxonomy-strict occurrence matrix when available, and explicit
  LOW-confidence trait exclusion.

Primary figure
--------------
figures/Figure_3_trait_partitioned_richness_turnover.svg

Panel A shows mean standardized genus richness and 2.5th–97.5th percentile
resampling intervals for each trait class in each band. It addresses whether
one trait class accumulates taxonomically in particular bands.

Panel B shows trait-specific adjacent-band Jaccard dissimilarity. It addresses
whether ballooning or non-ballooning genera contribute more strongly to genus
replacement across each boundary.

Important interpretation
------------------------
Primary classes are C3 = D1 + D2 + D3 and fixed N0; D4 is excluded.

This is a taxonomic trait-partitioning analysis. It tests whether the identities
of genera distributed across the peninsula are structured with respect to the
assigned ballooning trait. It does not demonstrate that an individual occurrence
or dispersal event resulted from ballooning.

Intervals
---------
All plotted intervals are the 2.5th and 97.5th percentiles of the Monte Carlo
resampling distribution, not confidence intervals from a fitted parametric model.

Key tables
----------
- 11G_trait_partitioned_richness_summary.csv
- 11G_trait_partitioned_adjacent_beta_summary.csv
- 11G_trait_richness_adjacent_band_contrasts_summary.csv
- 11G_trait_turnover_adjacent_class_contrasts_summary.csv
- 11G_main_figure_values.csv

Contrast signs
--------------
- Richness contrasts are band_2 minus band_1. Positive values mean higher
  standardized trait richness in the second/northern band.
- Turnover class contrasts are ballooning minus non-ballooning. Positive values
  mean higher dissimilarity/turnover among ballooning genera.
"""
            (output_dir / "README_11G.txt").write_text(readme, encoding="utf-8")

            provenance = {
                "script": SCRIPT_VERSION,
                "created_utc": utc_now(),
                "project_root": str(project_root),
                "iterations": args.iterations,
                "seed": args.seed,
                "equal_cells_per_band": sample_size,
                "available_cells_by_band": cell_counts,
                "band_order": BAND_ORDER,
                "analyses": {
                    analysis: ANALYSIS_LABELS[analysis] for analysis in analyses
                },
                "trait_classes": TRAIT_LABELS,
                "primary_definition": "C3 = D1 + D2 + D3 versus fixed N0; D4 excluded",
                "trait_evidence_field": traits["evidence_field"],
                "trait_evidence_counts": traits["evidence_counts"],
                "regional_trait_pools": {
                    f"{analysis}:{trait_class}": regional_pools[
                        (analysis, trait_class)
                    ]
                    for analysis in analyses
                    for trait_class in TRAIT_CLASSES
                },
                "trait_confidence_counts": traits["confidence_counts"],
                "formulas": {
                    "jaccard_dissimilarity": "(b+c)/(a+b+c)",
                    "sorensen_dissimilarity": "(b+c)/(2a+b+c)",
                    "simpson_turnover": "min(b,c)/(a+min(b,c))",
                    "sorensen_nestedness_resultant": (
                        "sorensen_dissimilarity - simpson_turnover"
                    ),
                    "richness_contrast": "band_2 richness - band_1 richness",
                    "turnover_class_contrast": (
                        "ballooning metric - non-ballooning metric"
                    ),
                },
                "interval_definition": (
                    "2.5th and 97.5th percentiles of the paired Monte Carlo "
                    "resampling distribution"
                ),
                "input_manifest": input_manifest,
                "validation_passed": True,
            }
            (output_dir / "11G_provenance.json").write_text(
                json.dumps(provenance, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )

            running_script = Path(__file__).resolve()
            shutil.copy2(
                running_script,
                output_dir / "11G_equal_cell_trait_partitioned_richness_turnover.py",
            )

            output_files = sorted(
                path
                for path in output_dir.rglob("*")
                if path.is_file() and path.name != "11G_output_file_manifest.csv"
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
                output_dir / "11G_output_file_manifest.csv",
                output_manifest,
                ["relative_path", "bytes", "sha256"],
            )

            print()
            print("=" * 78)
            print("STEP 11G COMPLETED SUCCESSFULLY")
            print("=" * 78)
            print(f"Iterations: {args.iterations:,}")
            print(f"Equal cells per latitude band: {sample_size}")
            print(f"Trait-richness iteration rows: {len(richness_rows):,}")
            print(f"Trait-beta iteration rows: {len(beta_rows):,}")
            print(f"Outputs: {output_dir}")
            print(
                "Main figure: "
                + str(
                    figure_dir
                    / "Figure_3_trait_partitioned_richness_turnover.svg"
                )
            )
            return 0

        finally:
            sys.stdout = original_stdout
            sys.stderr = original_stderr


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nStep 11G interrupted by user.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:
        print(f"\nSTEP 11G FAILED: {exc}", file=sys.stderr)
        raise
