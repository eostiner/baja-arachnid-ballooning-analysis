#!/usr/bin/env python3
"""
STEP 11E — CLEAN MANUSCRIPT JACCARD MAP
Baja Ballooning Publication

Creates a clean single-panel peninsula map with all five latitude bands and
the four calculated adjacent-band Jaccard dissimilarities from Step 11B.
Jaccard medians and 95% intervals are labeled directly at the four boundaries,
removing the duplicate side boxes and color bar from Step 11D.

Run:
    python3 11E_clean_jaccard_map.py \
        ~/Desktop/Baja_Ballooning_Pipeline

Outputs:
    <project>/04_analysis/11E_clean_jaccard_map/
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable
from urllib.request import Request, urlopen

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPolygon
    from matplotlib.ticker import FuncFormatter
except ImportError as exc:  # pragma: no cover - user environment check
    raise SystemExit(
        "Missing matplotlib. Install it with:\n"
        "  python3 -m pip install matplotlib\n"
        f"Original error: {exc}"
    )

try:
    from shapely.geometry import box, shape
    from shapely.ops import unary_union
except ImportError as exc:  # pragma: no cover - user environment check
    raise SystemExit(
        "Missing shapely. Install it with:\n"
        "  python3 -m pip install shapely\n"
        f"Original error: {exc}"
    )


SCRIPT_VERSION = "11E_v1_2026-07-13"
DEFAULT_PROJECT = Path.home() / "Desktop" / "Baja_Ballooning_Pipeline"

BAND_ORDER = ["23-24N", "24-26N", "26-28N", "28-30N", "30-32N"]
BAND_BOUNDS = {
    "23-24N": (23.0, 24.0),
    "24-26N": (24.0, 26.0),
    "26-28N": (26.0, 28.0),
    "28-30N": (28.0, 30.0),
    "30-32N": (30.0, 32.0),
}
BAND_LABELS = {
    "23-24N": "23–24°N",
    "24-26N": "24–26°N",
    "26-28N": "26–28°N",
    "28-30N": "28–30°N",
    "30-32N": "30–32°N",
}
PAIR_ORDER = [
    ("23-24N", "24-26N", 24.0),
    ("24-26N", "26-28N", 26.0),
    ("26-28N", "28-30N", 28.0),
    ("28-30N", "30-32N", 30.0),
]
ANALYSIS_LABELS = {
    "primary": "Primary analysis",
    "taxonomy_strict": "Taxonomy-strict sensitivity",
    "explicit_low_confidence_exclusion": "Explicit LOW-confidence exclusion",
}
TARGET_STATES = {"baja california", "baja california sur"}
GEOB_API = "https://www.geoboundaries.org/api/current/gbOpen/MEX/ADM1/"
USER_AGENT = "Baja-Ballooning-Pipeline/11E"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def normalize_text(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip().casefold())


def normalize_band(value: Any) -> str:
    text = str(value or "").strip()
    text = text.replace("–", "-").replace("—", "-").replace("°", "")
    text = re.sub(r"\s+", "", text).upper()
    if text.endswith("N"):
        text = text[:-1]
    return f"{text}N"


def parse_bool(value: Any) -> bool:
    return normalize_text(value) in {"true", "t", "1", "yes", "y"}


def parse_float(value: Any) -> float:
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def find_step11b_summary(step11b_dir: Path) -> Path:
    preferred = [
        step11b_dir / "11B_equal_cell_adjacent_band_turnover_summary.csv",
        step11b_dir / "11B_equal_cell_pairwise_beta_summary.csv",
    ]
    for path in preferred:
        if path.exists():
            return path

    required = {"analysis", "metric", "band_1", "band_2", "median"}
    for path in sorted(step11b_dir.rglob("*.csv")):
        try:
            with path.open("r", encoding="utf-8-sig", newline="") as source:
                header = next(csv.reader(source), [])
        except OSError:
            continue
        if required.issubset(set(header)):
            return path

    raise FileNotFoundError(
        "Could not locate a Step 11B beta-diversity summary CSV under:\n"
        f"  {step11b_dir}"
    )


def canonical_pair(first: str, second: str) -> tuple[str, str]:
    if first not in BAND_ORDER or second not in BAND_ORDER:
        raise ValueError(f"Unexpected latitude-band pair: {first!r}, {second!r}")
    if BAND_ORDER.index(first) < BAND_ORDER.index(second):
        return first, second
    return second, first


def load_jaccard_rows(summary_path: Path) -> dict[str, dict[tuple[str, str], dict[str, Any]]]:
    by_analysis: dict[str, dict[tuple[str, str], dict[str, Any]]] = {}

    with summary_path.open("r", encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        required = {"analysis", "metric", "band_1", "band_2", "median"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Summary file is missing required columns: {sorted(missing)}\n"
                f"File: {summary_path}"
            )

        for row in reader:
            if normalize_text(row.get("metric")) != "jaccard_dissimilarity":
                continue
            if "adjacent_bands" in row and not parse_bool(row.get("adjacent_bands")):
                continue

            analysis = normalize_text(row.get("analysis")).replace(" ", "_")
            first = normalize_band(row.get("band_1"))
            second = normalize_band(row.get("band_2"))
            pair = canonical_pair(first, second)

            record = {
                "analysis": analysis,
                "band_1": pair[0],
                "band_2": pair[1],
                "median": parse_float(row.get("median")),
                "p025": parse_float(row.get("p025")),
                "p975": parse_float(row.get("p975")),
                "mean": parse_float(row.get("mean")),
                "n_iterations": int(parse_float(row.get("n_iterations")))
                if math.isfinite(parse_float(row.get("n_iterations")))
                else None,
            }
            if not math.isfinite(record["median"]):
                continue
            by_analysis.setdefault(analysis, {})[pair] = record

    required_pairs = {(a, b) for a, b, _ in PAIR_ORDER}
    if "primary" not in by_analysis:
        raise ValueError("No primary Jaccard rows were found in the Step 11B summary.")
    missing_primary = required_pairs - set(by_analysis["primary"])
    if missing_primary:
        raise ValueError(
            "Primary analysis does not contain all four adjacent-band comparisons: "
            + ", ".join(f"{a} vs {b}" for a, b in sorted(missing_primary))
        )

    return by_analysis


def load_provenance(step11b_dir: Path, analyses: dict[str, Any]) -> tuple[int | None, int | None]:
    iterations: int | None = None
    cells_per_band: int | None = None
    provenance_path = step11b_dir / "11B_provenance.json"

    if provenance_path.exists():
        try:
            provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
            raw_iterations = provenance.get("iterations")
            if raw_iterations is not None:
                iterations = int(raw_iterations)
            raw_cells = provenance.get("sampling_design", {}).get("equal_cells_per_band")
            if raw_cells is not None:
                cells_per_band = int(raw_cells)
        except (OSError, ValueError, TypeError, json.JSONDecodeError):
            pass

    if iterations is None:
        values = [
            row.get("n_iterations")
            for pair_rows in analyses.values()
            for row in pair_rows.values()
            if row.get("n_iterations")
        ]
        if values:
            iterations = max(values)

    return iterations, cells_per_band


def open_url(url: str):
    return urlopen(Request(url, headers={"User-Agent": USER_AGENT}), timeout=180)


def locate_or_download_boundary(project_root: Path) -> Path:
    boundary_dir = project_root / "01_data_raw" / "boundaries"
    exact = boundary_dir / "geoBoundaries-MEX-ADM1.geojson"
    if exact.exists():
        return exact

    candidates = sorted(boundary_dir.glob("*.geojson")) if boundary_dir.exists() else []
    candidates += sorted(project_root.glob("**/*MEX*ADM1*.geojson"))
    for path in candidates:
        if path.is_file():
            return path

    boundary_dir.mkdir(parents=True, exist_ok=True)
    print("Boundary GeoJSON not found locally; downloading geoBoundaries MEX ADM1...")
    with open_url(GEOB_API) as response:
        metadata = json.load(response)
    download_url = metadata.get("gjDownloadURL")
    if not download_url:
        raise RuntimeError("geoBoundaries metadata did not include gjDownloadURL.")
    with open_url(download_url) as source, exact.open("wb") as destination:
        shutil.copyfileobj(source, destination)
    return exact


def load_baja_geometry(boundary_path: Path):
    geojson = json.loads(boundary_path.read_text(encoding="utf-8"))
    selected = []
    available = []

    for feature in geojson.get("features", []):
        properties = feature.get("properties", {})
        name = (
            properties.get("shapeName")
            or properties.get("name")
            or properties.get("NAME_1")
            or properties.get("admin1Name")
            or ""
        )
        available.append(name)
        if normalize_text(name) in TARGET_STATES:
            geom = shape(feature["geometry"])
            if not geom.is_valid:
                geom = geom.buffer(0)
            selected.append(geom)

    if len(selected) != 2:
        raise RuntimeError(
            "Could not identify both Baja California and Baja California Sur in:\n"
            f"  {boundary_path}\n"
            "Available unit names include: " + ", ".join(sorted(set(available))[:40])
        )

    geometry = unary_union(selected)
    if not geometry.is_valid:
        geometry = geometry.buffer(0)
    return geometry


def polygon_parts(geometry) -> Iterable[Any]:
    if geometry.is_empty:
        return []
    if geometry.geom_type == "Polygon":
        return [geometry]
    if geometry.geom_type == "MultiPolygon":
        return list(geometry.geoms)
    if geometry.geom_type == "GeometryCollection":
        parts = []
        for item in geometry.geoms:
            parts.extend(polygon_parts(item))
        return parts
    return []


def draw_geometry(ax, geometry, *, facecolor="none", edgecolor="black", linewidth=1.0, alpha=1.0, zorder=1):
    for polygon in polygon_parts(geometry):
        exterior = list(polygon.exterior.coords)
        patch = MplPolygon(
            exterior,
            closed=True,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
            alpha=alpha,
            zorder=zorder,
        )
        ax.add_patch(patch)

        # White-out holes so islands/coastline remain accurate.
        for interior in polygon.interiors:
            hole = MplPolygon(
                list(interior.coords),
                closed=True,
                facecolor="white",
                edgecolor="none",
                linewidth=0,
                zorder=zorder + 0.1,
            )
            ax.add_patch(hole)


def format_value(record: dict[str, Any], include_interval: bool = True) -> str:
    median = record["median"]
    p025 = record.get("p025", math.nan)
    p975 = record.get("p975", math.nan)
    if include_interval and math.isfinite(p025) and math.isfinite(p975):
        return f"{median:.3f} ({p025:.3f}–{p975:.3f})"
    return f"{median:.3f}"


def band_neighbor_records(pair_rows: dict[tuple[str, str], dict[str, Any]], band: str):
    index = BAND_ORDER.index(band)
    south = None
    north = None
    if index > 0:
        south_pair = (BAND_ORDER[index - 1], band)
        south = pair_rows[south_pair]
    if index < len(BAND_ORDER) - 1:
        north_pair = (band, BAND_ORDER[index + 1])
        north = pair_rows[north_pair]
    return south, north


def map_title(analysis: str) -> str:
    if analysis == "primary":
        return "Jaccard dissimilarity between adjacent Baja California latitude bands"
    return (
        "Jaccard dissimilarity between adjacent Baja California latitude bands\n"
        f"{ANALYSIS_LABELS.get(analysis, analysis)}"
    )


def create_map(
    geometry,
    pair_rows: dict[tuple[str, str], dict[str, Any]],
    analysis: str,
    figure_stem: Path,
    iterations: int | None,
    cells_per_band: int | None,
) -> list[dict[str, Any]]:
    """Create a clean single-panel map.

    Design choices:
    - five subtle alternating latitude-band fills;
    - one label per band;
    - Jaccard median and 95% interval placed directly on each shared boundary;
    - no duplicate right-hand boxes and no color bar;
    - explanatory details are moved to the caption rather than the plotting area.
    """
    fig, ax = plt.subplots(figsize=(8.4, 8.9))

    band_fills = ["#eef4f8", "#dbeaf3", "#eef4f8", "#dbeaf3", "#eef4f8"]
    boundary_color = "#1f5a85"

    # Draw five clipped latitude bands with restrained alternating fills.
    for band, fill in zip(BAND_ORDER, band_fills):
        lower, upper = BAND_BOUNDS[band]
        clipped = geometry.intersection(box(-120.5, lower, -107.5, upper))
        draw_geometry(
            ax,
            clipped,
            facecolor=fill,
            edgecolor="none",
            linewidth=0,
            alpha=1.0,
            zorder=1,
        )

    draw_geometry(
        ax,
        geometry,
        facecolor="none",
        edgecolor="#222222",
        linewidth=1.25,
        alpha=1.0,
        zorder=4,
    )

    # Label the latitude bands once, in open water east of the peninsula.
    for band in BAND_ORDER:
        lower, upper = BAND_BOUNDS[band]
        ax.text(
            -109.18,
            (lower + upper) / 2,
            BAND_LABELS[band],
            ha="left",
            va="center",
            fontsize=10.2,
            fontweight="bold",
            color="#202020",
            zorder=8,
        )

    # Four adjacent-band comparisons are shown directly at the four boundaries.
    # A white gap behind the text keeps labels legible without adding a second panel.
    for first, second, latitude in PAIR_ORDER:
        record = pair_rows[(first, second)]
        ax.axhline(latitude, color=boundary_color, linewidth=1.65, zorder=5)
        p025 = record.get("p025", math.nan)
        p975 = record.get("p975", math.nan)
        if math.isfinite(p025) and math.isfinite(p975):
            label = f"J = {record['median']:.3f}  [{p025:.3f}–{p975:.3f}]"
        else:
            label = f"J = {record['median']:.3f}"
        ax.text(
            -118.92,
            latitude + 0.09,
            label,
            ha="left",
            va="bottom",
            fontsize=9.1,
            color="#111111",
            bbox=dict(
                boxstyle="round,pad=0.18",
                facecolor="white",
                edgecolor="none",
                alpha=0.94,
            ),
            zorder=9,
        )

    ax.set_xlim(-119.30, -108.15)
    ax.set_ylim(22.70, 32.75)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([-118, -116, -114, -112, -110])
    ax.set_yticks([23, 24, 26, 28, 30, 32])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{abs(x):.0f}°W"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:.0f}°N"))
    ax.set_xlabel("Longitude", fontsize=11.2)
    ax.set_ylabel("Latitude", fontsize=11.2)
    ax.tick_params(labelsize=9.5)
    ax.grid(False)

    # Clean title hierarchy. The primary manuscript figure has no subtitle.
    title = "Adjacent-band Jaccard dissimilarity across the Baja California Peninsula"
    fig.suptitle(title, fontsize=14.5, y=0.965)
    if analysis != "primary":
        fig.text(
            0.50,
            0.925,
            ANALYSIS_LABELS.get(analysis, analysis),
            ha="center",
            va="center",
            fontsize=10.5,
            color="#444444",
        )

    fig.text(
        0.50,
        0.060,
        "Labels show median Jaccard dissimilarity [95% resampling interval].",
        ha="center",
        va="center",
        fontsize=9.2,
        color="#444444",
    )

    fig.subplots_adjust(left=0.12, right=0.94, bottom=0.11, top=0.88)

    export_status: list[dict[str, Any]] = []
    formats = [
        ("png", {"dpi": 600, "bbox_inches": "tight"}),
        ("pdf", {"bbox_inches": "tight"}),
        ("svg", {"bbox_inches": "tight"}),
        ("tiff", {"dpi": 600, "bbox_inches": "tight", "pil_kwargs": {"compression": "tiff_lzw"}}),
    ]

    for extension, kwargs in formats:
        output_path = figure_stem.with_suffix(f".{extension}")
        try:
            fig.savefig(output_path, **kwargs)
            success = output_path.exists() and output_path.stat().st_size > 0
            error = "" if success else "File was not created or is empty."
        except Exception as exc:
            success = False
            error = str(exc)
        export_status.append(
            {
                "figure": figure_stem.name,
                "analysis": analysis,
                "format": extension,
                "success": success,
                "path": str(output_path),
                "error": error,
            }
        )

    plt.close(fig)
    return export_status

def write_values_csv(path: Path, analyses: dict[str, dict[tuple[str, str], dict[str, Any]]]):
    fields = [
        "analysis",
        "band_1",
        "band_2",
        "boundary_latitude",
        "median_jaccard_dissimilarity",
        "p025",
        "p975",
        "n_iterations",
    ]
    with path.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=fields)
        writer.writeheader()
        for analysis in ANALYSIS_LABELS:
            if analysis not in analyses:
                continue
            for first, second, latitude in PAIR_ORDER:
                record = analyses[analysis].get((first, second))
                if record is None:
                    continue
                writer.writerow(
                    {
                        "analysis": analysis,
                        "band_1": first,
                        "band_2": second,
                        "boundary_latitude": latitude,
                        "median_jaccard_dissimilarity": record["median"],
                        "p025": record["p025"],
                        "p975": record["p975"],
                        "n_iterations": record["n_iterations"],
                    }
                )


def write_caption(
    path: Path,
    pair_rows: dict[tuple[str, str], dict[str, Any]],
    iterations: int | None,
    cells_per_band: int | None,
):
    pair_text = "; ".join(
        (
            f"{BAND_LABELS[first]}–{BAND_LABELS[second]}: "
            f"{format_value(pair_rows[(first, second)])}"
        )
        for first, second, _ in PAIR_ORDER
    )
    iteration_text = f"{iterations:,} Monte Carlo iterations" if iterations else "Monte Carlo resampling"
    cell_text = (
        f", each using {cells_per_band} occupied 25-km grid cells per latitude band"
        if cells_per_band
        else ""
    )
    caption = (
        "Figure 2. Equal-cell Jaccard dissimilarity between adjacent latitudinal bands "
        "of the Baja California Peninsula. Horizontal lines mark the four shared boundaries "
        "among the five latitude bands; labels report median Jaccard dissimilarity with "
        f"95% resampling intervals from {iteration_text}{cell_text}. "
        f"Adjacent-band estimates were {pair_text}. Higher values indicate greater "
        "compositional difference and stronger beta diversity between neighboring bands."
    )
    path.write_text(caption + "\n", encoding="utf-8")


def archive_existing(output_dir: Path, archive_root: Path):
    if not output_dir.exists() or not any(output_dir.iterdir()):
        return None
    timestamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    archive_dir = archive_root / f"11E_clean_jaccard_map_{timestamp}"
    archive_dir.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(output_dir, archive_dir)
    shutil.rmtree(output_dir)
    return archive_dir


def main() -> int:
    parser = argparse.ArgumentParser(description="Create clean manuscript Jaccard latitude-band map.")
    parser.add_argument(
        "project_root",
        nargs="?",
        default=str(DEFAULT_PROJECT),
        help="Baja Ballooning project root (default: ~/Desktop/Baja_Ballooning_Pipeline)",
    )
    args = parser.parse_args()

    project_root = Path(args.project_root).expanduser().resolve()
    if not project_root.exists():
        raise FileNotFoundError(f"Project folder not found: {project_root}")

    step11b_dir = project_root / "04_analysis" / "11B_equal_cell_resampling"
    if not step11b_dir.exists():
        raise FileNotFoundError(f"Step 11B output folder not found: {step11b_dir}")

    output_dir = project_root / "04_analysis" / "11E_clean_jaccard_map"
    figure_dir = output_dir / "figures"
    archive_root = project_root / "08_archive"
    archived = archive_existing(output_dir, archive_root)
    figure_dir.mkdir(parents=True, exist_ok=True)

    log_path = output_dir / "11E_analysis_log.txt"
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    with log_path.open("w", encoding="utf-8") as log:
        class Tee:
            def __init__(self, *streams):
                self.streams = streams

            def write(self, text):
                for stream in self.streams:
                    stream.write(text)
                    stream.flush()

            def flush(self):
                for stream in self.streams:
                    stream.flush()

        sys.stdout = Tee(original_stdout, log)
        sys.stderr = Tee(original_stderr, log)

        try:
            print("STEP 11E STARTED")
            print(f"Version: {SCRIPT_VERSION}")
            print(f"Project: {project_root}")
            if archived:
                print(f"Archived previous Step 11E output: {archived}")

            summary_path = find_step11b_summary(step11b_dir)
            print(f"Step 11B summary: {summary_path}")
            analyses = load_jaccard_rows(summary_path)
            iterations, cells_per_band = load_provenance(step11b_dir, analyses)
            print(f"Iterations: {iterations if iterations is not None else 'not found'}")
            print(f"Equal cells per band: {cells_per_band if cells_per_band is not None else 'not found'}")

            boundary_path = locate_or_download_boundary(project_root)
            print(f"Boundary: {boundary_path}")
            geometry = load_baja_geometry(boundary_path)

            export_rows: list[dict[str, Any]] = []
            required_pairs = {(a, b) for a, b, _ in PAIR_ORDER}
            for analysis in ANALYSIS_LABELS:
                pair_rows = analyses.get(analysis)
                if pair_rows is None:
                    print(f"Skipping unavailable analysis: {analysis}")
                    continue
                missing = required_pairs - set(pair_rows)
                if missing:
                    print(f"Skipping incomplete analysis {analysis}; missing {sorted(missing)}")
                    continue

                stem = figure_dir / f"Figure_2_Jaccard_map_clean_{analysis}"
                print(f"Creating: {stem.name}")
                export_rows.extend(
                    create_map(
                        geometry,
                        pair_rows,
                        analysis,
                        stem,
                        iterations,
                        cells_per_band,
                    )
                )

            write_values_csv(output_dir / "11E_adjacent_jaccard_map_values.csv", analyses)
            write_caption(
                output_dir / "Figure_2_caption.txt",
                analyses["primary"],
                iterations,
                cells_per_band,
            )

            status_path = output_dir / "11E_figure_export_status.csv"
            with status_path.open("w", encoding="utf-8", newline="") as destination:
                fields = ["figure", "analysis", "format", "success", "path", "error"]
                writer = csv.DictWriter(destination, fieldnames=fields)
                writer.writeheader()
                writer.writerows(export_rows)

            provenance = {
                "created_utc": utc_now(),
                "script_version": SCRIPT_VERSION,
                "script_path": str(Path(__file__).resolve()),
                "project_root": str(project_root),
                "step11b_summary": str(summary_path),
                "boundary_geojson": str(boundary_path),
                "iterations": iterations,
                "equal_cells_per_band": cells_per_band,
                "latitude_bands": BAND_ORDER,
                "adjacent_pairs": [[a, b] for a, b, _ in PAIR_ORDER],
                "interpretation": "Higher Jaccard dissimilarity means greater compositional difference.",
                "validation_passed": all(row["success"] for row in export_rows if row["format"] in {"png", "pdf"}),
            }
            (output_dir / "11E_provenance.json").write_text(
                json.dumps(provenance, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )

            shutil.copy2(Path(__file__).resolve(), output_dir / Path(__file__).name)

            manifest_rows = []
            for path in sorted(output_dir.rglob("*")):
                if path.is_file() and path.name != "11E_output_file_manifest.csv":
                    manifest_rows.append(
                        {
                            "relative_path": str(path.relative_to(output_dir)),
                            "bytes": path.stat().st_size,
                            "sha256": sha256(path),
                        }
                    )
            with (output_dir / "11E_output_file_manifest.csv").open(
                "w", encoding="utf-8", newline=""
            ) as destination:
                writer = csv.DictWriter(
                    destination, fieldnames=["relative_path", "bytes", "sha256"]
                )
                writer.writeheader()
                writer.writerows(manifest_rows)

            required_failures = [
                row for row in export_rows if row["format"] in {"png", "pdf"} and not row["success"]
            ]
            if required_failures:
                raise RuntimeError(
                    "One or more required PNG/PDF exports failed. See 11E_figure_export_status.csv"
                )

            print()
            print("=" * 78)
            print("STEP 11E COMPLETED SUCCESSFULLY")
            print("=" * 78)
            print("Primary manuscript map:")
            print(figure_dir / "Figure_2_Jaccard_map_clean_primary.png")
            print(f"Outputs: {output_dir}")
            return 0

        finally:
            sys.stdout = original_stdout
            sys.stderr = original_stderr


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nStep 11E interrupted by user.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:
        print(f"\nSTEP 11E FAILED: {exc}", file=sys.stderr)
        raise
