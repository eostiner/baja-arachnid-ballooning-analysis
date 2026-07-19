#!/usr/bin/env python3
"""Phase 15G — publication figure linking spatial turnover and temporal stress tracking.

The script uses only previously generated outputs. It does not refit any model.

Panel A
    Reconstructs adjacent-boundary and peninsula-wide C3 minus N0 turnover
    contrasts only from the frozen original 5,000-resample table
    (05_boundary_iteration_values.csv). The peninsula-wide value is
    computed within every Monte Carlo iteration as the mean of the four
    adjacent-boundary contrasts, then summarized across iterations.

Panel B
    Displays the 15 paired temporal observations and their equal-event
    resampling intervals against standardized recent-stress change. The fitted
    line and ribbon are calculated from the Phase 15 posterior draws while
    holding the regularizing covariates at their means and the cell effect at 0.

Panel C
    Reports the primary Bayesian coefficient using the current convention of a
    posterior point estimate with nested 50% and 95% credible intervals. Predictive
    validation, prior sensitivity, and posterior-density diagnostics remain in the
    supplementary Phase 15 outputs rather than the main synthesis figure.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

SCRIPT_VERSION = "15G_v0.5.0_2026-07-19"
ADJACENT_ORDER = [
    ("23-24N", "24-26N"),
    ("24-26N", "26-28N"),
    ("26-28N", "28-30N"),
    ("28-30N", "30-32N"),
]
BAND_LABELS = {
    "23-24N": "23–24°N",
    "24-26N": "24–26°N",
    "26-28N": "26–28°N",
    "28-30N": "28–30°N",
    "30-32N": "30–32°N",
}

# Frozen manuscript values from the original accepted 5,000-resample C3 run.
# The script aborts if a different Monte Carlo output is supplied.
FROZEN_SPATIAL_EXPECTED = {
    "simpson_turnover": {"median": -0.141745, "q025": -0.275743, "q975": -0.013406},
    "jaccard_dissimilarity": {"median": -0.007552, "q025": -0.084820, "q975": 0.072382},
}
FROZEN_CHECK_TOLERANCE = 1e-5
CANONICAL_SPATIAL_RELATIVE_PATHS = [
    "04_analysis/C3_pipeline_rebuild/04_C3_trait_turnover_5000/05_boundary_iteration_values.csv",
    "04_analysis_USE _THIS/C3_pipeline_rebuild/04_C3_trait_turnover_5000/05_boundary_iteration_values.csv",
    "C3_pipeline_rebuild/04_C3_trait_turnover_5000/05_boundary_iteration_values.csv",
]


SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
INKSCAPE_NS = "http://www.inkscape.org/namespaces/inkscape"
XML_NS = "http://www.w3.org/XML/1998/namespace"
ET.register_namespace("", SVG_NS)
ET.register_namespace("xlink", XLINK_NS)
ET.register_namespace("inkscape", INKSCAPE_NS)
_TRANSLATE_RE = re.compile(
    r"^\s*translate\(\s*([-+0-9.eE]+)[ ,]+([-+0-9.eE]+)\s*\)\s*$"
)


def _svg_text_content(element: ET.Element) -> str:
    """Return all textual content from an SVG text element."""
    return "".join(element.itertext())


def _svg_text_xy(element: ET.Element) -> tuple[str, str] | None:
    """Extract absolute x/y positions from Matplotlib SVG text output."""
    if "x" in element.attrib and "y" in element.attrib and not element.attrib.get("transform"):
        return element.attrib["x"], element.attrib["y"]
    match = _TRANSLATE_RE.match(element.attrib.get("transform", ""))
    if match:
        return match.group(1), match.group(2)
    return None


def make_svg_text_editable(svg_path: Path) -> dict[str, object]:
    """Preserve each Matplotlib text artist as one editable SVG text object.

    Matplotlib normally outlines fonts in SVG output. The rcParams below prevent
    outlining. This post-processor additionally collapses multiline labels emitted
    as sibling ``<text>`` elements into one ``<text>`` object containing ``<tspan>``
    lines. Thus a phrase or multiline label created by one Matplotlib text call is
    imported into Inkscape, Illustrator, or Affinity Designer as one grouped text
    object rather than as separate glyph paths or separate line objects.
    """
    tree = ET.parse(svg_path)
    root = tree.getroot()
    text_tag = f"{{{SVG_NS}}}text"
    tspan_tag = f"{{{SVG_NS}}}tspan"
    group_tag = f"{{{SVG_NS}}}g"
    use_tag = f"{{{SVG_NS}}}use"

    stats: dict[str, object] = {
        "svg": str(svg_path),
        "text_groups": 0,
        "single_line_text_objects": 0,
        "multiline_text_objects_collapsed": 0,
        "multiline_groups_not_collapsed": 0,
    }

    for group in root.iter(group_tag):
        group_id = group.attrib.get("id", "")
        if not group_id.startswith("text_"):
            continue
        stats["text_groups"] = int(stats["text_groups"]) + 1
        direct_text = [child for child in list(group) if child.tag == text_tag]
        if not direct_text:
            continue

        label_parts = [
            _svg_text_content(child).strip()
            for child in direct_text
            if _svg_text_content(child).strip()
        ]
        if label_parts:
            group.set(f"{{{INKSCAPE_NS}}}label", " / ".join(label_parts)[:240])

        if len(direct_text) == 1:
            direct_text[0].set("id", f"{group_id}_object")
            direct_text[0].set(f"{{{XML_NS}}}space", "preserve")
            stats["single_line_text_objects"] = int(stats["single_line_text_objects"]) + 1
            continue

        positions = [_svg_text_xy(child) for child in direct_text]
        if any(position is None for position in positions):
            # Keep the original Matplotlib group when a complex transform cannot be
            # represented safely as x/y tspans. The lines still remain editable text.
            stats["multiline_groups_not_collapsed"] = int(stats["multiline_groups_not_collapsed"]) + 1
            continue

        styles = [child.attrib.get("style") for child in direct_text]
        common_style = styles[0] if all(style == styles[0] for style in styles) else None
        parent_text = ET.Element(text_tag, {"id": f"{group_id}_object"})
        parent_text.set(f"{{{XML_NS}}}space", "preserve")
        if common_style:
            parent_text.set("style", common_style)

        for child, position in zip(direct_text, positions):
            assert position is not None
            tspan = ET.SubElement(parent_text, tspan_tag, {"x": position[0], "y": position[1]})
            if not common_style and child.attrib.get("style"):
                tspan.set("style", child.attrib["style"])
            for key, value in child.attrib.items():
                if key not in {"x", "y", "transform", "style", "id"}:
                    tspan.set(key, value)
            tspan.text = _svg_text_content(child)

        insertion_index = min(list(group).index(child) for child in direct_text)
        for child in direct_text:
            group.remove(child)
        group.insert(insertion_index, parent_text)
        stats["multiline_text_objects_collapsed"] = int(stats["multiline_text_objects_collapsed"]) + 1

    # Audit: every Matplotlib text group must contain one editable text object and
    # no glyph paths or glyph-use elements. Figure geometry may still use paths.
    problem_groups: list[str] = []
    total_text_elements = 0
    for group in root.iter(group_tag):
        group_id = group.attrib.get("id", "")
        if not group_id.startswith("text_"):
            continue
        descendants = list(group.iter())
        text_descendants = [element for element in descendants if element.tag == text_tag]
        total_text_elements += len(text_descendants)
        # Matplotlib may place a legitimate background-box path in a text group.
        # Outlined glyphs, however, are emitted as <use> elements and contain no
        # editable <text>. Require exactly one text object and no glyph uses.
        has_outlined_glyphs = any(element.tag == use_tag for element in descendants)
        if len(text_descendants) != 1 or has_outlined_glyphs:
            problem_groups.append(group_id)

    stats["editable_text_elements"] = total_text_elements
    stats["problem_text_groups"] = problem_groups
    stats["status"] = "PASS" if not problem_groups and total_text_elements > 0 else "FAIL"
    if stats["status"] != "PASS":
        raise RuntimeError(
            "Editable SVG text audit failed. Problem groups: " + ", ".join(problem_groups[:20])
        )

    tree.write(svg_path, encoding="utf-8", xml_declaration=True)
    return stats


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Create the Phase 15 spatial–temporal dispersal figure from existing outputs."
    )
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--output-dir", type=Path)
    p.add_argument("--canonical-spatial-iterations", "--step11g-iterations", dest="step11g_iterations", type=Path, help="Optional explicit path to the canonical frozen 05_boundary_iteration_values.csv. No alternative Step 11G file is accepted.")
    p.add_argument("--phase15-input", type=Path)
    p.add_argument("--posterior-draws", type=Path)
    p.add_argument("--primary-result", type=Path)
    p.add_argument("--dpi", type=int, default=600)
    p.add_argument(
        "--basename",
        default="Figure_4_spatial_temporal_turnover_bayesian_convention_editable_text",
        help="Output filename stem for the frozen-source revised figure.",
    )
    return p.parse_args()


def analysis_roots(project_root: Path) -> list[Path]:
    roots = []
    for name in ("04_analysis_USE _THIS", "04_analysis"):
        candidate = project_root / name
        if candidate.exists():
            roots.append(candidate)
    if not roots:
        roots = [project_root / "04_analysis_USE _THIS", project_root / "04_analysis"]
    return roots


def first_existing(paths: Iterable[Path], label: str) -> Path:
    tried = []
    for path in paths:
        path = path.expanduser().resolve()
        tried.append(path)
        if path.is_file():
            return path
    raise FileNotFoundError(f"Could not find {label}. Tried:\n" + "\n".join(map(str, tried)))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def locate_canonical_spatial_input(project_root: Path, override: Path | None) -> Path:
    if override is not None:
        path = override.expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Canonical spatial input not found: {path}")
        if path.name != "05_boundary_iteration_values.csv":
            raise RuntimeError(
                "Panel A must use the frozen canonical file named "
                "05_boundary_iteration_values.csv; refusing an alternative file."
            )
        return path

    candidates = [project_root / rel for rel in CANONICAL_SPATIAL_RELATIVE_PATHS]
    for path in candidates:
        if path.is_file():
            return path.resolve()
    raise FileNotFoundError(
        "Could not locate the frozen original 5,000-resample Panel A input. "
        "No fallback is allowed. Searched:\n" + "\n".join(str(p) for p in candidates)
    )


def discover(project_root: Path, override: Path | None, relative_candidates: list[str], filename: str, label: str) -> Path:
    if override is not None:
        return first_existing([override], label)
    roots = analysis_roots(project_root)
    candidates = [root / rel for root in roots for rel in relative_candidates]
    for path in candidates:
        if path.is_file():
            return path.resolve()
    # Final conservative fallback: exact filename only, under analysis roots.
    matches = []
    for root in roots:
        if root.exists():
            matches.extend(root.rglob(filename))
    matches = sorted({m.resolve() for m in matches if m.is_file()}, key=lambda p: (len(p.parts), str(p)))
    if matches:
        return matches[0]
    return first_existing(candidates, label)


def truthy(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})


def qsummary(values: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        raise ValueError("Cannot summarize an empty numeric vector.")
    return {
        "median": float(np.median(values)),
        "mean": float(np.mean(values)),
        "q025": float(np.quantile(values, 0.025)),
        "q055": float(np.quantile(values, 0.055)),
        "q945": float(np.quantile(values, 0.945)),
        "q975": float(np.quantile(values, 0.975)),
        "n": int(values.size),
    }


def build_spatial_contrasts(canonical_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Reconstruct Panel A only from the frozen original C3 Monte Carlo table."""
    df = pd.read_csv(canonical_path)
    required = {
        "iteration", "analysis", "scenario", "band_1", "band_2", "metric",
        "positive_value", "N0_value", "positive_minus_N0",
    }
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"Canonical turnover table is missing columns: {sorted(missing)}")

    df = df[
        (df["analysis"].astype(str) == "primary")
        & (df["scenario"].astype(str) == "C3")
        & (df["metric"].isin(["simpson_turnover", "jaccard_dissimilarity"]))
    ].copy()
    if df.empty:
        raise RuntimeError("No primary C3 rows were found in the canonical turnover table.")

    actual_boundaries = set(zip(df["band_1"].astype(str), df["band_2"].astype(str)))
    expected_boundaries = set(ADJACENT_ORDER)
    if actual_boundaries != expected_boundaries:
        raise RuntimeError(
            "Canonical turnover table boundaries differ from the four frozen adjacent contrasts. "
            f"Expected {sorted(expected_boundaries)}, found {sorted(actual_boundaries)}"
        )

    df["iteration"] = pd.to_numeric(df["iteration"], errors="raise").astype(int)
    df["positive_minus_N0"] = pd.to_numeric(df["positive_minus_N0"], errors="raise")
    key = ["iteration", "metric", "band_1", "band_2"]
    if df.duplicated(key, keep=False).any():
        raise RuntimeError("Duplicate iteration/metric/boundary rows in canonical turnover table.")

    counts = df.groupby(["metric", "band_1", "band_2"])["iteration"].nunique()
    if not (counts == 5000).all():
        raise RuntimeError(
            "Frozen Panel A requires exactly 5,000 iterations per metric and boundary. "
            f"Observed counts:\n{counts.to_string()}"
        )

    rows: list[dict[str, object]] = []
    sim_iteration = pd.DataFrame()
    for metric in ("simpson_turnover", "jaccard_dissimilarity"):
        metric_df = df[df["metric"] == metric].copy()
        for band_1, band_2 in ADJACENT_ORDER:
            vals = metric_df.loc[
                (metric_df["band_1"] == band_1) & (metric_df["band_2"] == band_2),
                "positive_minus_N0",
            ].to_numpy(float)
            s = qsummary(vals)
            rows.append({
                "section": "adjacent_boundary",
                "metric": metric,
                "band_1": band_1,
                "band_2": band_2,
                "label": f"{BAND_LABELS[band_1]} vs. {BAND_LABELS[band_2]}",
                **s,
            })

        global_by_iteration = (
            metric_df.groupby("iteration", as_index=False)["positive_minus_N0"]
            .mean()
            .rename(columns={"positive_minus_N0": "peninsula_mean_contrast"})
        )
        s = qsummary(global_by_iteration["peninsula_mean_contrast"].to_numpy(float))
        rows.append({
            "section": "peninsula_mean",
            "metric": metric,
            "band_1": "",
            "band_2": "",
            "label": "Peninsula-wide mean",
            **s,
        })
        if metric == "simpson_turnover":
            sim_iteration = global_by_iteration.copy()

        expected = FROZEN_SPATIAL_EXPECTED[metric]
        for field in ("median", "q025", "q975"):
            if abs(float(s[field]) - float(expected[field])) > FROZEN_CHECK_TOLERANCE:
                raise RuntimeError(
                    "The supplied Panel A file is not the frozen original accepted run. "
                    f"For {metric} {field}, expected {expected[field]:.6f} but reconstructed "
                    f"{s[field]:.6f}. Aborting rather than silently changing manuscript results."
                )

    return pd.DataFrame(rows), sim_iteration


def load_phase15(input_path: Path, draws_path: Path, result_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series, dict[str, float]]:
    dat = pd.read_csv(input_path)
    draws = pd.read_csv(draws_path)
    result = pd.read_csv(result_path)
    if len(result) != 1:
        raise RuntimeError(f"Expected one row in {result_path}; found {len(result)}")
    req_dat = {"pair_id", "grid_cell_id", "y_observed", "resampling_q025", "resampling_q975", "stress_change_raw"}
    req_draws = {"intercept", "stress_beta"}
    if missing := req_dat - set(dat.columns):
        raise RuntimeError(f"Phase 15A input missing columns: {sorted(missing)}")
    if missing := req_draws - set(draws.columns):
        raise RuntimeError(f"Posterior draws missing columns: {sorted(missing)}")
    stress_mean = float(dat["stress_change_raw"].mean())
    stress_sd = float(dat["stress_change_raw"].std(ddof=1))
    if not math.isfinite(stress_sd) or stress_sd <= 0:
        raise RuntimeError("Stress predictor has zero or invalid standard deviation.")
    dat["stress_change_z"] = (dat["stress_change_raw"] - stress_mean) / stress_sd
    scaling = {"stress_mean": stress_mean, "stress_sd": stress_sd}
    return dat, draws, result.iloc[0], scaling


def fitted_relationship(draws: pd.DataFrame, z_min: float, z_max: float, points: int = 200) -> pd.DataFrame:
    grid = np.linspace(z_min, z_max, points)
    intercept = draws["intercept"].to_numpy(float)[:, None]
    beta = draws["stress_beta"].to_numpy(float)[:, None]
    pred = intercept + beta * grid[None, :]
    return pd.DataFrame({
        "stress_change_z": grid,
        "median": np.median(pred, axis=0),
        "q025": np.quantile(pred, 0.025, axis=0),
        "q975": np.quantile(pred, 0.975, axis=0),
    })


def configure_matplotlib() -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10.5,
        "axes.titlesize": 12.5,
        "axes.labelsize": 11,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.titlesize": 18,
        "savefig.bbox": "tight",
        # Keep SVG lettering as editable text rather than outlining each glyph.
        "svg.fonttype": "none",
        "svg.hashsalt": "baja-phase15g-v0.5-editable-text",
        "text.usetex": False,
        # Preserve searchable/editable TrueType text in PDF as well.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def make_figure(
    spatial: pd.DataFrame,
    temporal: pd.DataFrame,
    draws: pd.DataFrame,
    result: pd.Series,
    fitted: pd.DataFrame,
    outdir: Path,
    dpi: int,
    basename: str,
) -> list[Path]:
    """Create the manuscript synthesis figure using conventional Bayesian reporting.

    Panel A preserves the frozen latitude-boundary profile and peninsula-wide
    contrasts. Panel B shows the observed temporal comparisons and posterior mean
    relationship. Panel C reports the primary standardized stress coefficient as a
    posterior median with nested 50% and 95% credible intervals. No model is refit.
    """
    configure_matplotlib()
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", ["C0", "C1", "C2"])
    primary = cycle[0]
    accent = cycle[1] if len(cycle) > 1 else "C1"
    secondary = cycle[2] if len(cycle) > 2 else "C2"
    muted = "0.42"

    fig = plt.figure(figsize=(16.6, 7.7), constrained_layout=False)
    outer = GridSpec(
        1, 3, figure=fig, width_ratios=[1.30, 1.22, 0.86],
        wspace=0.40, left=0.082, right=0.988, top=0.865, bottom=0.14,
    )
    ax_a = fig.add_subplot(outer[0, 0])
    ax_b = fig.add_subplot(outer[0, 1])
    ax_c = fig.add_subplot(outer[0, 2])

    # ------------------------------------------------------------------
    # PANEL A: frozen spatial boundary profile and peninsula summaries.
    # ------------------------------------------------------------------
    sim = spatial[(spatial.metric == "simpson_turnover") & (spatial.section == "adjacent_boundary")].copy()
    sim["order"] = sim.apply(lambda r: ADJACENT_ORDER.index((r.band_1, r.band_2)), axis=1)
    sim = sim.sort_values("order", ascending=False).reset_index(drop=True)
    pen_sim = spatial[(spatial.metric == "simpson_turnover") & (spatial.section == "peninsula_mean")].iloc[0]
    pen_jac = spatial[(spatial.metric == "jaccard_dissimilarity") & (spatial.section == "peninsula_mean")].iloc[0]

    boundary_names = {
        ("28-30N", "30-32N"): "30°N boundary\n(28–30 vs. 30–32°N)",
        ("26-28N", "28-30N"): "28°N boundary\n(26–28 vs. 28–30°N)",
        ("24-26N", "26-28N"): "26°N boundary\n(24–26 vs. 26–28°N)",
        ("23-24N", "24-26N"): "24°N boundary\n(23–24 vs. 24–26°N)",
    }
    y_boundary = np.array([5.2, 4.15, 3.10, 2.05])
    medians = sim["median"].to_numpy(float)

    ax_a.plot(medians, y_boundary, color=muted, linewidth=1.05, alpha=0.55, zorder=1)
    for idx, row in sim.iterrows():
        pair = (row["band_1"], row["band_2"])
        is_26 = pair == ("24-26N", "26-28N")
        color = accent if is_26 else primary
        ax_a.errorbar(
            row["median"], y_boundary[idx],
            xerr=[[row["median"] - row["q025"]], [row["q975"] - row["median"]]],
            fmt="o", color=color, ecolor=color, capsize=3.0,
            linewidth=1.65 if is_26 else 1.40,
            markersize=8.8 if is_26 else 6.8, zorder=3,
            markeredgecolor="white" if is_26 else color,
            markeredgewidth=0.9 if is_26 else 0.0,
        )

    y_sim, y_jac = 0.75, -0.15
    ax_a.errorbar(
        pen_sim["median"], y_sim,
        xerr=[[pen_sim["median"] - pen_sim["q025"]], [pen_sim["q975"] - pen_sim["median"]]],
        fmt="D", color=secondary, ecolor=secondary, capsize=3.5,
        linewidth=2.15, markersize=8.1, zorder=4,
    )
    ax_a.errorbar(
        pen_jac["median"], y_jac,
        xerr=[[pen_jac["median"] - pen_jac["q025"]], [pen_jac["q975"] - pen_jac["median"]]],
        fmt="s", color=muted, ecolor=muted, capsize=3.0,
        linewidth=1.35, markersize=6.4, zorder=3,
    )
    ax_a.axvline(0, linestyle="--", linewidth=1.1, color=muted, alpha=0.85)
    ax_a.axhline(1.40, linewidth=0.8, color="0.72")

    ylabels = [boundary_names[(r.band_1, r.band_2)] for _, r in sim.iterrows()]
    ax_a.set_yticks(
        list(y_boundary) + [y_sim, y_jac],
        ylabels + ["Peninsula-wide\nSimpson replacement", "Peninsula-wide\nJaccard dissimilarity"],
    )
    ax_a.set_ylim(-0.55, 5.65)
    ax_a.set_xlabel("Spatial turnover contrast (C3 − N0)")
    ax_a.set_title("A  Spatial turnover by latitude boundary", loc="left", fontweight="bold", pad=8, fontsize=11.5)
    ax_a.text(
        0.0, 0.985,
        "Negative values indicate lower C3 turnover than N0",
        transform=ax_a.transAxes, ha="left", va="top", fontsize=8.9, color=muted,
    )
    ax_a.grid(axis="x", linewidth=0.45, alpha=0.22)

    all_lo = np.r_[sim.q025.to_numpy(float), pen_sim.q025, pen_jac.q025]
    all_hi = np.r_[sim.q975.to_numpy(float), pen_sim.q975, pen_jac.q975]
    xmin, xmax = float(np.nanmin(all_lo)), float(np.nanmax(all_hi))
    span = max(xmax - xmin, 0.2)
    ax_a.set_xlim(xmin - 0.10 * span, xmax + 0.24 * span)

    row26_idx = int(sim.index[(sim["band_1"] == "24-26N") & (sim["band_2"] == "26-28N")][0])
    row26 = sim.loc[row26_idx]
    xlo, xhi = ax_a.get_xlim()
    xspan = xhi - xlo
    ax_a.annotate(
        "Local attenuation of the\nC3–N0 contrast near 26°N\n(descriptive)",
        xy=(float(row26["median"]), float(y_boundary[row26_idx])),
        xytext=(min(xhi - 0.31 * xspan, float(row26["median"]) + 0.08 * xspan), 3.60),
        ha="left", va="center", fontsize=8.4, color=accent,
        arrowprops=dict(arrowstyle="->", color=accent, linewidth=0.95),
    )
    ax_a.text(
        0.985, 0.145,
        f"Peninsula-wide Simpson contrast\n"
        f"{pen_sim['median']:+.3f} [95% interval {pen_sim['q025']:+.3f}, {pen_sim['q975']:+.3f}]",
        transform=ax_a.transAxes, ha="right", va="center", fontsize=8.7,
        bbox=dict(boxstyle="round,pad=0.28", facecolor="white", edgecolor="0.82", alpha=0.92),
    )
    ax_a.text(
        0.0, 0.945, "Profile line joins point estimates; no continuous trend was fitted.",
        transform=ax_a.transAxes, ha="left", va="top", fontsize=7.8, color=muted,
    )

    # ------------------------------------------------------------------
    # PANEL B: observed temporal contrasts and posterior mean relationship.
    # ------------------------------------------------------------------
    x = temporal["stress_change_z"].to_numpy(float)
    yy = temporal["y_observed"].to_numpy(float)
    yerr = np.vstack([
        yy - temporal["resampling_q025"].to_numpy(float),
        temporal["resampling_q975"].to_numpy(float) - yy,
    ])
    ax_b.fill_between(
        fitted["stress_change_z"], fitted["q025"], fitted["q975"],
        color=primary, alpha=0.16, linewidth=0,
        label="95% credible interval for mean response",
    )
    ax_b.plot(
        fitted["stress_change_z"], fitted["median"], color=primary, linewidth=2.1,
        label="Posterior median response",
    )
    ax_b.errorbar(
        x, yy, yerr=yerr, fmt="o", color=accent, ecolor=accent,
        capsize=2.5, linewidth=0.85, markersize=5.9, alpha=0.92,
        label="Temporal comparison (95% resampling interval)",
    )
    ax_b.axhline(0, linewidth=1.0, linestyle="--", color=muted)
    ax_b.axvline(0, linewidth=0.7, color="0.68")
    ax_b.set_xlabel("Standardized change in recent environmental stress (SD)")
    ax_b.set_ylabel("Temporal Simpson replacement contrast (C3 − N0)")
    ax_b.set_title("B  Temporal response to recent environmental stress", loc="left", fontweight="bold", pad=8, fontsize=11.5)
    ax_b.text(
        0.0, 0.985,
        "Positive values indicate greater temporal replacement in C3",
        transform=ax_b.transAxes, ha="left", va="top", fontsize=8.9, color=muted,
    )
    ax_b.grid(linewidth=0.45, alpha=0.20)
    n_pairs = len(temporal)
    n_cells = temporal["grid_cell_id"].nunique()
    ax_b.text(
        0.03, 0.90,
        f"{n_pairs} temporal comparisons across {n_cells} 25-km cells",
        transform=ax_b.transAxes, va="top", ha="left", fontsize=8.9,
        bbox=dict(boxstyle="round,pad=0.30", facecolor="white", edgecolor="0.84", alpha=0.90),
    )
    ax_b.legend(loc="lower right", frameon=False, fontsize=7.6, handlelength=2.0)

    # ------------------------------------------------------------------
    # PANEL C: conventional Bayesian coefficient display.
    # ------------------------------------------------------------------
    med = float(result["stress_beta_median"])
    lo95 = float(result["stress_beta_q025"])
    hi95 = float(result["stress_beta_q975"])
    beta_draws = pd.to_numeric(draws["stress_beta"], errors="coerce").dropna().to_numpy(float)
    lo50, hi50 = np.quantile(beta_draws, [0.25, 0.75])
    ppos = float(result["posterior_probability_beta_gt_0"])

    cmin = min(-0.08, lo95 - 0.06)
    cmax = max(0.42, hi95 + 0.06)
    y0 = 0.54
    ax_c.axvline(0, linewidth=1.15, color=muted)
    ax_c.hlines(y0, lo95, hi95, color=primary, linewidth=2.4, zorder=2)
    ax_c.hlines(y0, lo50, hi50, color=primary, linewidth=8.0, alpha=0.78, zorder=3)
    ax_c.plot(med, y0, marker="o", markersize=8.1, color=accent, zorder=4)
    ax_c.set_xlim(cmin, cmax)
    ax_c.set_ylim(0, 1)
    ax_c.set_yticks([])
    ax_c.set_xlabel("Standardized stress coefficient, β")
    ax_c.set_title("C  Bayesian stress-effect estimate", loc="left", fontweight="bold", pad=8, fontsize=11.5)
    ax_c.text(
        0.0, 0.985, "Posterior median with 50% and 95% credible intervals",
        transform=ax_c.transAxes, ha="left", va="top", fontsize=8.9, color=muted,
    )
    ax_c.grid(axis="x", linewidth=0.45, alpha=0.20)
    ax_c.text(
        0.04, 0.84,
        f"β = {med:+.3f}\n"
        f"95% CrI [{lo95:+.3f}, {hi95:+.3f}]\n"
        f"Pr(β > 0 | data) = {ppos:.3f}",
        transform=ax_c.transAxes, ha="left", va="top", fontsize=10.0,
    )
    ax_c.text(
        0.50, 0.36, "thick line: 50% CrI   |   thin line: 95% CrI",
        transform=ax_c.transAxes, ha="center", va="top", fontsize=8.1, color=muted,
    )
    ax_c.text(
        0.50, 0.11, "β > 0 indicates a greater temporal C3 response",
        transform=ax_c.transAxes, ha="center", va="bottom", fontsize=8.2, color=primary,
    )

    fig.suptitle(
        "Spatial and temporal turnover contrasts between ballooning-capable and non-ballooning arachnid assemblages",
        y=0.976, fontweight="bold", fontsize=16.6,
    )
    fig.text(
        0.5, 0.025,
        "Panel A uses the frozen 5,000-resample C3–N0 spatial analysis. Panels B–C summarize an exploratory Bayesian temporal analysis "
        "based on 15 comparisons across 12 cells; the shaded band is a credible interval for the mean response, not a prediction interval.",
        ha="center", va="bottom", fontsize=8.8,
    )

    outputs = []
    svg_audit: dict[str, object] | None = None
    metadata = {
        "Title": "Spatial and temporal turnover contrasts between ballooning-capable and non-ballooning arachnid assemblages",
        "Creator": SCRIPT_VERSION,
    }
    for ext in ("png", "svg", "pdf"):
        path = outdir / f"{basename}.{ext}"
        fig.savefig(path, dpi=dpi if ext == "png" else None, metadata=metadata)
        if ext == "svg":
            svg_audit = make_svg_text_editable(path)
            (outdir / "15G_svg_text_audit.json").write_text(
                json.dumps(svg_audit, indent=2) + "\n", encoding="utf-8"
            )
        outputs.append(path)
    plt.close(fig)
    return outputs

def main() -> int:
    args = parse_args()
    project_root = args.project_root.expanduser().resolve()
    if not project_root.is_dir():
        raise FileNotFoundError(f"Project root not found: {project_root}")
    roots = analysis_roots(project_root)
    default_out = roots[0] / "15_bayesian_h3_evidence" / "15G_spatial_temporal_figure"
    outdir = args.output_dir.expanduser().resolve() if args.output_dir else default_out
    outdir.mkdir(parents=True, exist_ok=True)

    canonical_spatial = locate_canonical_spatial_input(project_root, args.step11g_iterations)
    phase15_input = discover(
        project_root, args.phase15_input,
        ["15_bayesian_h3_evidence/15A_input_audit/15A_bayesian_model_input.csv"],
        "15A_bayesian_model_input.csv", "Phase 15A model input",
    )
    posterior_draws = discover(
        project_root, args.posterior_draws,
        ["15_bayesian_h3_evidence/15B_primary_bayesian_model/15B_posterior_draws.csv"],
        "15B_posterior_draws.csv", "Phase 15B posterior draws",
    )
    primary_result = discover(
        project_root, args.primary_result,
        ["15_bayesian_h3_evidence/15B_primary_bayesian_model/15B_primary_result.csv"],
        "15B_primary_result.csv", "Phase 15B primary result",
    )
    spatial, spatial_iteration = build_spatial_contrasts(canonical_spatial)
    temporal, draws, result, scaling = load_phase15(phase15_input, posterior_draws, primary_result)
    zmin = float(temporal["stress_change_z"].min()) - 0.15
    zmax = float(temporal["stress_change_z"].max()) + 0.15
    fitted = fitted_relationship(draws, zmin, zmax)
    outputs = make_figure(spatial, temporal, draws, result, fitted, outdir, args.dpi, args.basename)

    spatial.to_csv(outdir / "15G_spatial_contrast_values.csv", index=False)
    temporal.to_csv(outdir / "15G_temporal_plot_values.csv", index=False)
    fitted.to_csv(outdir / "15G_posterior_fitted_relationship.csv", index=False)
    provenance = {
        "phase": "15G",
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "project_root": str(project_root),
        "inputs": {
            "canonical_spatial_iterations": str(canonical_spatial),
            "canonical_spatial_sha256": sha256(canonical_spatial),
            "phase15_input": str(phase15_input),
            "posterior_draws": str(posterior_draws),
            "primary_result": str(primary_result),
        },
        "stress_scaling": scaling,
        "outputs": [str(p) for p in outputs],
        "notes": [
            "No model is refitted in Phase 15G.",
            "Panel A is locked to the original accepted 5,000-resample C3 file: 05_boundary_iteration_values.csv.",
            "The script aborts if the reconstructed frozen Simpson or Jaccard summaries differ beyond 1e-5.",
            "Peninsula-wide spatial contrasts are the mean of four adjacent contrasts within each canonical Monte Carlo iteration.",
            "The temporal fitted relationship holds log-event and transition-midpoint covariates at their means and cell effect at zero.",
            "The posterior ribbon is uncertainty in the mean relationship, not a predictive interval for new observations.",
            "The 26°N annotation describes the boundary point estimate; it is not a separately confirmed breakpoint.",
            "Panel C follows conventional Bayesian coefficient reporting: posterior median with nested 50% and 95% credible intervals.",
            "The SVG preserves fonts as editable text; each Matplotlib text artist is one SVG text object, with multiline strings represented by tspans.",
            "The SVG text audit must pass with no outlined glyph paths inside Matplotlib text groups.",
            "Predictive validation, prior sensitivity, and posterior-density diagnostics remain in the supplementary Phase 15 outputs.",
        ],
    }
    (outdir / "15G_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    readme = f"""PHASE 15G v0.5 — CONVENTIONAL BAYESIAN SYNTHESIS FIGURE WITH EDITABLE SVG TEXT
================================================
Status: COMPLETED

This figure combines two completed analyses without refitting them:
1. The original accepted 5,000-resample C3/N0 spatial turnover run;
2. Phase 15A/15B temporal observations and primary Bayesian stress coefficient.

The main figure reports the Bayesian result using a posterior point estimate with
nested 50% and 95% credible intervals. The SVG keeps each label as editable text;
multiline labels created from one string are collapsed into one text object with tspans. Prediction, prior sensitivity, posterior-density,
and power diagnostics remain in the supplementary Phase 15 outputs.

Primary output:
{outputs[0]}

Editable SVG:
{next(path for path in outputs if path.suffix == ".svg")}

SVG text audit:
{outdir / "15G_svg_text_audit.json"}

Frozen Panel A source:
{canonical_spatial}
SHA-256: {sha256(canonical_spatial)}

Frozen peninsula-wide targets:
Simpson = -0.141745 [-0.275743, -0.013406]
Jaccard = -0.007552 [-0.084820, 0.072382]

Interpretive guardrail:
The boundary profile is descriptive and the peninsula-wide Simpson contrast is the primary spatial result. The temporal result is
an exploratory association based on {len(temporal)} comparisons across
{temporal['grid_cell_id'].nunique()} cells. They answer complementary but distinct
questions and must not be described as a single causal process.
"""
    (outdir / "15G_README.txt").write_text(readme, encoding="utf-8")
    print(readme)
    print("OUTPUTS:")
    for path in outputs:
        print(path)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"PHASE 15G FAILED: {exc}", file=sys.stderr)
        raise
