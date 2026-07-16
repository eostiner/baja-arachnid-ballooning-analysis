#!/usr/bin/env python3
"""
Step 10G — final biogeographic figure refinement and equal-cell
ballooning:non-ballooning ratio by latitude band.

Primary comparison:
  Ballooning = C3 = D1 + D2 + D3
  Non-ballooning = fixed N0 reference
  D4 is excluded.

Outputs:
  * map-only ecoregion synthesis figure (no representative-taxon table),
  * equal-cell ballooning:non-ballooning richness ratio by latitude band,
  * combined map + ratio figure,
  * 5,000-iteration draw table and adjacent-band random-split tests.

The ratio is computed from pooled UNIQUE genera after sampling the same number
of occupied 25-km cells (22) from each latitude band. The plotted statistic is
log2(ballooning richness / non-ballooning richness), which is symmetric around
zero: 0 = equal richness, +1 = 2:1, -1 = 1:2. Direct B:N ratios and ballooning
shares are also written to CSV.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import traceback
from pathlib import Path
from typing import Any

try:
    import numpy as np
    import pandas as pd
    import geopandas as gpd
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, Wedge, Circle
    from matplotlib import patheffects
except Exception as exc:
    raise SystemExit(
        "Missing packages. Activate the Baja environment and install: "
        "pandas numpy matplotlib geopandas shapely pyproj\n"
        f"Original error: {exc}"
    )

EXPECTED_CELLS = 205
EXPECTED_GENERA = 267
EQUAL_CELLS = 22
FORMAL_MIN_CELLS = 5
DEFAULT_ITERATIONS = 5000

BANDS = [
    ("23–24°N", 23.0, 24.0),
    ("24–26°N", 24.0, 26.0),
    ("26–28°N", 26.0, 28.0),
    ("28–30°N", 28.0, 30.0),
    ("30–32°N", 30.0, 32.0000001),
]
BOUNDARIES = [24.0, 26.0, 28.0, 30.0]

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

DONUT_OFFSETS = {
    "Matorral Costero": (-0.20, 0.22),
    "Chaparral": (0.18, 0.08),
    "Matorral Costero Rosetófilo": (-0.18, -0.05),
    "Desierto de San Felipe": (0.18, 0.08),
    "Matorrales Tropicales": (0.13, -0.05),
    "Planicies de Magdalena": (-0.12, -0.04),
    "Sierra de la Giganta": (0.16, 0.02),
    "Costa Central del Golfo": (0.24, 0.05),
    "Desierto de Vizcaíno": (-0.15, 0.00),
    "Desierto Central": (0.00, 0.08),
}


def norm(x: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(x).lower())


def clean_str(x: Any) -> str:
    return "" if pd.isna(x) else str(x).strip()


def bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    return s.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes", "y"})


def md5_file(path: Path, block: int = 1024 * 1024) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        while True:
            b = fh.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


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
    for p in preferred:
        if norm(p) in lookup:
            return lookup[norm(p)]
    if contains:
        for c in df.columns:
            if contains in norm(c):
                return str(c)
    raise KeyError(f"Could not identify required column from {preferred}; columns={list(df.columns)}")


def find_existing(candidates: list[Path], role: str) -> Path:
    for p in candidates:
        if p.exists() and p.is_file():
            return p
    raise FileNotFoundError(f"Could not locate {role}. Checked:\n" + "\n".join(str(p) for p in candidates))


def input_paths(project_root: Path) -> dict[str, Path]:
    bio = project_root / "04_analysis" / "C3_pipeline_rebuild" / "09_C3_biogeographic_concordance"
    incidence_candidates = [
        project_root / "02_data_clean" / "08_grid25km_incidence" / "10_ballooning_final_genus_grid25km_incidence_long.csv",
        project_root / "02_data_clean" / "08_grid25km_incidence" / "10_genus_cell_incidence_long.csv",
    ]
    trait_candidates = [
        bio / "10C_equal_cell_ecoregion_richness" / "10C_genus_trait_lookup_used.csv",
        project_root / "04_analysis" / "C3_pipeline_rebuild" / "01_trait_merge" / "C3_authoritative_trait_master.csv",
    ]
    for manifest in [
        bio / "10E_published_break_concordance" / "10E_input_manifest.csv",
        bio / "10D_ecoregion_boundary_turnover" / "10D_input_manifest.csv",
        bio / "10C_equal_cell_ecoregion_richness" / "10C_input_manifest.csv",
    ]:
        if manifest.exists():
            m = read_csv(manifest)
            if {"role", "file"}.issubset(m.columns):
                for role_name in ["genus_cell_incidence", "incidence"]:
                    rows = m.loc[m["role"].astype(str) == role_name, "file"]
                    if len(rows):
                        incidence_candidates.insert(0, Path(str(rows.iloc[0])))
                for role_name in ["genus_trait_lookup", "trait"]:
                    rows = m.loc[m["role"].astype(str) == role_name, "file"]
                    if len(rows):
                        trait_candidates.insert(0, Path(str(rows.iloc[0])))
    return {
        "ecoregions": find_existing([
            bio / "10A_ecoregion_gis_audit" / "10A_ecoregions_validated_mainland_only.gpkg",
        ], "Step 10A mainland ecoregion GIS"),
        "crosswalk": find_existing([
            bio / "10B_cell_ecoregion_crosswalk" / "10B_cell_ecoregion_crosswalk.csv",
        ], "Step 10B cell crosswalk"),
        "ecoregion_summary": find_existing([
            bio / "10C_equal_cell_ecoregion_richness" / "10C_equal_cell_summary.csv",
        ], "Step 10C equal-cell ecoregion summary"),
        "breaks": find_existing([
            bio / "10E_published_break_concordance" / "10E_break_registry_used.csv",
            bio / "10E_published_break_concordance" / "10E_PUBLISHED_BREAK_REGISTRY.csv",
        ], "Step 10E break registry"),
        "incidence": find_existing(incidence_candidates, "final genus-by-cell incidence"),
        "traits": find_existing(trait_candidates, "C3/N0 trait lookup"),
    }


def standardize_traits(path: Path) -> pd.DataFrame:
    raw = read_csv(path)
    genus_col = detect_col(raw, ["genus", "accepted_genus", "final_genus"], "genus")
    evidence_col = None
    for preferred, contains in [
        (["evidence_class", "exclusive_tier", "ballooning_evidence_class", "evidence_tier"], "evidence"),
        (["exclusive_tier"], "tier"),
    ]:
        try:
            evidence_col = detect_col(raw, preferred, contains)
            break
        except KeyError:
            pass
    class_col = None
    for preferred, contains in [
        (["analysis_class", "c3_analysis_class", "primary_class", "primary_C3_group"], "analysisclass"),
        (["primary_C3_group"], "primaryc3group"),
    ]:
        try:
            class_col = detect_col(raw, preferred, contains)
            break
        except KeyError:
            pass
    if evidence_col is None and class_col is None:
        raise KeyError(f"Trait table lacks a recognized evidence or primary-class column: {list(raw.columns)}")

    out = pd.DataFrame({"genus": raw[genus_col].map(clean_str)})
    if evidence_col is not None:
        ev = raw[evidence_col].map(clean_str).str.upper()
        ev = ev.str.extract(r"(?<![A-Z0-9])(D1|D2|D3|D4|N0|C3)(?![A-Z0-9])", expand=False).fillna(ev)
        out["evidence_class"] = ev
        out["analysis_class"] = np.where(ev.isin(["D1", "D2", "D3", "C3"]), "C3",
                                  np.where(ev.eq("N0"), "N0",
                                  np.where(ev.eq("D4"), "D4_excluded", "")))
    else:
        out["evidence_class"] = ""
        cls = raw[class_col].map(clean_str).str.upper()
        out["analysis_class"] = np.where(cls.str.contains("N0|NON"), "N0",
                                  np.where(cls.str.contains("C3|BALLOON"), "C3",
                                  np.where(cls.str.contains("D4"), "D4_excluded", "")))
    bad = out[~out["analysis_class"].isin(["C3", "N0", "D4_excluded"])]
    if len(bad):
        raise RuntimeError(f"Unresolved trait assignments for {len(bad)} rows; examples={bad.head(10).to_dict('records')}")
    conflicts = out.groupby("genus")["analysis_class"].nunique()
    if (conflicts > 1).any():
        raise RuntimeError(f"Conflicting classes for genera: {conflicts[conflicts > 1].index.tolist()[:20]}")
    return out[["genus", "analysis_class"]].drop_duplicates("genus")


def standardize_incidence(path: Path, valid_cells: set[str]) -> pd.DataFrame:
    d = read_csv(path)
    cell_col = detect_col(d, ["grid25km_id", "step10b_cell_id", "cell_id", "grid_id"], "cellid")
    genus_col = detect_col(d, ["genus", "accepted_genus", "final_genus"], "genus")
    presence_col = None
    for c in d.columns:
        if norm(c) in {"presence", "present", "incidence", "occupied", "pa"}:
            presence_col = str(c)
            break
    keep = [cell_col, genus_col] + ([presence_col] if presence_col else [])
    x = d[keep].copy()
    x.columns = ["cell_id", "genus"] + (["presence"] if presence_col else [])
    x["cell_id"] = x["cell_id"].map(clean_str)
    x["genus"] = x["genus"].map(clean_str)
    if "presence" in x:
        p = pd.to_numeric(x["presence"], errors="coerce")
        x = x[p.fillna(0) > 0] if p.notna().mean() >= 0.8 else x[bool_series(x["presence"])]
    x = x[x["cell_id"].isin(valid_cells)]
    x = x[(x["genus"] != "") & (~x["genus"].str.lower().isin({"na", "nan", "unknown", "unidentified"}))]
    return x[["cell_id", "genus"]].drop_duplicates()


def pooled(cell_ids: np.ndarray | list[str], sets_by_cell: dict[str, set[str]]) -> set[str]:
    ans: set[str] = set()
    for cid in cell_ids:
        ans.update(sets_by_cell.get(str(cid), set()))
    return ans


def q(x: pd.Series | np.ndarray, prob: float) -> float:
    z = np.asarray(x, dtype=float)
    z = z[np.isfinite(z)]
    return float(np.quantile(z, prob)) if len(z) else np.nan


def mean_finite(x: pd.Series | np.ndarray) -> float:
    z = np.asarray(x, dtype=float)
    z = z[np.isfinite(z)]
    return float(np.mean(z)) if len(z) else np.nan


def bh_adjust(pvals: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvals, errors="coerce").to_numpy(float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    vals = p[ok]
    if not len(vals):
        return pd.Series(out, index=pvals.index)
    order = np.argsort(vals)
    ranked = vals[order]
    adj = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    rev = np.empty_like(adj)
    rev[order] = adj
    out[np.where(ok)[0]] = rev
    return pd.Series(out, index=pvals.index)


def prepare_data(paths: dict[str, Path]) -> dict[str, Any]:
    cross0 = read_csv(paths["crosswalk"])
    id_col = detect_col(cross0, ["step10b_cell_id", "grid25km_id", "cell_id"], "cellid")
    lat_col = detect_col(cross0, ["centroid_latitude", "latitude", "lat"], "latitude")
    lon_col = detect_col(cross0, ["centroid_longitude", "longitude", "lon"], "longitude")
    c3_col = detect_col(cross0, ["C3_positive_genera", "C3"], "c3")
    n0_col = detect_col(cross0, ["N0_reference_genera", "N0"], "n0")
    cross = pd.DataFrame({
        "cell_id": cross0[id_col].map(clean_str),
        "centroid_latitude": pd.to_numeric(cross0[lat_col], errors="coerce"),
        "centroid_longitude": pd.to_numeric(cross0[lon_col], errors="coerce"),
        "C3_frozen": pd.to_numeric(cross0[c3_col], errors="coerce").fillna(0).astype(int),
        "N0_frozen": pd.to_numeric(cross0[n0_col], errors="coerce").fillna(0).astype(int),
    }).drop_duplicates("cell_id")
    if len(cross) != EXPECTED_CELLS:
        raise RuntimeError(f"Frozen cell count failed: expected {EXPECTED_CELLS}, found {len(cross)}")
    if cross[["centroid_latitude", "centroid_longitude"]].isna().any().any():
        raise RuntimeError("Missing cell centroid coordinates.")

    traits = standardize_traits(paths["traits"])
    incidence = standardize_incidence(paths["incidence"], set(cross.cell_id))
    inc = incidence.merge(traits, on="genus", how="left")
    if inc["analysis_class"].isna().any():
        missing = inc.loc[inc["analysis_class"].isna(), "genus"].drop_duplicates().tolist()
        raise RuntimeError(f"Incidence genera missing from trait table: {missing[:20]}")
    if inc["genus"].nunique() != EXPECTED_GENERA:
        raise RuntimeError(f"Frozen genus count failed: expected {EXPECTED_GENERA}, found {inc['genus'].nunique()}")

    set_maps: dict[str, dict[str, set[str]]] = {}
    for cls in ["C3", "N0"]:
        z = inc[inc.analysis_class == cls]
        groups = z.groupby("cell_id")["genus"].apply(lambda s: set(s.astype(str))).to_dict()
        set_maps[cls] = {cid: groups.get(cid, set()) for cid in cross.cell_id}
    validation = cross.copy()
    validation["C3_reconstructed"] = validation.cell_id.map(lambda c: len(set_maps["C3"][c]))
    validation["N0_reconstructed"] = validation.cell_id.map(lambda c: len(set_maps["N0"][c]))
    validation["C3_exact"] = validation.C3_frozen == validation.C3_reconstructed
    validation["N0_exact"] = validation.N0_frozen == validation.N0_reconstructed
    if not validation.C3_exact.all() or not validation.N0_exact.all():
        raise RuntimeError("Frozen per-cell C3/N0 count validation failed.")

    return {
        "cross": cross,
        "incidence": inc,
        "set_maps": set_maps,
        "validation": validation,
        "ecoregions": gpd.read_file(paths["ecoregions"]),
        "ecoregion_summary": read_csv(paths["ecoregion_summary"]),
        "breaks": read_csv(paths["breaks"]),
    }


def build_ecoregion_values(summary: pd.DataFrame, crosswalk_path: Path) -> pd.DataFrame:
    cross = read_csv(crosswalk_path)
    id_col = detect_col(cross, ["step10b_cell_id", "grid25km_id", "cell_id"], "cellid")
    cross = cross.rename(columns={id_col: "cell_id"})
    eligible_col = detect_col(cross, ["primary_assignment_eligible"], "primaryassignmenteligible")
    eco_col = detect_col(cross, ["dominant_ecoregion", "ecoregion"], "ecoregion")
    cross[eligible_col] = bool_series(cross[eligible_col])
    counts = (cross[cross[eligible_col]].groupby(eco_col)["cell_id"].nunique()
              .rename("n_primary_cells").reset_index().rename(columns={eco_col: "ecoregion"}))
    s = summary[summary["analysis_set"].astype(str).str.lower().eq("primary")].copy()
    piv = s.pivot_table(index="ecoregion", columns="metric", values="mean", aggfunc="first").reset_index()
    v = counts.merge(piv, on="ecoregion", how="left")
    v["formal_comparison"] = v.n_primary_cells >= FORMAL_MIN_CELLS
    for c in ["C3_richness", "N0_richness"]:
        if c not in v:
            v[c] = np.nan
    v["primary_contrast_richness"] = v.C3_richness + v.N0_richness
    v["ballooning_share"] = np.where(v.primary_contrast_richness > 0,
                                      v.C3_richness / v.primary_contrast_richness, np.nan)
    return v


def band_ratio_analysis(cross: pd.DataFrame, set_maps: dict[str, dict[str, set[str]]],
                        seed: int, iterations: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ids_by_band: dict[str, np.ndarray] = {}
    count_rows = []
    for name, lo, hi in BANDS:
        ids = cross.loc[(cross.centroid_latitude >= lo) & (cross.centroid_latitude < hi), "cell_id"].to_numpy(str)
        if len(ids) < EQUAL_CELLS:
            raise RuntimeError(f"Band {name} has {len(ids)} cells; frozen equal n is {EQUAL_CELLS}")
        ids_by_band[name] = ids
        count_rows.append({"latitude_band": name, "available_cells": len(ids), "equal_cells": EQUAL_CELLS})

    rng = np.random.default_rng(seed)
    rows: list[dict[str, Any]] = []
    for it in range(iterations):
        for name, _, _ in BANDS:
            ids = rng.choice(ids_by_band[name], size=EQUAL_CELLS, replace=False)
            c3 = len(pooled(ids, set_maps["C3"]))
            n0 = len(pooled(ids, set_maps["N0"]))
            if c3 == 0 or n0 == 0:
                ratio = (c3 + 0.5) / (n0 + 0.5)
                correction = True
            else:
                ratio = c3 / n0
                correction = False
            rows.append({
                "iteration": it + 1,
                "latitude_band": name,
                "equal_cells": EQUAL_CELLS,
                "ballooning_richness": c3,
                "non_ballooning_richness": n0,
                "ballooning_to_non_ballooning_ratio": ratio,
                "log2_ratio": math.log2(ratio),
                "ballooning_share": c3 / (c3 + n0) if (c3 + n0) > 0 else np.nan,
                "zero_correction_used": correction,
            })
    draws = pd.DataFrame(rows)
    sums = []
    for name, _, _ in BANDS:
        g = draws[draws.latitude_band == name]
        sums.append({
            "latitude_band": name,
            "available_cells": len(ids_by_band[name]),
            "equal_cells": EQUAL_CELLS,
            "iterations": iterations,
            "mean_ballooning_richness": mean_finite(g.ballooning_richness),
            "mean_non_ballooning_richness": mean_finite(g.non_ballooning_richness),
            "mean_log2_ratio": mean_finite(g.log2_ratio),
            "median_log2_ratio": float(np.median(g.log2_ratio.to_numpy(float))),
            "log2_ratio_q025": q(g.log2_ratio, 0.025),
            "log2_ratio_q975": q(g.log2_ratio, 0.975),
            "geometric_mean_B_to_N_ratio": 2 ** mean_finite(g.log2_ratio),
            "median_B_to_N_ratio": 2 ** float(np.median(g.log2_ratio.to_numpy(float))),
            "ratio_q025": 2 ** q(g.log2_ratio, 0.025),
            "ratio_q975": 2 ** q(g.log2_ratio, 0.975),
            "mean_ballooning_share": mean_finite(g.ballooning_share),
            "ballooning_share_q025": q(g.ballooning_share, 0.025),
            "ballooning_share_q975": q(g.ballooning_share, 0.975),
            "zero_correction_fraction": float(g.zero_correction_used.mean()),
        })
    summaries = pd.DataFrame(sums)

    # Adjacent-band tests: each iteration samples 22 actual cells from each band,
    # then compares the observed north-south log-ratio difference with a random
    # 22/22 split of the same local 44-cell pool. This mirrors Step 10E's local
    # random-split logic and avoids comparing unequal sampling effort.
    rng2 = np.random.default_rng(seed + 1000003)
    test_draws: list[dict[str, Any]] = []
    for idx, boundary in enumerate(BOUNDARIES):
        south_name = BANDS[idx][0]
        north_name = BANDS[idx + 1][0]
        south_ids_all = ids_by_band[south_name]
        north_ids_all = ids_by_band[north_name]
        for it in range(iterations):
            south_ids = rng2.choice(south_ids_all, size=EQUAL_CELLS, replace=False)
            north_ids = rng2.choice(north_ids_all, size=EQUAL_CELLS, replace=False)
            selected = np.concatenate([south_ids, north_ids])
            rng2.shuffle(selected)
            null_s = selected[:EQUAL_CELLS]
            null_n = selected[EQUAL_CELLS:]

            def lr(ids: np.ndarray) -> float:
                c = len(pooled(ids, set_maps["C3"]))
                n = len(pooled(ids, set_maps["N0"]))
                rr = (c + 0.5) / (n + 0.5) if (c == 0 or n == 0) else c / n
                return math.log2(rr)

            obs_diff = lr(north_ids) - lr(south_ids)
            null_diff = lr(null_n) - lr(null_s)
            test_draws.append({
                "boundary_latitude": boundary,
                "south_band": south_name,
                "north_band": north_name,
                "iteration": it + 1,
                "observed_north_minus_south_log2_ratio": obs_diff,
                "null_random_split_log2_ratio_difference": null_diff,
            })
    contrast_draws = pd.DataFrame(test_draws)
    tests = []
    for (boundary, south_name, north_name), g in contrast_draws.groupby(
        ["boundary_latitude", "south_band", "north_band"], sort=False
    ):
        obs = g.observed_north_minus_south_log2_ratio.to_numpy(float)
        null = g.null_random_split_log2_ratio_difference.to_numpy(float)
        obs_mean = mean_finite(obs)
        p = float((1 + np.sum(np.abs(null) >= abs(obs_mean))) / (len(null) + 1))
        tests.append({
            "boundary_latitude": boundary,
            "south_band": south_name,
            "north_band": north_name,
            "observed_mean_north_minus_south_log2_ratio": obs_mean,
            "observed_q025": q(obs, 0.025),
            "observed_q975": q(obs, 0.975),
            "fold_change_in_B_to_N_ratio_north_vs_south": 2 ** obs_mean,
            "null_mean": mean_finite(null),
            "random_split_two_sided_p": p,
        })
    tests = pd.DataFrame(tests)
    tests["BH_q_four_adjacent_boundaries"] = bh_adjust(tests.random_split_two_sided_p)
    return draws, summaries, contrast_draws, tests, pd.DataFrame(count_rows)


def draw_donut(ax, x: float, y: float, total: float, ballooning: float, non_ballooning: float,
               max_total: float, colors: dict[str, Any]) -> None:
    if not np.isfinite(total) or total <= 0:
        return
    radius = 0.16 + 0.20 * math.sqrt(total / max_total)
    vals = [max(ballooning, 0), max(non_ballooning, 0)]
    s = sum(vals)
    start = 90.0
    for val, cls in zip(vals, ["C3", "N0"]):
        theta = 360.0 * val / s if s > 0 else 0
        ax.add_patch(Wedge((x, y), radius, start, start + theta,
                           facecolor=colors[cls], edgecolor="white", linewidth=0.6, zorder=8))
        start += theta
    ax.add_patch(Circle((x, y), radius, facecolor="none", edgecolor="0.15", linewidth=0.75, zorder=9))
    ax.add_patch(Circle((x, y), radius * 0.44, facecolor="white", edgecolor="0.25", linewidth=0.45, zorder=10))
    ax.text(x, y, f"{total:.0f}", ha="center", va="center", fontsize=7.2, fontweight="bold", zorder=11)


def plot_map(ax, data: dict[str, Any], values: pd.DataFrame, show_title: bool = True) -> None:
    eco = data["ecoregions"].copy()
    label_col = "ecoregion_label" if "ecoregion_label" in eco else detect_col(eco, ["ecoregion_label", "Nombre", "name"], "ecoregion")
    eco = eco.rename(columns={label_col: "ecoregion"})
    labels = sorted(eco.ecoregion.astype(str).unique())
    cmap = plt.get_cmap("tab20")
    eco_colors = {label: cmap(i % 20) for i, label in enumerate(labels)}
    class_colors = {"C3": plt.get_cmap("tab10")(1), "N0": plt.get_cmap("tab10")(2)}

    for label in labels:
        g = eco[eco.ecoregion.astype(str) == label]
        vv = values[values.ecoregion == label]
        formal = bool(vv.formal_comparison.iloc[0]) if len(vv) else False
        face = eco_colors[label] if formal else (0.90, 0.90, 0.90, 1)
        g.plot(ax=ax, facecolor=face, edgecolor="0.25", linewidth=0.55,
               hatch=None if formal else "///", zorder=1)

    breaks = data["breaks"].sort_values("anchor_latitude")
    for (_, r), letter in zip(breaks.iterrows(), ["D", "C", "B", "A"]):
        lat = float(r.anchor_latitude)
        ax.axhline(lat, color="0.12", linestyle=(0, (5, 3)), linewidth=1.15, zorder=4)
        ax.text(-117.23, lat + 0.03, letter, fontsize=8.6, fontweight="bold", ha="center", va="center",
                bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="0.2", lw=0.8), zorder=10)

    max_total = float(values.primary_contrast_richness.max(skipna=True))
    for _, er in eco.iterrows():
        label = str(er.ecoregion)
        p = er.geometry.representative_point()
        dx, dy = DONUT_OFFSETS.get(label, (0.0, 0.0))
        x, y = p.x + dx, p.y + dy
        vv = values[values.ecoregion == label]
        if len(vv) and bool(vv.formal_comparison.iloc[0]) and np.isfinite(vv.primary_contrast_richness.iloc[0]):
            row = vv.iloc[0]
            draw_donut(ax, x, y, float(row.primary_contrast_richness),
                       float(row.C3_richness), float(row.N0_richness), max_total, class_colors)
            ax.text(x, y - 0.31, ABBREVIATIONS.get(label, label[:4]), fontsize=7.0, fontweight="bold",
                    ha="center", va="top",
                    path_effects=[patheffects.withStroke(linewidth=2.3, foreground="white")], zorder=12)
        else:
            ax.text(x, y, ABBREVIATIONS.get(label, label[:4]), fontsize=6.9, fontweight="bold",
                    ha="center", va="center", color="0.35",
                    path_effects=[patheffects.withStroke(linewidth=2.3, foreground="white")], zorder=12)

    ax.set_xlim(-117.6, -109.0)
    ax.set_ylim(22.75, 32.9)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(color="0.90", linewidth=0.45, zorder=0)
    if show_title:
        ax.set_title("Baja ecoregions, published breaks, and standardized dispersal-group richness",
                     fontsize=12.5, pad=10)

    # Compact legend inside blank northeast margin.
    x0, y0 = -110.1, 32.55
    ax.add_patch(Rectangle((x0 - 0.15, y0 - 1.45), 1.15, 1.35, facecolor="white", edgecolor="0.45",
                           linewidth=0.7, alpha=0.94, zorder=20))
    ax.add_patch(Rectangle((x0, y0 - 0.25), 0.18, 0.16, facecolor=class_colors["C3"], edgecolor="0.2", zorder=21))
    ax.text(x0 + 0.24, y0 - 0.17, "Ballooning (C3: D1–D3)", fontsize=7.4, va="center", zorder=21)
    ax.add_patch(Rectangle((x0, y0 - 0.52), 0.18, 0.16, facecolor=class_colors["N0"], edgecolor="0.2", zorder=21))
    ax.text(x0 + 0.24, y0 - 0.44, "Non-ballooning (fixed N0)", fontsize=7.4, va="center", zorder=21)
    ax.text(x0, y0 - 0.73, "Donut area = expected C3 + N0\nrichness from 8 occupied cells", fontsize=7.0, va="top", zorder=21)
    ax.text(x0, y0 - 1.12, "Dashed lines = a priori breaks\n24°, 26°, 28°, and 30°N", fontsize=7.0, va="top", zorder=21)


def plot_ratio(ax, summaries: pd.DataFrame, tests: pd.DataFrame, annotate_tests: bool = True) -> None:
    order = [b[0] for b in BANDS]
    s = summaries.set_index("latitude_band").loc[order].reset_index()
    x = np.arange(len(s))
    y = s.median_log2_ratio.to_numpy(float)
    lo = s.log2_ratio_q025.to_numpy(float)
    hi = s.log2_ratio_q975.to_numpy(float)
    ax.errorbar(x, y, yerr=[y - lo, hi - y], fmt="o", capsize=4, linewidth=1.4, markersize=6)
    ax.axhline(0, color="0.25", linewidth=1.0, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels(order)
    ax.set_ylabel("Ballooning : non-ballooning richness ratio (log₂ scale)")
    ax.set_xlabel("Latitude band")
    ax.set_title("Equal-cell dispersal-group balance by latitude band", fontsize=12.5, pad=10)
    ax.grid(axis="y", color="0.90", linewidth=0.6)

    # Ratio-equivalent y ticks.
    ticks = np.array([-2, -1, 0, 1, 2], dtype=float)
    ymin = min(float(lo.min()) - 0.25, -1.25)
    ymax = max(float(hi.max()) + 0.35, 0.75)
    ax.set_ylim(ymin, ymax)
    ticks = ticks[(ticks >= ymin) & (ticks <= ymax)]
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{2**t:.2g}:1" if t >= 0 else f"1:{2**(-t):.2g}" for t in ticks])

    for i, row in s.iterrows():
        ratio = float(row.median_B_to_N_ratio)
        share = 100 * float(row.mean_ballooning_share)
        ax.text(i, float(row.log2_ratio_q975) + 0.08, f"{ratio:.2f}:1\n({share:.0f}% B)",
                ha="center", va="bottom", fontsize=8.0)

    if annotate_tests and len(tests):
        # Show only corrected significance; full raw and adjusted results are in CSV.
        sig = tests[pd.to_numeric(tests.BH_q_four_adjacent_boundaries, errors="coerce") < 0.05]
        base = ymax - 0.10
        for j, r in enumerate(sig.itertuples(index=False)):
            idx = BOUNDARIES.index(float(r.boundary_latitude))
            yy = base - 0.16 * j
            ax.plot([idx, idx, idx + 1, idx + 1], [yy - 0.03, yy, yy, yy - 0.03], color="0.2", lw=0.8)
            ax.text(idx + 0.5, yy + 0.015, f"q={r.BH_q_four_adjacent_boundaries:.3f}", ha="center", va="bottom", fontsize=7.5)

    ax.text(0.01, 0.02,
            "Each point: 5,000 draws of 22 occupied cells per band; unique genera pooled within each draw.\n"
            "0 (1:1) = equal richness; positive values favor ballooning; negative values favor non-ballooning.",
            transform=ax.transAxes, fontsize=7.5, va="bottom")


def make_figures(out_dir: Path, data: dict[str, Any], values: pd.DataFrame,
                 summaries: pd.DataFrame, tests: pd.DataFrame, dpi: int) -> None:
    pub = out_dir / "publication_outputs"
    pub.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.6, 11.0))
    plot_map(ax, data, values, show_title=False)
    fig.suptitle("Baja California arachnid dispersal biogeography", fontsize=17, fontweight="bold", y=0.985)
    fig.text(0.5, 0.958, "Ballooning (C3: D1–D3) versus non-ballooning (fixed N0) across independently mapped ecoregions",
             ha="center", fontsize=9.6)
    fig.tight_layout(rect=[0.02, 0.02, 0.98, 0.925])
    fig.savefig(pub / "10G_biogeographic_map_only.png", dpi=dpi, bbox_inches="tight")
    fig.savefig(pub / "10G_biogeographic_map_only.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.5, 6.2))
    plot_ratio(ax, summaries, tests, annotate_tests=True)
    fig.tight_layout()
    fig.savefig(pub / "10G_equal_cell_ballooning_nonballooning_ratio_by_band.png", dpi=dpi, bbox_inches="tight")
    fig.savefig(pub / "10G_equal_cell_ballooning_nonballooning_ratio_by_band.pdf", bbox_inches="tight")
    plt.close(fig)

    fig = plt.figure(figsize=(14.0, 9.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], left=0.055, right=0.98,
                          top=0.91, bottom=0.08, wspace=0.18)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    plot_map(ax1, data, values, show_title=False)
    plot_ratio(ax2, summaries, tests, annotate_tests=True)
    ax1.set_title("A. Ecoregion synthesis map", fontsize=12.5, loc="left", pad=10)
    ax2.set_title("B. Equal-cell ballooning : non-ballooning ratio", fontsize=12.5, loc="left", pad=10)
    fig.suptitle("Baja arachnid biogeography and dispersal-group balance", fontsize=17, fontweight="bold")
    fig.savefig(pub / "10G_map_and_band_ratio_combined.png", dpi=dpi, bbox_inches="tight")
    fig.savefig(pub / "10G_map_and_band_ratio_combined.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root")
    parser.add_argument("seed", nargs="?", type=int, default=20260715)
    parser.add_argument("iterations", nargs="?", type=int, default=DEFAULT_ITERATIONS)
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    root = Path(args.project_root).expanduser().resolve()
    out_dir = root / "04_analysis" / "C3_pipeline_rebuild" / "09_C3_biogeographic_concordance" / "10G_band_ratio_synthesis"
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        paths = input_paths(root)
        manifest = pd.DataFrame([
            {"role": role, "file": str(path), "md5": md5_file(path), "size_bytes": path.stat().st_size}
            for role, path in paths.items()
        ])
        manifest.to_csv(out_dir / "10G_input_manifest.csv", index=False)

        data = prepare_data(paths)
        data["validation"].to_csv(out_dir / "10G_frozen_count_validation.csv", index=False)
        values = build_ecoregion_values(data["ecoregion_summary"], paths["crosswalk"])
        values.to_csv(out_dir / "10G_ecoregion_map_values.csv", index=False)

        draws, summaries, contrast_draws, tests, counts = band_ratio_analysis(
            data["cross"], data["set_maps"], args.seed, args.iterations
        )
        counts.to_csv(out_dir / "10G_band_cell_counts.csv", index=False)
        draws.to_csv(out_dir / "10G_band_ratio_draws.csv", index=False)
        summaries.to_csv(out_dir / "10G_band_ratio_summaries.csv", index=False)
        contrast_draws.to_csv(out_dir / "10G_adjacent_band_ratio_null_draws.csv", index=False)
        tests.to_csv(out_dir / "10G_adjacent_band_ratio_tests.csv", index=False)
        make_figures(out_dir, data, values, summaries, tests, args.dpi)

        summary = {
            "step": "10G",
            "audit_status": "PASS_FROZEN_COUNTS_REPRODUCED_AND_BAND_RATIO_COMPLETE",
            "seed": args.seed,
            "iterations": args.iterations,
            "occupied_cells": int(len(data["cross"])),
            "total_genera": int(data["incidence"].genus.nunique()),
            "equal_cells_per_band": EQUAL_CELLS,
            "bands": [b[0] for b in BANDS],
            "primary_ratio_statistic": "log2(unique pooled C3 richness / unique pooled N0 richness)",
            "adjacent_test": "two-sided local random 22/22 split of the same sampled 44-cell two-band pool",
            "multiple_test_correction": "Benjamini-Hochberg across four adjacent boundaries",
            "D4": "excluded",
            "representative_taxa_table": "removed",
        }
        (out_dir / "10G_analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        (out_dir / "README_RESULTS_FIRST.txt").write_text(
            "STEP 10G COMPLETE\n\n"
            "Start with:\n"
            "  publication_outputs/10G_biogeographic_map_only.png\n"
            "  publication_outputs/10G_equal_cell_ballooning_nonballooning_ratio_by_band.png\n"
            "  publication_outputs/10G_map_and_band_ratio_combined.png\n"
            "  10G_band_ratio_summaries.csv\n"
            "  10G_adjacent_band_ratio_tests.csv\n\n"
            "Interpretation: the ratio uses only ballooning C3 and fixed non-ballooning N0 genera. "
            "D4 and the representative-taxon table are excluded. The log2 scale is symmetric around a 1:1 ratio.\n",
            encoding="utf-8",
        )
        print("STEP 10G COMPLETE")
        print("AUDIT_STATUS=PASS_FROZEN_COUNTS_REPRODUCED_AND_BAND_RATIO_COMPLETE")
        print(f"OUTPUT_DIR={out_dir}")
        print(f"OCCUPIED_CELLS={len(data['cross'])}")
        print(f"TOTAL_GENERA={data['incidence'].genus.nunique()}")
        print(f"EQUAL_CELLS_PER_BAND={EQUAL_CELLS}")
        print(f"ITERATIONS={args.iterations}")
        return 0
    except Exception as exc:
        (out_dir / "10G_FAILURE.txt").write_text(
            f"STEP 10G FAILED\n{type(exc).__name__}: {exc}\n\n{traceback.format_exc()}", encoding="utf-8"
        )
        print(f"STEP 10G FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    import sys
    raise SystemExit(main())
