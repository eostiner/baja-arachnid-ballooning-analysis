#!/usr/bin/env python3
"""
Step 10D — neighboring-cell turnover across independently mapped Baja ecoregion boundaries.

Primary inferential set:
  - Step 10B primary-assignment-eligible occupied cells
  - both cells belong to the 10 Step 10C formal ecoregions (>=5 cells)
  - rook-adjacent 25-km grid cells (shared grid edge; no diagonals)

Sensitivity set:
  - Step 10B unambiguous cells
  - same formal 10 ecoregions and rook adjacency

Assemblages:
  - all genera
  - C3 = D1 + D2 + D3
  - N0 = fixed non-ballooning reference

Turnover metrics:
  - Jaccard dissimilarity = (b + c) / (a + b + c)
  - Simpson replacement = min(b, c) / (a + min(b, c))

Inference:
  - observed mean turnover difference: cross-boundary minus within-ecoregion
  - ecoregion labels are shuffled among cells within the five 2-degree latitude bands
  - cell assemblages, adjacency network, edge dependencies, and within-band label counts remain fixed
  - one-sided and two-sided permutation P values, plus BH-FDR within each analysis set/filter

The script refuses to proceed unless it reproduces the frozen 205-cell, 267-genus,
C3, and N0 counts used in Steps 10B–10C.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
import traceback
from collections import defaultdict
from pathlib import Path
from typing import Any

try:
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
except Exception as exc:
    raise SystemExit(
        "Missing Python packages. Activate the Baja environment and install: "
        "pandas numpy matplotlib\nOriginal error: %s" % exc
    )

STEP = "10D"
EXPECTED_CELLS = 205
EXPECTED_GENERA = 267
DEFAULT_PERMUTATIONS = 5000
FORMAL_MIN_CELLS = 5


def clean_str(x: Any) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip()


def bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    return s.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes", "y"})


def md5_file(path: Path, block: int = 1024 * 1024) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        while True:
            b = f.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def read_csv_flexible(path: Path) -> pd.DataFrame:
    last = None
    for enc in ("utf-8-sig", "utf-8", "latin1"):
        try:
            return pd.read_csv(path, encoding=enc, low_memory=False)
        except Exception as exc:
            last = exc
    raise RuntimeError(f"Could not read {path}: {last}")


def find_existing(paths: list[Path], role: str) -> Path:
    for p in paths:
        if p.exists() and p.is_file():
            return p
    raise FileNotFoundError(f"Could not locate required {role}. Tried:\n" + "\n".join(str(p) for p in paths))


def source_paths(project_root: Path) -> dict[str, Path]:
    bio_root = project_root / "04_analysis/C3_pipeline_rebuild/09_C3_biogeographic_concordance"
    crosswalk = find_existing([
        bio_root / "10B_cell_ecoregion_crosswalk/10B_cell_ecoregion_crosswalk.csv",
    ], "Step 10B crosswalk")

    incidence_candidates = [
        project_root / "02_data_clean/08_grid25km_incidence/10_ballooning_final_genus_grid25km_incidence_long.csv",
    ]
    trait_candidates = [
        project_root / "04_analysis/C3_pipeline_rebuild/01_trait_merge/C3_authoritative_trait_master.csv",
    ]

    # Reuse the source paths already validated by Step 10C when available.
    manifest = bio_root / "10C_equal_cell_ecoregion_richness/10C_input_manifest.csv"
    if manifest.exists():
        m = read_csv_flexible(manifest)
        if {"role", "file"}.issubset(m.columns):
            for role, target in [("genus_cell_incidence", incidence_candidates),
                                 ("genus_trait_lookup", trait_candidates)]:
                rows = m.loc[m["role"].astype(str) == role, "file"]
                if len(rows):
                    target.insert(0, Path(str(rows.iloc[0])))

    incidence = find_existing(incidence_candidates, "genus-by-cell incidence table")
    trait = find_existing(trait_candidates, "C3/N0 genus trait lookup")

    sample_sizes = find_existing([
        bio_root / "10C_equal_cell_ecoregion_richness/10C_ecoregion_sample_sizes.csv",
    ], "Step 10C ecoregion sample-size table")

    return {
        "crosswalk": crosswalk,
        "incidence": incidence,
        "trait": trait,
        "sample_sizes": sample_sizes,
    }


def identify_column(df: pd.DataFrame, candidates: list[str], contains: list[str] | None = None) -> str:
    norm = {re.sub(r"[^a-z0-9]+", "", str(c).lower()): str(c) for c in df.columns}
    for c in candidates:
        k = re.sub(r"[^a-z0-9]+", "", c.lower())
        if k in norm:
            return norm[k]
    if contains:
        for col in df.columns:
            n = re.sub(r"[^a-z0-9]+", "", str(col).lower())
            if all(x in n for x in contains):
                return str(col)
    raise KeyError(f"Could not identify required column from {candidates}. Available: {list(df.columns)}")


def load_incidence(path: Path, valid_ids: set[str]) -> pd.DataFrame:
    d = read_csv_flexible(path)
    cell_col = identify_column(d, ["grid25km_id", "cell_id", "step10b_cell_id", "grid_id"], ["cell", "id"])
    genus_col = identify_column(d, ["genus", "final_genus", "accepted_genus"], ["genus"])
    keep = [cell_col, genus_col]
    presence_col = None
    for c in d.columns:
        n = re.sub(r"[^a-z0-9]+", "", str(c).lower())
        if n in {"presence", "present", "incidence", "occupied", "pa"}:
            presence_col = str(c)
            keep.append(presence_col)
            break
    x = d[keep].copy()
    x.columns = ["cell_id", "genus"] + (["presence"] if presence_col else [])
    x["cell_id"] = x["cell_id"].map(clean_str)
    x["genus"] = x["genus"].map(clean_str)
    if "presence" in x.columns:
        p = pd.to_numeric(x["presence"], errors="coerce")
        if p.notna().mean() >= 0.8:
            x = x[p.fillna(0) > 0]
        else:
            x = x[bool_series(x["presence"])]
    x = x[x["cell_id"].isin(valid_ids)]
    x = x[(x["genus"] != "") & (~x["genus"].str.lower().isin({"na", "nan", "unknown", "unidentified"}))]
    return x[["cell_id", "genus"]].drop_duplicates().reset_index(drop=True)


TRAIT_TOKEN_RE = re.compile(r"(?<![A-Z0-9])(D1|D2|D3|D4|N0|C3)(?![A-Z0-9])", re.I)


def normalize_trait_class(value: Any) -> str | None:
    text = clean_str(value).upper()
    if not text:
        return None
    tokens = {token.upper() for token in TRAIT_TOKEN_RE.findall(text)}
    if tokens and tokens.issubset({"D1", "D2", "D3", "C3"}):
        return "C3"
    if tokens == {"N0"}:
        return "N0"
    if tokens == {"D4"}:
        return "D4_excluded"
    normalized = re.sub(r"[^a-z0-9]+", "", text.casefold())
    if normalized in {
        "c3", "primaryc3", "d1d2d3", "d1tod3", "ballooning",
        "ballooningc3", "c3ballooning",
    }:
        return "C3"
    if normalized in {
        "n0", "nonballooning", "fixednonballooning",
        "referencenonballooning", "nonballooningreference", "noballooning",
    }:
        return "N0"
    if normalized in {"d4", "d4excluded", "excludedd4"}:
        return "D4_excluded"
    return None


def load_traits(path: Path) -> pd.DataFrame:
    d = read_csv_flexible(path)
    genus_col = identify_column(d, ["genus", "final_genus", "accepted_genus"], ["genus"])
    analysis_col = None
    evidence_col = None
    for c in d.columns:
        n = re.sub(r"[^a-z0-9]+", "", str(c).lower())
        if n in {"analysisclass", "c3analysisclass", "traitanalysisclass", "primaryclass"}:
            analysis_col = str(c)
        if n in {
            "exclusivetier", "evidenceclass", "evidencetier",
            "ballooningevidencetier", "ballooningevidenceclass",
        }:
            evidence_col = str(c)
    source_col = evidence_col or analysis_col
    if source_col is None:
        raise KeyError(
            "Trait table has neither an explicit analysis-class nor a D1/D2/D3/D4/N0 evidence-tier column."
        )
    out = pd.DataFrame({
        "genus": d[genus_col].map(clean_str),
        "analysis_class": d[source_col].map(normalize_trait_class),
    })
    unresolved = out[(out["genus"] != "") & out["analysis_class"].isna()]
    if len(unresolved):
        raise RuntimeError(
            f"Unresolved explicit trait classes for {len(unresolved)} rows from {source_col}; "
            f"examples={unresolved.head(10).to_dict('records')}"
        )
    out = out[out["genus"] != ""].copy()
    conflicts = out.groupby("genus")["analysis_class"].nunique()
    if (conflicts > 1).any():
        raise RuntimeError(
            "Conflicting trait classes for genera: "
            + ", ".join(conflicts[conflicts > 1].index.astype(str).tolist()[:20])
        )
    return out.drop_duplicates("genus").reset_index(drop=True)


GRID_RE = re.compile(r"(?:^|_)C(?P<c>[+-]?\d+)_R(?P<r>[+-]?\d+)(?:$|_)", re.IGNORECASE)


def parse_grid_id(cell_id: str) -> tuple[int, int] | None:
    m = GRID_RE.search(cell_id)
    if not m:
        return None
    return int(m.group("c")), int(m.group("r"))


def haversine_km(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    r = 6371.0088
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


def build_rook_edges(cw: pd.DataFrame) -> pd.DataFrame:
    coords: dict[tuple[int, int], str] = {}
    parsed = {}
    for cid in cw["step10b_cell_id"].astype(str):
        cr = parse_grid_id(cid)
        parsed[cid] = cr
        if cr is not None:
            if cr in coords:
                raise RuntimeError(f"Duplicate parsed grid coordinate {cr}: {coords[cr]} and {cid}")
            coords[cr] = cid
    failed = [k for k, v in parsed.items() if v is None]
    if failed:
        raise RuntimeError(f"Could not parse C/R grid coordinates from {len(failed)} cell IDs. Examples: {failed[:5]}")

    meta = cw.set_index("step10b_cell_id", drop=False)
    rows = []
    for (c, r), cid1 in sorted(coords.items()):
        for dc, dr in [(1, 0), (0, 1)]:
            cid2 = coords.get((c + dc, r + dr))
            if cid2 is None:
                continue
            a, b = sorted([cid1, cid2])
            m1, m2 = meta.loc[a], meta.loc[b]
            dist = haversine_km(float(m1["centroid_longitude"]), float(m1["centroid_latitude"]),
                                float(m2["centroid_longitude"]), float(m2["centroid_latitude"]))
            rows.append({
                "cell_i": a,
                "cell_j": b,
                "grid_step": "rook",
                "centroid_distance_km": dist,
                "latitude_i": float(m1["centroid_latitude"]),
                "latitude_j": float(m2["centroid_latitude"]),
                "longitude_i": float(m1["centroid_longitude"]),
                "longitude_j": float(m2["centroid_longitude"]),
                "latitude_band_i": clean_str(m1["latitude_band"]),
                "latitude_band_j": clean_str(m2["latitude_band"]),
                "ecoregion_i": clean_str(m1["dominant_ecoregion"]),
                "ecoregion_j": clean_str(m2["dominant_ecoregion"]),
                "primary_i": bool(m1["primary_assignment_eligible"]),
                "primary_j": bool(m2["primary_assignment_eligible"]),
                "unambiguous_i": bool(m1["sensitivity_unambiguous_eligible"]),
                "unambiguous_j": bool(m2["sensitivity_unambiguous_eligible"]),
            })
    e = pd.DataFrame(rows).drop_duplicates(["cell_i", "cell_j"]).reset_index(drop=True)
    if e.empty:
        raise RuntimeError("No rook-adjacent occupied-cell pairs were reconstructed.")
    return e


def beta_components(s1: set[str], s2: set[str]) -> tuple[int, int, int, float, float, float]:
    a = len(s1 & s2)
    b = len(s1 - s2)
    c = len(s2 - s1)
    union = a + b + c
    jac = (b + c) / union if union > 0 else np.nan
    den = a + min(b, c)
    sim = min(b, c) / den if den > 0 else np.nan
    nes = jac - sim if np.isfinite(jac) and np.isfinite(sim) else np.nan
    return a, b, c, jac, sim, nes


def add_turnover(edges: pd.DataFrame, incidence: pd.DataFrame, traits: pd.DataFrame) -> pd.DataFrame:
    trait_map = traits.set_index("genus")["analysis_class"].to_dict()
    sets_all = incidence.groupby("cell_id")["genus"].apply(set).to_dict()
    sets_c3 = incidence[incidence["genus"].map(trait_map).eq("C3")].groupby("cell_id")["genus"].apply(set).to_dict()
    sets_n0 = incidence[incidence["genus"].map(trait_map).eq("N0")].groupby("cell_id")["genus"].apply(set).to_dict()
    group_sets = {"all": sets_all, "C3": sets_c3, "N0": sets_n0}

    out_rows = []
    for _, row in edges.iterrows():
        rec = row.to_dict()
        rec["crosses_ecoregion_boundary"] = rec["ecoregion_i"] != rec["ecoregion_j"]
        rec["ecoregion_pair"] = " | ".join(sorted([rec["ecoregion_i"], rec["ecoregion_j"]]))
        rec["pair_mean_latitude"] = (rec["latitude_i"] + rec["latitude_j"]) / 2
        rec["pair_crosses_latitude_band"] = rec["latitude_band_i"] != rec["latitude_band_j"]
        for g, ss in group_sets.items():
            s1 = ss.get(rec["cell_i"], set())
            s2 = ss.get(rec["cell_j"], set())
            a, b, c, jac, sim, nes = beta_components(s1, s2)
            prefix = g
            rec[f"{prefix}_richness_i"] = len(s1)
            rec[f"{prefix}_richness_j"] = len(s2)
            rec[f"{prefix}_shared_a"] = a
            rec[f"{prefix}_unique_i_b"] = b
            rec[f"{prefix}_unique_j_c"] = c
            rec[f"{prefix}_union_positive"] = (a + b + c) > 0
            rec[f"{prefix}_both_nonempty"] = len(s1) > 0 and len(s2) > 0
            rec[f"{prefix}_jaccard"] = jac
            rec[f"{prefix}_simpson"] = sim
            rec[f"{prefix}_nestedness_residual"] = nes
        out_rows.append(rec)
    return pd.DataFrame(out_rows)


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) == 0 or len(y) == 0:
        return np.nan
    # Edge counts are small; exact pairwise comparison is transparent and deterministic.
    comp = np.sign(x[:, None] - y[None, :])
    return float(comp.mean())


def bh_fdr(p: pd.Series) -> pd.Series:
    vals = pd.to_numeric(p, errors="coerce").to_numpy(float)
    out = np.full(len(vals), np.nan)
    ok = np.isfinite(vals)
    idx = np.where(ok)[0]
    if not len(idx):
        return pd.Series(out, index=p.index)
    order = idx[np.argsort(vals[idx])]
    ranked = vals[order] * len(order) / np.arange(1, len(order) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out[order] = np.minimum(ranked, 1.0)
    return pd.Series(out, index=p.index)


def run_permutations(
    pair_df: pd.DataFrame,
    cells: pd.DataFrame,
    analysis_set: str,
    eligible_col: str,
    formal_ecoregions: set[str],
    permutations: int,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, np.ndarray]]:
    eligible_cells = cells[bool_series(cells[eligible_col]) & cells["dominant_ecoregion"].isin(formal_ecoregions)].copy()
    ids = set(eligible_cells["step10b_cell_id"].astype(str))
    e = pair_df[pair_df["cell_i"].isin(ids) & pair_df["cell_j"].isin(ids)].copy().reset_index(drop=True)
    if e.empty:
        raise RuntimeError(f"No adjacent pairs in analysis set {analysis_set}.")

    node_ids = eligible_cells["step10b_cell_id"].astype(str).tolist()
    node_index = {c: i for i, c in enumerate(node_ids)}
    i_idx = e["cell_i"].map(node_index).to_numpy(int)
    j_idx = e["cell_j"].map(node_index).to_numpy(int)
    labels = eligible_cells["dominant_ecoregion"].astype(str).to_numpy(object)
    bands = eligible_cells["latitude_band"].astype(str).to_numpy(object)
    band_indices = [np.where(bands == b)[0] for b in sorted(set(bands))]
    observed_cross = labels[i_idx] != labels[j_idx]

    specs = []
    for assemblage in ["all", "C3", "N0"]:
        filters = ["union_positive"] if assemblage == "all" else ["union_positive", "both_nonempty"]
        for pair_filter in filters:
            valid_filter = e[f"{assemblage}_{pair_filter}"].astype(bool).to_numpy()
            for metric in ["jaccard", "simpson"]:
                values = pd.to_numeric(e[f"{assemblage}_{metric}"], errors="coerce").to_numpy(float)
                valid = valid_filter & np.isfinite(values)
                specs.append((assemblage, pair_filter, metric, values, valid))

    nulls: dict[str, np.ndarray] = {}
    for assemblage, pair_filter, metric, values, valid in specs:
        key = f"{analysis_set}::{assemblage}::{pair_filter}::{metric}"
        nulls[key] = np.full(permutations, np.nan)

    for pidx in range(permutations):
        perm_labels = labels.copy()
        for idx in band_indices:
            if len(idx) > 1:
                perm_labels[idx] = rng.permutation(perm_labels[idx])
        cross = perm_labels[i_idx] != perm_labels[j_idx]
        for assemblage, pair_filter, metric, values, valid in specs:
            key = f"{analysis_set}::{assemblage}::{pair_filter}::{metric}"
            cx = valid & cross
            wi = valid & (~cross)
            if cx.any() and wi.any():
                nulls[key][pidx] = float(np.nanmean(values[cx]) - np.nanmean(values[wi]))

    test_rows = []
    summary_rows = []
    for assemblage, pair_filter, metric, values, valid in specs:
        key = f"{analysis_set}::{assemblage}::{pair_filter}::{metric}"
        cx = valid & observed_cross
        wi = valid & (~observed_cross)
        x = values[cx]
        y = values[wi]
        obs = float(np.mean(x) - np.mean(y)) if len(x) and len(y) else np.nan
        null = nulls[key][np.isfinite(nulls[key])]
        p_greater = (1 + np.sum(null >= obs)) / (len(null) + 1) if np.isfinite(obs) and len(null) else np.nan
        p_two = (1 + np.sum(np.abs(null) >= abs(obs))) / (len(null) + 1) if np.isfinite(obs) and len(null) else np.nan
        test_rows.append({
            "analysis_set": analysis_set,
            "assemblage": assemblage,
            "pair_filter": pair_filter,
            "metric": metric,
            "eligible_cells": len(eligible_cells),
            "adjacent_pairs": len(e),
            "valid_pairs": int(valid.sum()),
            "cross_boundary_pairs": int(cx.sum()),
            "within_ecoregion_pairs": int(wi.sum()),
            "cross_mean": float(np.mean(x)) if len(x) else np.nan,
            "within_mean": float(np.mean(y)) if len(y) else np.nan,
            "mean_difference_cross_minus_within": obs,
            "cross_median": float(np.median(x)) if len(x) else np.nan,
            "within_median": float(np.median(y)) if len(y) else np.nan,
            "median_difference_cross_minus_within": float(np.median(x) - np.median(y)) if len(x) and len(y) else np.nan,
            "cliffs_delta_cross_vs_within": cliffs_delta(x, y),
            "null_mean": float(np.mean(null)) if len(null) else np.nan,
            "null_q025": float(np.quantile(null, 0.025)) if len(null) else np.nan,
            "null_q975": float(np.quantile(null, 0.975)) if len(null) else np.nan,
            "permutations_valid": len(null),
            "p_one_sided_cross_greater": p_greater,
            "p_two_sided": p_two,
        })
        for boundary, arr in [("cross_boundary", x), ("within_ecoregion", y)]:
            summary_rows.append({
                "analysis_set": analysis_set,
                "assemblage": assemblage,
                "pair_filter": pair_filter,
                "metric": metric,
                "pair_type": boundary,
                "n_pairs": len(arr),
                "mean": float(np.mean(arr)) if len(arr) else np.nan,
                "median": float(np.median(arr)) if len(arr) else np.nan,
                "sd": float(np.std(arr, ddof=1)) if len(arr) > 1 else np.nan,
                "q025": float(np.quantile(arr, 0.025)) if len(arr) else np.nan,
                "q975": float(np.quantile(arr, 0.975)) if len(arr) else np.nan,
            })

    tests = pd.DataFrame(test_rows)
    # Multiplicity adjustment within each analysis set and pair-filter family over 6 metric/assemblage tests.
    tests["q_one_sided_BH"] = np.nan
    tests["q_two_sided_BH"] = np.nan
    for pf, idx in tests.groupby("pair_filter").groups.items():
        tests.loc[idx, "q_one_sided_BH"] = bh_fdr(tests.loc[idx, "p_one_sided_cross_greater"]).values
        tests.loc[idx, "q_two_sided_BH"] = bh_fdr(tests.loc[idx, "p_two_sided"]).values

    e["analysis_set"] = analysis_set
    e["formal_analysis_pair"] = True
    e["crosses_ecoregion_boundary"] = observed_cross
    return tests, pd.DataFrame(summary_rows), nulls


def boundary_pair_summary(edges: pd.DataFrame, formal_ecoregions: set[str], cells: pd.DataFrame) -> pd.DataFrame:
    eligible = cells[bool_series(cells["primary_assignment_eligible"]) & cells["dominant_ecoregion"].isin(formal_ecoregions)]
    ids = set(eligible["step10b_cell_id"].astype(str))
    e = edges[edges["cell_i"].isin(ids) & edges["cell_j"].isin(ids) & edges["crosses_ecoregion_boundary"]].copy()
    rows = []
    for ep, g in e.groupby("ecoregion_pair"):
        base = {
            "ecoregion_pair": ep,
            "n_cross_boundary_edges": len(g),
            "mean_latitude": g["pair_mean_latitude"].mean(),
            "min_latitude": g["pair_mean_latitude"].min(),
            "max_latitude": g["pair_mean_latitude"].max(),
        }
        for assemblage in ["all", "C3", "N0"]:
            for metric in ["jaccard", "simpson"]:
                vals = pd.to_numeric(g[f"{assemblage}_{metric}"], errors="coerce").dropna()
                base[f"{assemblage}_{metric}_n"] = len(vals)
                base[f"{assemblage}_{metric}_mean"] = vals.mean() if len(vals) else np.nan
                base[f"{assemblage}_{metric}_median"] = vals.median() if len(vals) else np.nan
        rows.append(base)
    return pd.DataFrame(rows).sort_values(["all_simpson_mean", "n_cross_boundary_edges"], ascending=[False, False]) if rows else pd.DataFrame()


def make_effect_plot(tests: pd.DataFrame, metric: str, out_png: Path, out_pdf: Path) -> None:
    d = tests[(tests["metric"] == metric) & (tests["pair_filter"] == "union_positive")].copy()
    order = []
    labels = []
    for aset in ["primary_formal10", "unambiguous_formal10"]:
        for assemblage in ["all", "C3", "N0"]:
            order.append((aset, assemblage))
            labels.append(f"{assemblage}\n{'Primary' if aset.startswith('primary') else 'Unambiguous'}")
    d["_key"] = list(zip(d["analysis_set"], d["assemblage"]))
    d = d.set_index("_key").reindex(order).reset_index()
    x = np.arange(len(d))
    y = d["mean_difference_cross_minus_within"].to_numpy(float)
    lo = d["null_q025"].to_numpy(float)
    hi = d["null_q975"].to_numpy(float)
    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    ax.vlines(x, lo, hi, linewidth=3, alpha=0.65, label="95% latitude-stratified null interval")
    ax.scatter(x, y, s=60, zorder=3, label="Observed cross − within mean")
    ax.axhline(0, linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel(f"Across-boundary minus within-ecoregion mean {metric.title()} dissimilarity")
    ax.set_title(f"Step 10D: ecoregion-boundary concordance — {metric.title()}")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)


def make_edge_map(edges: pd.DataFrame, cells: pd.DataFrame, formal_ecoregions: set[str], out_png: Path, out_pdf: Path) -> None:
    eligible = cells[bool_series(cells["primary_assignment_eligible"]) & cells["dominant_ecoregion"].isin(formal_ecoregions)].copy()
    ids = set(eligible["step10b_cell_id"].astype(str))
    e = edges[edges["cell_i"].isin(ids) & edges["cell_j"].isin(ids) & edges["crosses_ecoregion_boundary"]].copy()
    fig, ax = plt.subplots(figsize=(6.5, 10))
    ax.scatter(eligible["centroid_longitude"], eligible["centroid_latitude"], s=10, alpha=0.5)
    vals = pd.to_numeric(e["all_simpson"], errors="coerce")
    finite = vals[np.isfinite(vals)]
    vmax = float(finite.max()) if len(finite) else 1.0
    for _, r in e.iterrows():
        v = float(r["all_simpson"]) if np.isfinite(r["all_simpson"]) else 0.0
        lw = 0.5 + 3.0 * (v / vmax if vmax > 0 else 0)
        ax.plot([r["longitude_i"], r["longitude_j"]], [r["latitude_i"], r["latitude_j"]], linewidth=lw, alpha=0.75)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("Cross-ecoregion rook edges\nLine width = total Simpson replacement")
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)


def write_results_readme(out_dir: Path, summary: dict[str, Any], tests: pd.DataFrame) -> None:
    primary = tests[(tests["analysis_set"] == "primary_formal10") & (tests["pair_filter"] == "union_positive")]
    lines = [
        "STEP 10D RESULTS — READ THIS FIRST",
        "===================================",
        "",
        f"Audit status: {summary['audit_status']}",
        f"Occupied cells validated: {summary['occupied_cells']}",
        f"Genera validated: {summary['total_genera']}",
        f"All occupied rook-neighbor pairs: {summary['all_rook_pairs']}",
        f"Primary formal-10 adjacent pairs: {summary['primary_formal10_pairs']}",
        f"Unambiguous formal-10 adjacent pairs: {summary['unambiguous_formal10_pairs']}",
        f"Permutations: {summary['permutations']}",
        "",
        "PRIMARY UNION-POSITIVE TESTS",
        "----------------------------",
    ]
    for _, r in primary.iterrows():
        lines.append(
            f"{r['assemblage']} {r['metric']}: cross mean={r['cross_mean']:.4f}; "
            f"within mean={r['within_mean']:.4f}; difference={r['mean_difference_cross_minus_within']:.4f}; "
            f"P(one-sided)={r['p_one_sided_cross_greater']:.5f}; q(BH)={r['q_one_sided_BH']:.5f}"
        )
    lines += [
        "",
        "Interpret the unambiguous set and both-nonempty C3/N0 tests as sensitivity analyses.",
        "The permutation shuffles ecoregion labels among cells only within latitude bands; assemblages and the adjacency network remain fixed.",
        "Sparse ecoregions (<5 primary cells) are retained in the complete pair table but excluded from formal inference.",
        "",
        "Key files:",
        "  10D_permutation_tests.csv",
        "  10D_turnover_group_summaries.csv",
        "  10D_neighbor_pair_metrics.csv",
        "  10D_cross_boundary_ecoregion_pair_summary.csv",
        "  publication_outputs/10D_*_boundary_effect.*",
        "  publication_outputs/10D_cross_boundary_simpson_map.*",
    ]
    (out_dir / "README_RESULTS_FIRST.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("project_root", type=Path)
    ap.add_argument("seed", type=int, nargs="?", default=20260715)
    ap.add_argument("permutations", type=int, nargs="?", default=DEFAULT_PERMUTATIONS)
    args = ap.parse_args()

    project_root = args.project_root.expanduser().resolve()
    out_dir = project_root / "04_analysis/C3_pipeline_rebuild/09_C3_biogeographic_concordance/10D_ecoregion_boundary_turnover"
    pub_dir = out_dir / "publication_outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    pub_dir.mkdir(parents=True, exist_ok=True)

    paths = source_paths(project_root)
    manifest_rows = [{"role": k, "file": str(v), "md5": md5_file(v)} for k, v in paths.items()]
    pd.DataFrame(manifest_rows).to_csv(out_dir / "10D_input_manifest.csv", index=False)

    cw = read_csv_flexible(paths["crosswalk"])
    required = {
        "step10b_cell_id", "dominant_ecoregion", "primary_assignment_eligible",
        "sensitivity_unambiguous_eligible", "latitude_band", "centroid_latitude",
        "centroid_longitude", "C3_positive_genera", "N0_reference_genera"
    }
    missing = required - set(cw.columns)
    if missing:
        raise RuntimeError(f"Step 10B crosswalk is missing required fields: {sorted(missing)}")
    cw["step10b_cell_id"] = cw["step10b_cell_id"].astype(str)
    cw["primary_assignment_eligible"] = bool_series(cw["primary_assignment_eligible"])
    cw["sensitivity_unambiguous_eligible"] = bool_series(cw["sensitivity_unambiguous_eligible"])
    if len(cw) != EXPECTED_CELLS or cw["step10b_cell_id"].nunique() != EXPECTED_CELLS:
        raise RuntimeError(f"Frozen occupied-cell validation failed: expected {EXPECTED_CELLS}, found {len(cw)} rows / {cw['step10b_cell_id'].nunique()} IDs")
    if int(cw["primary_assignment_eligible"].sum()) != 193:
        raise RuntimeError(f"Expected 193 primary-eligible cells, found {int(cw['primary_assignment_eligible'].sum())}")
    if int(cw["sensitivity_unambiguous_eligible"].sum()) != 183:
        raise RuntimeError(f"Expected 183 unambiguous cells, found {int(cw['sensitivity_unambiguous_eligible'].sum())}")

    incidence = load_incidence(paths["incidence"], set(cw["step10b_cell_id"]))
    traits = load_traits(paths["trait"])
    if incidence["cell_id"].nunique() != EXPECTED_CELLS:
        raise RuntimeError(f"Incidence table has {incidence['cell_id'].nunique()} occupied cells; expected {EXPECTED_CELLS}")
    if incidence["genus"].nunique() != EXPECTED_GENERA:
        raise RuntimeError(f"Incidence table has {incidence['genus'].nunique()} genera; expected {EXPECTED_GENERA}")
    if traits["genus"].nunique() != EXPECTED_GENERA:
        raise RuntimeError(f"Trait lookup has {traits['genus'].nunique()} genera; expected {EXPECTED_GENERA}")
    missing_traits = sorted(set(incidence["genus"]) - set(traits["genus"]))
    if missing_traits:
        raise RuntimeError(f"Incidence genera missing from trait lookup: {missing_traits[:20]}")

    trait_map = traits.set_index("genus")["analysis_class"].to_dict()
    c3_counts = incidence[incidence["genus"].map(trait_map).eq("C3")].groupby("cell_id")["genus"].nunique()
    n0_counts = incidence[incidence["genus"].map(trait_map).eq("N0")].groupby("cell_id")["genus"].nunique()
    validation = cw[["step10b_cell_id", "C3_positive_genera", "N0_reference_genera"]].copy()
    validation["reconstructed_C3"] = validation["step10b_cell_id"].map(c3_counts).fillna(0).astype(int)
    validation["reconstructed_N0"] = validation["step10b_cell_id"].map(n0_counts).fillna(0).astype(int)
    validation["C3_exact"] = validation["reconstructed_C3"] == pd.to_numeric(validation["C3_positive_genera"], errors="coerce").fillna(0).astype(int)
    validation["N0_exact"] = validation["reconstructed_N0"] == pd.to_numeric(validation["N0_reference_genera"], errors="coerce").fillna(0).astype(int)
    validation.to_csv(out_dir / "10D_frozen_count_validation.csv", index=False)
    if not validation["C3_exact"].all() or not validation["N0_exact"].all():
        raise RuntimeError(
            f"Frozen trait counts failed: C3 exact={validation['C3_exact'].mean():.4f}, "
            f"N0 exact={validation['N0_exact'].mean():.4f}"
        )

    ss = read_csv_flexible(paths["sample_sizes"])
    if not {"ecoregion", "primary_cells", "formal_comparison"}.issubset(ss.columns):
        raise RuntimeError("Step 10C sample-size table does not have expected columns.")
    formal = set(ss.loc[bool_series(ss["formal_comparison"]), "ecoregion"].astype(str))
    if len(formal) != 10:
        raise RuntimeError(f"Expected 10 formal ecoregions, found {len(formal)}: {sorted(formal)}")

    edges = build_rook_edges(cw)
    # The 25-km grid should generate approximately 25-km rook-neighbor centroid distances.
    q025, med, q975 = edges["centroid_distance_km"].quantile([0.025, 0.5, 0.975])
    if not (20 <= med <= 30 and 15 <= q025 <= 32 and 18 <= q975 <= 35):
        raise RuntimeError(f"Rook adjacency distance audit failed: q2.5={q025:.2f}, median={med:.2f}, q97.5={q975:.2f} km")
    pair_metrics = add_turnover(edges, incidence, traits)

    # Mark complete descriptive eligibility, including sparse ecoregions.
    primary_ids = set(cw.loc[cw["primary_assignment_eligible"], "step10b_cell_id"])
    unambig_ids = set(cw.loc[cw["sensitivity_unambiguous_eligible"], "step10b_cell_id"])
    formal_ids_primary = set(cw.loc[cw["primary_assignment_eligible"] & cw["dominant_ecoregion"].isin(formal), "step10b_cell_id"])
    formal_ids_unambig = set(cw.loc[cw["sensitivity_unambiguous_eligible"] & cw["dominant_ecoregion"].isin(formal), "step10b_cell_id"])
    pair_metrics["both_primary_eligible"] = pair_metrics["cell_i"].isin(primary_ids) & pair_metrics["cell_j"].isin(primary_ids)
    pair_metrics["both_unambiguous_eligible"] = pair_metrics["cell_i"].isin(unambig_ids) & pair_metrics["cell_j"].isin(unambig_ids)
    pair_metrics["both_primary_formal10"] = pair_metrics["cell_i"].isin(formal_ids_primary) & pair_metrics["cell_j"].isin(formal_ids_primary)
    pair_metrics["both_unambiguous_formal10"] = pair_metrics["cell_i"].isin(formal_ids_unambig) & pair_metrics["cell_j"].isin(formal_ids_unambig)
    pair_metrics.to_csv(out_dir / "10D_neighbor_pair_metrics.csv", index=False)

    rng = np.random.default_rng(args.seed)
    all_tests = []
    all_summaries = []
    all_nulls = {}
    for aset, col in [
        ("primary_formal10", "primary_assignment_eligible"),
        ("unambiguous_formal10", "sensitivity_unambiguous_eligible"),
    ]:
        t, s, n = run_permutations(pair_metrics, cw, aset, col, formal, args.permutations, rng)
        all_tests.append(t)
        all_summaries.append(s)
        all_nulls.update(n)
    tests = pd.concat(all_tests, ignore_index=True)
    summaries = pd.concat(all_summaries, ignore_index=True)
    tests.to_csv(out_dir / "10D_permutation_tests.csv", index=False)
    summaries.to_csv(out_dir / "10D_turnover_group_summaries.csv", index=False)

    # Compact null-distribution summary rather than a very large raw matrix.
    null_rows = []
    for key, arr in all_nulls.items():
        x = arr[np.isfinite(arr)]
        aset, assemblage, pair_filter, metric = key.split("::")
        null_rows.append({
            "analysis_set": aset, "assemblage": assemblage, "pair_filter": pair_filter, "metric": metric,
            "n": len(x), "mean": np.mean(x) if len(x) else np.nan, "sd": np.std(x, ddof=1) if len(x) > 1 else np.nan,
            "q005": np.quantile(x, .005) if len(x) else np.nan,
            "q025": np.quantile(x, .025) if len(x) else np.nan,
            "q50": np.quantile(x, .5) if len(x) else np.nan,
            "q975": np.quantile(x, .975) if len(x) else np.nan,
            "q995": np.quantile(x, .995) if len(x) else np.nan,
        })
    pd.DataFrame(null_rows).to_csv(out_dir / "10D_permutation_null_summaries.csv", index=False)

    bsum = boundary_pair_summary(pair_metrics, formal, cw)
    bsum.to_csv(out_dir / "10D_cross_boundary_ecoregion_pair_summary.csv", index=False)

    top = pair_metrics[pair_metrics["both_primary_formal10"]].sort_values("all_simpson", ascending=False).head(50)
    top.to_csv(out_dir / "10D_top_50_total_simpson_edges.csv", index=False)

    make_effect_plot(tests, "jaccard", pub_dir / "10D_jaccard_boundary_effect.png", pub_dir / "10D_jaccard_boundary_effect.pdf")
    make_effect_plot(tests, "simpson", pub_dir / "10D_simpson_boundary_effect.png", pub_dir / "10D_simpson_boundary_effect.pdf")
    make_edge_map(pair_metrics, cw, formal, pub_dir / "10D_cross_boundary_simpson_map.png", pub_dir / "10D_cross_boundary_simpson_map.pdf")

    summary = {
        "step": STEP,
        "audit_status": "PASS_FROZEN_COUNTS_REPRODUCED_AND_BOUNDARY_TURNOVER_COMPLETE",
        "seed": args.seed,
        "permutations": args.permutations,
        "project_root": str(project_root),
        "output_dir": str(out_dir),
        "occupied_cells": int(cw["step10b_cell_id"].nunique()),
        "total_genera": int(incidence["genus"].nunique()),
        "C3_exact_fraction": float(validation["C3_exact"].mean()),
        "N0_exact_fraction": float(validation["N0_exact"].mean()),
        "formal_ecoregions": len(formal),
        "formal_ecoregion_names": sorted(formal),
        "all_rook_pairs": int(len(pair_metrics)),
        "primary_all13_pairs": int(pair_metrics["both_primary_eligible"].sum()),
        "unambiguous_all13_pairs": int(pair_metrics["both_unambiguous_eligible"].sum()),
        "primary_formal10_pairs": int(pair_metrics["both_primary_formal10"].sum()),
        "unambiguous_formal10_pairs": int(pair_metrics["both_unambiguous_formal10"].sum()),
        "rook_distance_q025_km": float(q025),
        "rook_distance_median_km": float(med),
        "rook_distance_q975_km": float(q975),
        "primary_cross_boundary_pairs_all_metric": int(tests.loc[(tests["analysis_set"] == "primary_formal10") & (tests["assemblage"] == "all") & (tests["pair_filter"] == "union_positive") & (tests["metric"] == "simpson"), "cross_boundary_pairs"].iloc[0]),
        "primary_within_ecoregion_pairs_all_metric": int(tests.loc[(tests["analysis_set"] == "primary_formal10") & (tests["assemblage"] == "all") & (tests["pair_filter"] == "union_positive") & (tests["metric"] == "simpson"), "within_ecoregion_pairs"].iloc[0]),
        "next_step": "Interpret whether independent ecoregion boundaries align with total, C3, and N0 turnover; then decide whether a predefined published-break analysis is warranted as Step 10E.",
    }
    (out_dir / "10D_analysis_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    write_results_readme(out_dir, summary, tests)

    print("STEP 10D COMPLETE")
    print(f"AUDIT_STATUS={summary['audit_status']}")
    print(f"OUTPUT_DIR={out_dir}")
    print(f"OCCUPIED_CELLS={summary['occupied_cells']}")
    print(f"TOTAL_GENERA={summary['total_genera']}")
    print(f"ALL_ROOK_PAIRS={summary['all_rook_pairs']}")
    print(f"PRIMARY_FORMAL10_PAIRS={summary['primary_formal10_pairs']}")
    print(f"PRIMARY_CROSS_BOUNDARY_PAIRS={summary['primary_cross_boundary_pairs_all_metric']}")
    print(f"PRIMARY_WITHIN_ECOREGION_PAIRS={summary['primary_within_ecoregion_pairs_all_metric']}")
    print(f"PERMUTATIONS={args.permutations}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception as exc:
        # Print a full traceback so the user can return a useful error, while the runner archives diagnostics.
        print("STEP 10D FAILED", file=sys.stderr)
        print(str(exc), file=sys.stderr)
        traceback.print_exc()
        raise SystemExit(1)
