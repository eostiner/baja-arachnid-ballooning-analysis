#!/usr/bin/env python3
"""
Step 10C — independently defined ecoregion concordance for richness and C3 representation.

Primary spatial set: Step 10B cells with primary_assignment_eligible == TRUE.
Sensitivity set: Step 10B cells with sensitivity_unambiguous_eligible == TRUE.
Formal ecoregions: primary-set ecoregions represented by >= 5 occupied cells.
Equal-cell standardization: minimum formal-ecoregion cell count (expected 8),
5,000 iterations by default.

The script discovers and validates the final genus-by-cell incidence data and the
final D1–D4/N0 genus trait lookup. It will not proceed unless reconstructed C3 and
N0 cell counts closely reproduce the frozen Step 10B counts.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import sys
import traceback
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

try:
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
except Exception as exc:
    raise SystemExit(
        "Missing Python packages. Activate the Baja environment and install: "
        "pandas numpy matplotlib openpyxl\nOriginal error: %s" % exc
    )

STEP = "10C"
EXPECTED_TOTAL_GENERA = 267
EXPECTED_OCCUPIED_CELLS = 205
FORMAL_MIN_CELLS = 5
DEFAULT_ITERATIONS = 5000


def norm_name(x: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(x).lower())


def clean_str(x: Any) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip()


def md5_file(path: Path, block: int = 1024 * 1024) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        while True:
            b = f.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    x = s.astype(str).str.strip().str.lower()
    return x.isin({"true", "t", "1", "yes", "y"})


def read_delimited(path: Path, nrows: int | None = None) -> pd.DataFrame:
    encodings = ["utf-8-sig", "utf-8", "latin1"]
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab"} else None
    last: Exception | None = None
    for enc in encodings:
        try:
            if sep is None:
                return pd.read_csv(path, sep=None, engine="python", encoding=enc,
                                   nrows=nrows)
            return pd.read_csv(path, sep=sep, encoding=enc, nrows=nrows,
                               low_memory=False)
        except Exception as exc:
            last = exc
    raise RuntimeError(f"Could not read {path}: {last}")


def read_excel_candidates(path: Path) -> list[tuple[str, pd.DataFrame]]:
    try:
        xls = pd.ExcelFile(path, engine="openpyxl")
    except Exception as exc:
        return [("<read_error>", pd.DataFrame({"_error": [str(exc)]}))]
    out: list[tuple[str, pd.DataFrame]] = []
    for sheet in xls.sheet_names[:20]:
        try:
            d = pd.read_excel(path, sheet_name=sheet, engine="openpyxl")
            out.append((sheet, d))
        except Exception:
            continue
    return out


def candidate_tabular_files(project_root: Path) -> list[Path]:
    allowed = {".csv", ".tsv", ".tab", ".xlsx", ".xls"}
    files: list[Path] = []
    for root, dirs, names in os.walk(project_root):
        # Avoid recursing into source-download caches and this step's outputs.
        dirs[:] = [d for d in dirs if d not in {".git", "__pycache__", "00_source_download"}]
        if "10C_equal_cell_ecoregion_richness" in root:
            continue
        for name in names:
            p = Path(root) / name
            if p.suffix.lower() not in allowed:
                continue
            try:
                size = p.stat().st_size
            except OSError:
                continue
            if size <= 0 or size > 600_000_000:
                continue
            files.append(p)
    # Known reaudited workbook, if it is stored beside the project rather than inside it.
    desktop = project_root.parent
    for known in [
        "USE_GOOD_BalloonID_Baja_Arachnid_Species_Long_REAUDITED.xlsx",
        "USE_GOOD_BalloonID_Baja_Arachnid_Species_Long.xlsx",
    ]:
        p = desktop / known
        if p.exists():
            files.append(p)
    return sorted(set(files))


def find_cell_column(df: pd.DataFrame, crosswalk_ids: set[str]) -> tuple[str | None, float]:
    preferred = {
        "grid25kmid", "cellid", "gridid", "gridcellid", "step10bcellid",
        "siteid", "cell", "gridcell"
    }
    best_col: str | None = None
    best_score = -1.0
    for col in df.columns:
        vals = df[col].dropna().astype(str).str.strip()
        if vals.empty:
            continue
        unique = set(vals.unique())
        overlap = len(unique & crosswalk_ids)
        recall = overlap / max(1, len(crosswalk_ids))
        precision = overlap / max(1, len(unique))
        score = 100 * recall + 20 * precision
        if norm_name(col) in preferred:
            score += 25
        if score > best_score:
            best_score = score
            best_col = str(col)
    recall = max(0.0, (best_score - (25 if best_col and norm_name(best_col) in preferred else 0)) / 100)
    return best_col, recall


def find_genus_column(df: pd.DataFrame) -> str | None:
    exact = [
        "genus", "acceptedgenus", "finalgenus", "genusclean", "genusname",
        "canonicalgenus", "accepted_genus", "final_genus"
    ]
    normalized = {norm_name(c): str(c) for c in df.columns}
    for name in exact:
        if norm_name(name) in normalized:
            return normalized[norm_name(name)]
    # Conservative fallback: a column whose name contains genus but not richness/count.
    for col in df.columns:
        n = norm_name(col)
        if "genus" in n and not any(z in n for z in ("count", "richness", "number", "nunique")):
            return str(col)
    return None


def detect_presence_column(df: pd.DataFrame) -> str | None:
    normalized = {norm_name(c): str(c) for c in df.columns}
    for key in ["presence", "present", "incidence", "occupied", "pa"]:
        if key in normalized:
            return normalized[key]
    return None


def extract_long_incidence(df: pd.DataFrame, cell_col: str, genus_col: str,
                           crosswalk_ids: set[str]) -> pd.DataFrame:
    d = df[[cell_col, genus_col] + ([detect_presence_column(df)] if detect_presence_column(df) else [])].copy()
    d.columns = ["cell_id", "genus"] + (["presence"] if d.shape[1] == 3 else [])
    d["cell_id"] = d["cell_id"].map(clean_str)
    d["genus"] = d["genus"].map(clean_str)
    if "presence" in d.columns:
        p = pd.to_numeric(d["presence"], errors="coerce")
        if p.notna().mean() > 0.8:
            d = d[p.fillna(0) > 0]
        else:
            d = d[bool_series(d["presence"])]
    d = d[d["cell_id"].isin(crosswalk_ids)]
    d = d[(d["genus"] != "") & (~d["genus"].str.lower().isin({"na", "nan", "unknown", "unidentified"}))]
    return d[["cell_id", "genus"]].drop_duplicates().reset_index(drop=True)


def extract_wide_incidence(df: pd.DataFrame, cell_col: str,
                           crosswalk_ids: set[str]) -> pd.DataFrame | None:
    if len(df) < 150 or len(df) > 500:
        return None
    ids = df[cell_col].astype(str).str.strip()
    if ids.isin(crosswalk_ids).mean() < 0.7:
        return None
    genus_cols: list[str] = []
    for col in df.columns:
        if str(col) == cell_col:
            continue
        x = pd.to_numeric(df[col], errors="coerce")
        good = x.dropna()
        if len(good) < 0.9 * len(df):
            continue
        vals = set(np.unique(good))
        if vals.issubset({0, 1}) and len(vals) >= 1:
            n = norm_name(col)
            if not any(k in n for k in ["validfraction", "eligible", "complete", "ge5", "ge10", "candidate"]):
                genus_cols.append(str(col))
    if len(genus_cols) < 100:
        return None
    rows: list[tuple[str, str]] = []
    for _, row in df[[cell_col] + genus_cols].iterrows():
        cid = clean_str(row[cell_col])
        if cid not in crosswalk_ids:
            continue
        for genus in genus_cols:
            try:
                present = float(row[genus]) > 0
            except Exception:
                present = False
            if present:
                rows.append((cid, genus))
    return pd.DataFrame(rows, columns=["cell_id", "genus"]).drop_duplicates()


def discover_incidence(files: list[Path], crosswalk: pd.DataFrame,
                       discovery_rows: list[dict[str, Any]]) -> tuple[pd.DataFrame, Path, str]:
    ids = set(crosswalk["step10b_cell_id"].astype(str))
    candidates: list[tuple[float, pd.DataFrame, Path, str]] = []
    priority_terms = ["incidence", "genus_cell", "genusbycell", "grid25km", "community", "presence"]
    for path in files:
        if path.suffix.lower() in {".xlsx", ".xls"}:
            continue
        try:
            sample = read_delimited(path, nrows=2500)
        except Exception as exc:
            discovery_rows.append({"role": "incidence", "file": str(path), "status": "read_error", "detail": str(exc)[:300]})
            continue
        if sample.empty or sample.shape[1] < 2:
            continue
        cell_col, _ = find_cell_column(sample, ids)
        if cell_col is None:
            continue
        genus_col = find_genus_column(sample)
        mode = "long" if genus_col else "wide"
        path_score = sum(12 for t in priority_terms if t in norm_name(path.name))
        if "08grid25kmincidence" in norm_name(str(path.parent)):
            path_score += 80
        try:
            full = read_delimited(path)
            if genus_col and genus_col in full.columns and cell_col in full.columns:
                inc = extract_long_incidence(full, cell_col, genus_col, ids)
            elif cell_col in full.columns:
                inc = extract_wide_incidence(full, cell_col, ids)
                if inc is None:
                    continue
            else:
                continue
        except Exception as exc:
            discovery_rows.append({"role": "incidence", "file": str(path), "status": "parse_error", "detail": str(exc)[:300]})
            continue
        n_cells = inc["cell_id"].nunique()
        n_genera = inc["genus"].nunique()
        score = path_score - abs(n_cells - EXPECTED_OCCUPIED_CELLS) * 3 - abs(n_genera - EXPECTED_TOTAL_GENERA) * 2
        if 190 <= n_cells <= 210:
            score += 150
        if 250 <= n_genera <= 285:
            score += 150
        discovery_rows.append({
            "role": "incidence", "file": str(path), "status": "candidate", "detail": mode,
            "n_rows_deduplicated": len(inc), "n_cells": n_cells, "n_genera": n_genera, "score": score,
            "cell_column": cell_col, "genus_column": genus_col or "<wide columns>"
        })
        if n_cells >= 180 and n_genera >= 150:
            candidates.append((score, inc, path, mode))
    if not candidates:
        raise RuntimeError("Could not find a usable genus-by-cell incidence table. See 10C_discovery_candidates.csv.")
    candidates.sort(key=lambda x: x[0], reverse=True)
    score, inc, path, mode = candidates[0]
    if inc["cell_id"].nunique() != EXPECTED_OCCUPIED_CELLS:
        raise RuntimeError(
            f"Best incidence file has {inc['cell_id'].nunique()} occupied cells, expected {EXPECTED_OCCUPIED_CELLS}: {path}"
        )
    if inc["genus"].nunique() != EXPECTED_TOTAL_GENERA:
        raise RuntimeError(
            f"Best incidence file has {inc['genus'].nunique()} genera, expected {EXPECTED_TOTAL_GENERA}: {path}"
        )
    return inc, path, mode


TOKEN_RE = re.compile(r"(?<![A-Z0-9])(D1|D2|D3|D4|N0|C3)(?![A-Z0-9])", re.I)


def parse_evidence_value(value: Any) -> str | None:
    if pd.isna(value):
        return None
    s = str(value).strip().upper()
    if not s:
        return None
    toks = set(t.upper() for t in TOKEN_RE.findall(s))
    # Explicit single class or tier.
    if toks == {"N0"}:
        return "N0"
    if toks == {"C3"}:
        return "C3"
    if len(toks) == 1:
        t = next(iter(toks))
        if t in {"D1", "D2", "D3", "D4"}:
            return t
    # Common exact text labels not carrying explicit codes.
    n = norm_name(s)
    if n in {"nonballooning", "fixednonballooning", "reference_nonballooning", "noballooning"}:
        return "N0"
    return None


def evidence_columns(df: pd.DataFrame) -> list[tuple[float, str, float]]:
    out: list[tuple[float, str, float]] = []
    for col in df.columns:
        vals = df[col].dropna()
        if vals.empty:
            continue
        parsed = vals.map(parse_evidence_value)
        frac = parsed.notna().mean()
        unique_classes = set(parsed.dropna())
        n = norm_name(col)
        name_bonus = 0
        if any(k in n for k in ["evidence", "tier", "category", "class", "decision", "designation", "balloon"]):
            name_bonus += 25
        if n in {"dlevel", "evidencecategory", "finalevidencecategory", "traitclass", "primaryclass"}:
            name_bonus += 40
        if frac >= 0.25 and len(unique_classes) >= 2:
            out.append((100 * frac + name_bonus + 5 * len(unique_classes), str(col), frac))
    return sorted(out, reverse=True)


def trait_tables_from_file(path: Path) -> Iterable[tuple[str, pd.DataFrame]]:
    if path.suffix.lower() in {".xlsx", ".xls"}:
        yield from read_excel_candidates(path)
    else:
        try:
            yield "<table>", read_delimited(path)
        except Exception:
            return


def evaluate_trait_table(df: pd.DataFrame, source: str, sheet: str,
                         incidence: pd.DataFrame, crosswalk: pd.DataFrame,
                         discovery_rows: list[dict[str, Any]]) -> tuple[float, pd.DataFrame, dict[str, Any]] | None:
    if df.empty:
        return None
    genus_col = find_genus_column(df)
    if genus_col is None:
        return None
    ev_cols = evidence_columns(df)
    if not ev_cols:
        return None
    best: tuple[float, pd.DataFrame, dict[str, Any]] | None = None
    for ev_score, ev_col, parse_frac in ev_cols[:5]:
        d = df[[genus_col, ev_col]].copy()
        d.columns = ["genus", "raw_evidence"]
        d["genus"] = d["genus"].map(clean_str)
        d["evidence_class"] = d["raw_evidence"].map(parse_evidence_value)
        d = d[(d["genus"] != "") & d["evidence_class"].notna()]
        if d.empty:
            continue
        # Collapse duplicates only when all classified rows for a genus agree.
        grouped = d.groupby("genus")["evidence_class"].agg(lambda x: sorted(set(x)))
        conflict_genera = grouped[grouped.map(len) > 1].index.tolist()
        clean = grouped[grouped.map(len) == 1].map(lambda x: x[0]).reset_index()
        clean["analysis_class"] = np.where(clean["evidence_class"].isin(["D1", "D2", "D3", "C3"]), "C3",
                                             np.where(clean["evidence_class"] == "N0", "N0", "D4_excluded"))
        n_overlap = clean["genus"].isin(set(incidence["genus"])).sum()
        # Validate against frozen cell counts.
        m = incidence.merge(clean[["genus", "analysis_class"]], on="genus", how="left")
        counts = (m[m["analysis_class"].isin(["C3", "N0"])]
                  .groupby(["cell_id", "analysis_class"])["genus"].nunique()
                  .unstack(fill_value=0).reset_index())
        for col in ["C3", "N0"]:
            if col not in counts.columns:
                counts[col] = 0
        val = crosswalk[["step10b_cell_id", "C3_positive_genera", "N0_reference_genera"]].merge(
            counts, left_on="step10b_cell_id", right_on="cell_id", how="left"
        ).fillna({"C3": 0, "N0": 0})
        c3_exact = float((val["C3"].astype(int) == val["C3_positive_genera"].astype(int)).mean())
        n0_exact = float((val["N0"].astype(int) == val["N0_reference_genera"].astype(int)).mean())
        mae = float((val["C3"].astype(float) - val["C3_positive_genera"].astype(float)).abs().mean() +
                    (val["N0"].astype(float) - val["N0_reference_genera"].astype(float)).abs().mean())
        score = 500 * (c3_exact + n0_exact) + n_overlap + ev_score - 20 * len(conflict_genera) - 30 * mae
        meta = {
            "source": source, "sheet": sheet, "genus_column": genus_col, "evidence_column": ev_col,
            "parse_fraction": parse_frac, "classified_genera": len(clean), "incidence_genus_overlap": int(n_overlap),
            "conflict_genera": len(conflict_genera), "c3_cell_exact_fraction": c3_exact,
            "n0_cell_exact_fraction": n0_exact, "combined_cell_mae": mae, "score": score,
        }
        discovery_rows.append({"role": "trait", "file": source, "status": "candidate", "detail": sheet, **meta})
        if best is None or score > best[0]:
            best = (score, clean, meta)
    return best


def discover_traits(files: list[Path], incidence: pd.DataFrame, crosswalk: pd.DataFrame,
                    incidence_source: Path, discovery_rows: list[dict[str, Any]]) -> tuple[pd.DataFrame, dict[str, Any]]:
    # First evaluate the incidence source itself in case traits were carried through.
    priority: list[Path] = [incidence_source]
    keywords = ["trait", "balloon", "reaudit", "lookup", "evidence", "designation", "merge", "genus"]
    for p in files:
        if p == incidence_source:
            continue
        n = norm_name(p.name)
        parent = norm_name(str(p.parent))
        if any(k in n or k in parent for k in keywords):
            priority.append(p)
    # Deduplicate and cap obviously irrelevant large workbook search.
    seen: set[Path] = set()
    best: tuple[float, pd.DataFrame, dict[str, Any]] | None = None
    for path in priority:
        if path in seen:
            continue
        seen.add(path)
        try:
            for sheet, df in trait_tables_from_file(path):
                result = evaluate_trait_table(df, str(path), sheet, incidence, crosswalk, discovery_rows)
                if result is not None and (best is None or result[0] > best[0]):
                    best = result
        except Exception as exc:
            discovery_rows.append({"role": "trait", "file": str(path), "status": "read_error", "detail": str(exc)[:300]})
    if best is None:
        raise RuntimeError("Could not identify a genus trait table with explicit D1–D4/N0 or C3/N0 classes.")
    score, lookup, meta = best
    if meta["c3_cell_exact_fraction"] < 0.98 or meta["n0_cell_exact_fraction"] < 0.98:
        raise RuntimeError(
            "Best trait lookup failed frozen-count validation: "
            f"C3 exact={meta['c3_cell_exact_fraction']:.3f}, N0 exact={meta['n0_cell_exact_fraction']:.3f}, "
            f"source={meta['source']} sheet={meta['sheet']} column={meta['evidence_column']}. "
            "See 10C_discovery_candidates.csv."
        )
    return lookup, meta


def union_richness(sampled_cells: np.ndarray, cell_sets: dict[str, set[str]]) -> int:
    u: set[str] = set()
    for c in sampled_cells:
        u.update(cell_sets.get(str(c), set()))
    return len(u)


def summarize_iterations(iter_df: pd.DataFrame) -> pd.DataFrame:
    metrics = ["total_richness", "C3_richness", "N0_richness", "C3_representation"]
    rows: list[dict[str, Any]] = []
    for (aset, eco), g in iter_df.groupby(["analysis_set", "ecoregion"]):
        for metric in metrics:
            x = g[metric].dropna().to_numpy(float)
            rows.append({
                "analysis_set": aset, "ecoregion": eco, "metric": metric,
                "iterations": len(x), "equal_cells": int(g["equal_cells"].iloc[0]),
                "mean": float(np.mean(x)) if len(x) else np.nan,
                "median": float(np.median(x)) if len(x) else np.nan,
                "ci_low_2.5": float(np.quantile(x, 0.025)) if len(x) else np.nan,
                "ci_high_97.5": float(np.quantile(x, 0.975)) if len(x) else np.nan,
                "sd": float(np.std(x, ddof=1)) if len(x) > 1 else np.nan,
            })
    return pd.DataFrame(rows)


def pairwise_contrasts(iter_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    metrics = ["total_richness", "C3_richness", "N0_richness", "C3_representation"]
    for aset, ad in iter_df.groupby("analysis_set"):
        ecos = sorted(ad["ecoregion"].unique())
        by_eco = {e: ad[ad["ecoregion"] == e].sort_values("iteration") for e in ecos}
        for metric in metrics:
            for i, e1 in enumerate(ecos):
                for e2 in ecos[i + 1:]:
                    x = by_eco[e1][metric].to_numpy(float)
                    y = by_eco[e2][metric].to_numpy(float)
                    d = x - y
                    d = d[np.isfinite(d)]
                    rows.append({
                        "analysis_set": aset, "metric": metric,
                        "ecoregion_1": e1, "ecoregion_2": e2,
                        "mean_difference_1_minus_2": float(np.mean(d)),
                        "ci_low_2.5": float(np.quantile(d, 0.025)),
                        "ci_high_97.5": float(np.quantile(d, 0.975)),
                        "probability_ecoregion_1_greater": float(np.mean(d > 0)),
                        "interval_excludes_zero": bool(np.quantile(d, 0.025) > 0 or np.quantile(d, 0.975) < 0),
                    })
    return pd.DataFrame(rows)


def eta_squared(values: np.ndarray, groups: np.ndarray) -> float:
    good = np.isfinite(values)
    y = values[good]
    g = groups[good]
    if len(y) < 3:
        return np.nan
    grand = y.mean()
    ss_total = np.sum((y - grand) ** 2)
    if ss_total <= 0:
        return 0.0
    ss_between = 0.0
    for label in np.unique(g):
        z = y[g == label]
        if len(z):
            ss_between += len(z) * (z.mean() - grand) ** 2
    return float(ss_between / ss_total)


def permutation_test(values: np.ndarray, groups: np.ndarray, strata: np.ndarray | None,
                     rng: np.random.Generator, n_perm: int) -> tuple[float, float, int]:
    good = np.isfinite(values) & pd.notna(groups)
    if strata is not None:
        good &= pd.notna(strata)
    y = values[good]
    g = groups[good].copy()
    s = strata[good] if strata is not None else None
    obs = eta_squared(y, g)
    ge = 0
    for _ in range(n_perm):
        gp = g.copy()
        if s is None:
            rng.shuffle(gp)
        else:
            for st in np.unique(s):
                idx = np.flatnonzero(s == st)
                gp[idx] = rng.permutation(gp[idx])
        stat = eta_squared(y, gp)
        if stat >= obs - 1e-15:
            ge += 1
    return obs, (ge + 1) / (n_perm + 1), len(y)


def safe_label(s: str) -> str:
    return s.replace(" de ", " de\n").replace(" del ", " del\n")


def plot_richness(summary: pd.DataFrame, out_dir: Path) -> None:
    primary = summary[(summary["analysis_set"] == "primary") &
                      (summary["metric"].isin(["total_richness", "C3_richness", "N0_richness"]))].copy()
    order_df = primary[primary["metric"] == "total_richness"].sort_values("mean")
    order = order_df["ecoregion"].tolist()
    metrics = [("total_richness", "All genera"), ("C3_richness", "C3 genera"), ("N0_richness", "N0 genera")]
    fig, axes = plt.subplots(1, 3, figsize=(15, 7), sharey=True)
    for ax, (metric, title) in zip(axes, metrics):
        d = primary[primary["metric"] == metric].set_index("ecoregion").reindex(order).reset_index()
        y = np.arange(len(d))
        ax.errorbar(d["mean"], y,
                    xerr=np.vstack([d["mean"] - d["ci_low_2.5"], d["ci_high_97.5"] - d["mean"]]),
                    fmt="o", capsize=3)
        ax.set_title(title)
        ax.set_xlabel("Expected pooled genus richness\n(8 occupied cells)")
        ax.grid(axis="x", alpha=0.25)
        ax.set_yticks(y)
        ax.set_yticklabels([safe_label(x) for x in d["ecoregion"]])
    fig.suptitle("Equal-cell richness across independently mapped Baja ecoregions")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_dir / "10C_equal_cell_ecoregion_richness.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "10C_equal_cell_ecoregion_richness.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_representation(summary: pd.DataFrame, out_dir: Path) -> None:
    d = summary[(summary["analysis_set"] == "primary") & (summary["metric"] == "C3_representation")].copy()
    d = d.sort_values("mean")
    y = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.errorbar(d["mean"], y,
                xerr=np.vstack([d["mean"] - d["ci_low_2.5"], d["ci_high_97.5"] - d["mean"]]),
                fmt="o", capsize=3)
    ax.set_yticks(y)
    ax.set_yticklabels([safe_label(x) for x in d["ecoregion"]])
    ax.set_xlabel("C3 / (C3 + N0) pooled richness\n(8 occupied cells)")
    ax.set_title("C3 representation across independently mapped Baja ecoregions")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_dir / "10C_equal_cell_C3_representation.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "10C_equal_cell_C3_representation.pdf", bbox_inches="tight")
    plt.close(fig)


def write_readme(out_dir: Path, summary_json: dict[str, Any]) -> None:
    text = f"""STEP 10C RESULTS — EQUAL-CELL ECOREGION RICHNESS AND C3 REPRESENTATION

STATUS
{summary_json['audit_status']}

DESIGN
- Independent spatial units: validated González-Abraham et al. mainland ecoregions from Steps 10A–10B.
- Primary cells: {summary_json['primary_cells']}.
- Unambiguous sensitivity cells: {summary_json['unambiguous_cells']}.
- Formal ecoregions: {summary_json['formal_ecoregions']} with at least {FORMAL_MIN_CELLS} primary occupied cells.
- Equal-cell sample: {summary_json['equal_cells']} occupied cells per ecoregion.
- Monte Carlo iterations: {summary_json['iterations']}.
- Primary ballooning definition: C3 = D1 + D2 + D3.
- Fixed non-ballooning reference: N0.

INTERPRETATION RULES
1. Use 10C_equal_cell_summary.csv for expected pooled richness and 95% Monte Carlo intervals.
2. Use 10C_pairwise_contrasts.csv for paired iteration-level contrasts; do not treat every pairwise interval as an independent hypothesis test.
3. Use 10C_global_permutation_tests.csv to assess overall ecoregion structure. The latitude-stratified permutation is the more conservative test because it asks whether ecoregion identity explains variation beyond the five 2-degree latitude bands.
4. C3 representation is C3/(C3+N0); D4 is excluded from this denominator.
5. The denominator >=5 and >=10 analyses are cell-level sensitivities and may have sparse coverage in some ecoregions. They do not replace the equal-cell pooled estimate.
6. Step 10C tests richness and representation. Cross-boundary turnover is reserved for Step 10D.

VALIDATION
- Incidence source: {summary_json['incidence_source']}
- Trait source: {summary_json['trait_source']}
- C3 frozen-count exact match: {summary_json['c3_exact_fraction']:.3f}
- N0 frozen-count exact match: {summary_json['n0_exact_fraction']:.3f}

Do not interpret geographic concordance as causal environmental control. Ecoregions are an independently defined biogeographic hypothesis, and spatial autocorrelation remains relevant.
"""
    (out_dir / "README_RESULTS_FIRST.txt").write_text(text, encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("project_root", nargs="?", default="~/Desktop/Baja_Ballooning_Pipeline")
    ap.add_argument("seed", nargs="?", type=int, default=20260715)
    ap.add_argument("iterations", nargs="?", type=int, default=DEFAULT_ITERATIONS)
    args = ap.parse_args()

    project_root = Path(os.path.expanduser(args.project_root)).resolve()
    pipeline_root = project_root / "04_analysis" / "C3_pipeline_rebuild"
    step10_root = pipeline_root / "09_C3_biogeographic_concordance"
    step10b_dir = step10_root / "10B_cell_ecoregion_crosswalk"
    out_dir = step10_root / "10C_equal_cell_ecoregion_richness"
    pub_dir = out_dir / "publication_outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    pub_dir.mkdir(parents=True, exist_ok=True)

    crosswalk_path = step10b_dir / "10B_cell_ecoregion_crosswalk.csv"
    if not crosswalk_path.exists():
        raise FileNotFoundError(f"Missing Step 10B crosswalk: {crosswalk_path}")
    crosswalk = pd.read_csv(crosswalk_path, low_memory=False)
    required = {
        "step10b_cell_id", "dominant_ecoregion", "primary_assignment_eligible",
        "sensitivity_unambiguous_eligible", "latitude_band", "C3_positive_genera",
        "N0_reference_genera", "C3_classified_denominator"
    }
    missing = required - set(crosswalk.columns)
    if missing:
        raise RuntimeError(f"Step 10B crosswalk is missing columns: {sorted(missing)}")
    if len(crosswalk) != EXPECTED_OCCUPIED_CELLS:
        raise RuntimeError(f"Step 10B crosswalk has {len(crosswalk)} rows, expected {EXPECTED_OCCUPIED_CELLS}.")
    crosswalk["step10b_cell_id"] = crosswalk["step10b_cell_id"].astype(str).str.strip()
    crosswalk["primary_assignment_eligible"] = bool_series(crosswalk["primary_assignment_eligible"])
    crosswalk["sensitivity_unambiguous_eligible"] = bool_series(crosswalk["sensitivity_unambiguous_eligible"])

    discovery_rows: list[dict[str, Any]] = []
    files = candidate_tabular_files(project_root)
    try:
        incidence, incidence_source, incidence_mode = discover_incidence(files, crosswalk, discovery_rows)
    except Exception:
        pd.DataFrame(discovery_rows).to_csv(out_dir / "10C_discovery_candidates.csv", index=False)
        raise
    # Save discovery before trait selection so failures remain diagnosable.
    pd.DataFrame(discovery_rows).to_csv(out_dir / "10C_discovery_candidates.csv", index=False)

    try:
        lookup, trait_meta = discover_traits(files, incidence, crosswalk, incidence_source, discovery_rows)
    except Exception:
        pd.DataFrame(discovery_rows).to_csv(out_dir / "10C_discovery_candidates.csv", index=False)
        raise
    pd.DataFrame(discovery_rows).to_csv(out_dir / "10C_discovery_candidates.csv", index=False)
    lookup.to_csv(out_dir / "10C_genus_trait_lookup_used.csv", index=False)

    # Construct and validate per-cell response values.
    inc = incidence.merge(lookup[["genus", "analysis_class", "evidence_class"]], on="genus", how="left")
    total_counts = inc.groupby("cell_id")["genus"].nunique().rename("total_genera")
    class_counts = (inc[inc["analysis_class"].isin(["C3", "N0"])]
                    .groupby(["cell_id", "analysis_class"])["genus"].nunique()
                    .unstack(fill_value=0))
    for col in ["C3", "N0"]:
        if col not in class_counts.columns:
            class_counts[col] = 0
    response = crosswalk.merge(total_counts, left_on="step10b_cell_id", right_index=True, how="left")
    response = response.merge(class_counts[["C3", "N0"]], left_on="step10b_cell_id", right_index=True, how="left")
    response[["total_genera", "C3", "N0"]] = response[["total_genera", "C3", "N0"]].fillna(0).astype(int)
    response["C3_count_matches_frozen"] = response["C3"] == response["C3_positive_genera"].astype(int)
    response["N0_count_matches_frozen"] = response["N0"] == response["N0_reference_genera"].astype(int)
    response["C3_N0_denominator"] = response["C3"] + response["N0"]
    response["C3_representation"] = np.where(response["C3_N0_denominator"] > 0,
                                               response["C3"] / response["C3_N0_denominator"], np.nan)
    response.to_csv(out_dir / "10C_cell_response_validation.csv", index=False)
    c3_exact = float(response["C3_count_matches_frozen"].mean())
    n0_exact = float(response["N0_count_matches_frozen"].mean())
    if c3_exact < 0.98 or n0_exact < 0.98:
        raise RuntimeError(f"Frozen count validation failed after source selection: C3={c3_exact:.3f}, N0={n0_exact:.3f}")

    primary = response[response["primary_assignment_eligible"] & response["dominant_ecoregion"].notna()].copy()
    sensitivity = response[response["sensitivity_unambiguous_eligible"] & response["dominant_ecoregion"].notna()].copy()
    primary_counts = primary.groupby("dominant_ecoregion").size().sort_values(ascending=False)
    formal_ecos = sorted(primary_counts[primary_counts >= FORMAL_MIN_CELLS].index.tolist())
    if len(formal_ecos) != 10:
        raise RuntimeError(f"Expected 10 formal ecoregions with >=5 cells; found {len(formal_ecos)}: {formal_ecos}")
    equal_n = int(primary_counts.loc[formal_ecos].min())
    sens_counts = sensitivity.groupby("dominant_ecoregion").size().reindex(formal_ecos)
    if sens_counts.min() < equal_n:
        raise RuntimeError(
            f"Unambiguous sensitivity set cannot preserve the locked equal-cell size {equal_n}; minimum is {sens_counts.min()}."
        )
    if equal_n != 8:
        raise RuntimeError(f"Expected locked equal-cell size 8, found {equal_n}.")

    sample_sizes = pd.DataFrame({
        "ecoregion": sorted(response["dominant_ecoregion"].dropna().unique()),
    })
    sample_sizes["primary_cells"] = sample_sizes["ecoregion"].map(primary_counts).fillna(0).astype(int)
    sample_sizes["unambiguous_cells"] = sample_sizes["ecoregion"].map(sensitivity.groupby("dominant_ecoregion").size()).fillna(0).astype(int)
    sample_sizes["formal_comparison"] = sample_sizes["ecoregion"].isin(formal_ecos)
    for threshold in [1, 5, 10]:
        counts = primary[primary["C3_N0_denominator"] >= threshold].groupby("dominant_ecoregion").size()
        sample_sizes[f"primary_cells_denominator_ge{threshold}"] = sample_sizes["ecoregion"].map(counts).fillna(0).astype(int)
    sample_sizes.to_csv(out_dir / "10C_ecoregion_sample_sizes.csv", index=False)

    # Cell sets for pooled richness.
    cell_total = {cid: set(g["genus"]) for cid, g in inc.groupby("cell_id")}
    cell_c3 = {cid: set(g.loc[g["analysis_class"] == "C3", "genus"]) for cid, g in inc.groupby("cell_id")}
    cell_n0 = {cid: set(g.loc[g["analysis_class"] == "N0", "genus"]) for cid, g in inc.groupby("cell_id")}

    rng = np.random.default_rng(args.seed)
    iter_rows: list[dict[str, Any]] = []
    for aset, dset in [("primary", primary), ("unambiguous_sensitivity", sensitivity)]:
        pools = {eco: dset.loc[dset["dominant_ecoregion"] == eco, "step10b_cell_id"].astype(str).to_numpy()
                 for eco in formal_ecos}
        for iteration in range(1, args.iterations + 1):
            for eco in formal_ecos:
                sampled = rng.choice(pools[eco], size=equal_n, replace=False)
                rt = union_richness(sampled, cell_total)
                rc3 = union_richness(sampled, cell_c3)
                rn0 = union_richness(sampled, cell_n0)
                denom = rc3 + rn0
                iter_rows.append({
                    "analysis_set": aset, "iteration": iteration, "ecoregion": eco,
                    "equal_cells": equal_n, "total_richness": rt, "C3_richness": rc3,
                    "N0_richness": rn0, "C3_N0_denominator": denom,
                    "C3_representation": rc3 / denom if denom else np.nan,
                })
    iter_df = pd.DataFrame(iter_rows)
    iter_df.to_csv(out_dir / "10C_equal_cell_iteration_results.csv", index=False)
    summary = summarize_iterations(iter_df)
    summary.to_csv(out_dir / "10C_equal_cell_summary.csv", index=False)
    contrasts = pairwise_contrasts(iter_df)
    contrasts.to_csv(out_dir / "10C_pairwise_contrasts.csv", index=False)

    # Cell-level descriptive summaries and global ecoregion label-permutation tests.
    cell_summary_rows: list[dict[str, Any]] = []
    perm_rows: list[dict[str, Any]] = []
    perm_rng = np.random.default_rng(args.seed + 1009)
    n_perm = args.iterations
    for aset, dset in [("primary", primary), ("unambiguous_sensitivity", sensitivity)]:
        dformal = dset[dset["dominant_ecoregion"].isin(formal_ecos)].copy()
        response_specs = [
            ("total_cell_richness", "total_genera", 1),
            ("C3_cell_richness", "C3", 1),
            ("N0_cell_richness", "N0", 1),
            ("C3_representation_denominator_ge1", "C3_representation", 1),
            ("C3_representation_denominator_ge5", "C3_representation", 5),
            ("C3_representation_denominator_ge10", "C3_representation", 10),
        ]
        for response_name, col, threshold in response_specs:
            z = dformal[dformal["C3_N0_denominator"] >= threshold].copy() if "representation" in response_name else dformal.copy()
            for eco, g in z.groupby("dominant_ecoregion"):
                x = g[col].dropna().to_numpy(float)
                cell_summary_rows.append({
                    "analysis_set": aset, "response": response_name, "ecoregion": eco,
                    "n_cells": len(x), "mean": float(np.mean(x)) if len(x) else np.nan,
                    "median": float(np.median(x)) if len(x) else np.nan,
                    "sd": float(np.std(x, ddof=1)) if len(x) > 1 else np.nan,
                    "min": float(np.min(x)) if len(x) else np.nan,
                    "max": float(np.max(x)) if len(x) else np.nan,
                })
            values = z[col].to_numpy(float)
            groups = z["dominant_ecoregion"].astype(str).to_numpy()
            strata = z["latitude_band"].astype(str).to_numpy()
            obs_u, p_u, n_u = permutation_test(values, groups, None, perm_rng, n_perm)
            obs_s, p_s, n_s = permutation_test(values, groups, strata, perm_rng, n_perm)
            perm_rows.extend([
                {"analysis_set": aset, "response": response_name, "permutation_scheme": "unrestricted",
                 "n_cells": n_u, "eta_squared": obs_u, "permutations": n_perm, "p_value": p_u},
                {"analysis_set": aset, "response": response_name, "permutation_scheme": "within_latitude_band",
                 "n_cells": n_s, "eta_squared": obs_s, "permutations": n_perm, "p_value": p_s},
            ])
    pd.DataFrame(cell_summary_rows).to_csv(out_dir / "10C_cell_level_ecoregion_summary.csv", index=False)
    pd.DataFrame(perm_rows).to_csv(out_dir / "10C_global_permutation_tests.csv", index=False)

    # Manifest and figures.
    trait_source = Path(trait_meta["source"])
    manifest = pd.DataFrame([
        {"role": "step10b_crosswalk", "file": str(crosswalk_path), "md5": md5_file(crosswalk_path)},
        {"role": "genus_cell_incidence", "file": str(incidence_source), "md5": md5_file(incidence_source)},
        {"role": "genus_trait_lookup", "file": str(trait_source), "md5": md5_file(trait_source)},
    ])
    manifest.to_csv(out_dir / "10C_input_manifest.csv", index=False)
    plot_richness(summary, pub_dir)
    plot_representation(summary, pub_dir)

    summary_json = {
        "step": STEP,
        "audit_status": "PASS_FROZEN_COUNTS_REPRODUCED_AND_EQUAL_CELL_ANALYSIS_COMPLETE",
        "seed": args.seed,
        "iterations": args.iterations,
        "project_root": str(project_root),
        "output_dir": str(out_dir),
        "incidence_source": str(incidence_source),
        "incidence_mode": incidence_mode,
        "incidence_rows": int(len(incidence)),
        "occupied_cells": int(incidence["cell_id"].nunique()),
        "total_genera": int(incidence["genus"].nunique()),
        "trait_source": trait_meta["source"],
        "trait_sheet": trait_meta["sheet"],
        "trait_evidence_column": trait_meta["evidence_column"],
        "c3_exact_fraction": c3_exact,
        "n0_exact_fraction": n0_exact,
        "primary_cells": int(len(primary)),
        "unambiguous_cells": int(len(sensitivity)),
        "formal_ecoregions": int(len(formal_ecos)),
        "formal_ecoregion_names": formal_ecos,
        "equal_cells": equal_n,
        "global_permutations": n_perm,
        "next_step": "Review richness and representation concordance, then run Step 10D neighboring-cell cross-boundary turnover for total, C3, and N0.",
    }
    (out_dir / "10C_analysis_summary.json").write_text(json.dumps(summary_json, indent=2, ensure_ascii=False), encoding="utf-8")
    write_readme(out_dir, summary_json)

    print("STEP 10C COMPLETE")
    print(f"AUDIT_STATUS={summary_json['audit_status']}")
    print(f"OUTPUT_DIR={out_dir}")
    print(f"INCIDENCE_SOURCE={incidence_source}")
    print(f"TRAIT_SOURCE={trait_meta['source']}")
    print(f"FORMAL_ECOREGIONS={len(formal_ecos)}")
    print(f"EQUAL_CELLS={equal_n}")
    print(f"ITERATIONS={args.iterations}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception as exc:
        print("STEP 10C FAILED", file=sys.stderr)
        print(str(exc), file=sys.stderr)
        traceback.print_exc()
        raise SystemExit(1)
