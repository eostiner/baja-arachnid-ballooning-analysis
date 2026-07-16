#!/usr/bin/env python3
"""
Step 10E — a priori published-break concordance for the Baja C3 pipeline.

Frozen biological definition:
  C3 = D1 + D2 + D3
  N0 = fixed non-ballooning reference

A priori latitude anchors (defined from literature before testing):
  24 N — Isthmus of La Paz
  26 N — Loreto
  28 N — mid-peninsular / Vizcaino transition
  30 N — Mediterranean-desert transition

Primary local test:
  - one-degree flank south and north of each anchor
  - 16 occupied 25-km cells drawn without replacement from each flank
  - 5,000 iterations
  - total, C3, and N0 assemblages analyzed separately
  - pooled richness discontinuity, Jaccard dissimilarity, and Simpson replacement
  - permutation null randomly splits 32 cells from the same two-degree local pool

Secondary consistency analysis:
  - the five existing 2-degree latitude bands
  - 22 cells sampled from each band per iteration
  - adjacent-band turnover at 24, 26, 28, and 30 N
  - boundary rank frequencies and C3-N0 paired contrasts

The script refuses to proceed unless it reproduces the frozen 205 occupied cells,
267 genera, and exact per-cell C3 and N0 counts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import traceback
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

STEP = "10E"
EXPECTED_CELLS = 205
EXPECTED_GENERA = 267
EXPECTED_LOCAL_EQUAL_CELLS = 16
EXPECTED_BAND_EQUAL_CELLS = 22
DEFAULT_ITERATIONS = 5000


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
    ], "Step 10B cell table")
    incidence_candidates = [
        project_root / "02_data_clean/08_grid25km_incidence/10_ballooning_final_genus_grid25km_incidence_long.csv",
    ]
    trait_candidates = [
        project_root / "04_analysis/C3_pipeline_rebuild/01_trait_merge/C3_authoritative_trait_master.csv",
    ]
    # Reuse source paths independently validated by Steps 10C/10D where available.
    for manifest in [
        bio_root / "10D_ecoregion_boundary_turnover/10D_input_manifest.csv",
        bio_root / "10C_equal_cell_ecoregion_richness/10C_input_manifest.csv",
    ]:
        if manifest.exists():
            m = read_csv_flexible(manifest)
            if {"role", "file"}.issubset(m.columns):
                rows = m.loc[m["role"].astype(str) == "genus_cell_incidence", "file"]
                if len(rows):
                    incidence_candidates.insert(0, Path(str(rows.iloc[0])))
                rows = m.loc[m["role"].astype(str) == "genus_trait_lookup", "file"]
                if len(rows):
                    trait_candidates.insert(0, Path(str(rows.iloc[0])))
    return {
        "crosswalk": crosswalk,
        "incidence": find_existing(incidence_candidates, "genus-by-cell incidence table"),
        "trait": find_existing(trait_candidates, "C3/N0 trait table"),
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


def pooled(cell_ids: np.ndarray | list[str], sets_by_cell: dict[str, set[str]]) -> set[str]:
    ans: set[str] = set()
    for cid in cell_ids:
        ans.update(sets_by_cell.get(str(cid), set()))
    return ans


def beta_metrics(a_set: set[str], b_set: set[str]) -> tuple[float, float, int, int, int]:
    a = len(a_set & b_set)
    b = len(a_set - b_set)
    c = len(b_set - a_set)
    union = a + b + c
    j = (b + c) / union if union > 0 else np.nan
    mn = min(b, c)
    den = a + mn
    s = mn / den if den > 0 else np.nan
    return j, s, a, b, c


def q025(x: pd.Series | np.ndarray) -> float:
    z = np.asarray(x, dtype=float)
    z = z[np.isfinite(z)]
    return float(np.quantile(z, 0.025)) if len(z) else np.nan


def q975(x: pd.Series | np.ndarray) -> float:
    z = np.asarray(x, dtype=float)
    z = z[np.isfinite(z)]
    return float(np.quantile(z, 0.975)) if len(z) else np.nan


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


def p_upper(null: np.ndarray, obs: float) -> float:
    z = np.asarray(null, dtype=float)
    z = z[np.isfinite(z)]
    if not np.isfinite(obs) or not len(z):
        return np.nan
    return float((1 + np.sum(z >= obs)) / (len(z) + 1))


def p_two_sign(draws: np.ndarray) -> float:
    z = np.asarray(draws, dtype=float)
    z = z[np.isfinite(z)]
    if not len(z):
        return np.nan
    p = 2 * min(np.mean(z <= 0), np.mean(z >= 0))
    return float(min(1.0, max(1 / (len(z) + 1), p)))


def local_analysis(cw: pd.DataFrame, set_maps: dict[str, dict[str, set[str]]], registry: pd.DataFrame,
                   rng: np.random.Generator, iterations: int, equal_n: int):
    draw_rows: list[dict[str, Any]] = []
    null_rows: list[dict[str, Any]] = []
    count_rows: list[dict[str, Any]] = []
    for br in registry.itertuples(index=False):
        anchor = float(br.anchor_latitude)
        south = cw.loc[(cw["centroid_latitude"] >= anchor - 1) & (cw["centroid_latitude"] < anchor), "cell_id"].to_numpy(str)
        north = cw.loc[(cw["centroid_latitude"] >= anchor) & (cw["centroid_latitude"] < anchor + 1), "cell_id"].to_numpy(str)
        pool_ids = np.concatenate([south, north])
        count_rows.append({
            "break_id": br.break_id, "break_name": br.break_name, "anchor_latitude": anchor,
            "south_available_cells": len(south), "north_available_cells": len(north),
            "local_pool_cells": len(pool_ids), "equal_cells_per_flank": equal_n,
            "south_window": f"[{anchor-1:.1f}, {anchor:.1f})",
            "north_window": f"[{anchor:.1f}, {anchor+1:.1f})",
        })
        if len(south) < equal_n or len(north) < equal_n or len(pool_ids) < 2 * equal_n:
            raise RuntimeError(f"Break {br.break_id} has insufficient cells for frozen n={equal_n}: south={len(south)}, north={len(north)}")
        for it in range(iterations):
            sids = rng.choice(south, size=equal_n, replace=False)
            nids = rng.choice(north, size=equal_n, replace=False)
            selected = rng.choice(pool_ids, size=2 * equal_n, replace=False)
            rng.shuffle(selected)
            null_s = selected[:equal_n]
            null_n = selected[equal_n:]
            for assemblage, smap in set_maps.items():
                ss = pooled(sids, smap); nn = pooled(nids, smap)
                j, sim, a, b, c = beta_metrics(ss, nn)
                draw_rows.append({
                    "break_id": br.break_id, "break_name": br.break_name, "anchor_latitude": anchor,
                    "iteration": it + 1, "assemblage": assemblage,
                    "south_richness": len(ss), "north_richness": len(nn),
                    "north_minus_south_richness": len(nn) - len(ss),
                    "absolute_richness_difference": abs(len(nn) - len(ss)),
                    "jaccard": j, "simpson": sim,
                    "shared_a": a, "south_unique_b": b, "north_unique_c": c,
                })
                ns = pooled(null_s, smap); nnn = pooled(null_n, smap)
                nj, nsim, _, _, _ = beta_metrics(ns, nnn)
                null_rows.append({
                    "break_id": br.break_id, "iteration": it + 1, "assemblage": assemblage,
                    "absolute_richness_difference": abs(len(nnn) - len(ns)),
                    "jaccard": nj, "simpson": nsim,
                })
    return pd.DataFrame(draw_rows), pd.DataFrame(null_rows), pd.DataFrame(count_rows)


def summarize_local(draws: pd.DataFrame, nulls: pd.DataFrame):
    turnover_rows = []
    richness_rows = []
    trait_rows = []
    for (bid, bname, anchor, ass), g in draws.groupby(["break_id", "break_name", "anchor_latitude", "assemblage"], sort=False):
        ng = nulls[(nulls["break_id"] == bid) & (nulls["assemblage"] == ass)]
        for metric in ["jaccard", "simpson"]:
            obs = mean_finite(g[metric])
            null = ng[metric].to_numpy(float)
            turnover_rows.append({
                "break_id": bid, "break_name": bname, "anchor_latitude": anchor,
                "assemblage": ass, "metric": metric,
                "observed_mean": obs, "observed_q025": q025(g[metric]), "observed_q975": q975(g[metric]),
                "null_mean": mean_finite(null), "null_q025": q025(null), "null_q975": q975(null),
                "observed_minus_null_mean": obs - mean_finite(null),
                "permutation_p_one_sided_greater": p_upper(null, obs),
            })
        diff = g["north_minus_south_richness"].to_numpy(float)
        abs_obs = mean_finite(g["absolute_richness_difference"])
        abs_null = ng["absolute_richness_difference"].to_numpy(float)
        richness_rows.append({
            "break_id": bid, "break_name": bname, "anchor_latitude": anchor, "assemblage": ass,
            "south_mean_richness": mean_finite(g["south_richness"]),
            "north_mean_richness": mean_finite(g["north_richness"]),
            "north_minus_south_mean": mean_finite(diff),
            "north_minus_south_q025": q025(diff), "north_minus_south_q975": q975(diff),
            "absolute_difference_mean": abs_obs,
            "null_absolute_difference_mean": mean_finite(abs_null),
            "permutation_p_absolute_difference": p_upper(abs_null, abs_obs),
        })
    turn = pd.DataFrame(turnover_rows)
    rich = pd.DataFrame(richness_rows)
    turn["BH_q_all_local_turnover_tests"] = bh_adjust(turn["permutation_p_one_sided_greater"])
    rich["BH_q_all_local_richness_tests"] = bh_adjust(rich["permutation_p_absolute_difference"])

    for (bid, bname, anchor), g in draws.groupby(["break_id", "break_name", "anchor_latitude"], sort=False):
        c3 = g[g["assemblage"] == "C3"].sort_values("iteration")
        n0 = g[g["assemblage"] == "N0"].sort_values("iteration")
        if len(c3) != len(n0):
            continue
        for metric in ["jaccard", "simpson"]:
            d = c3[metric].to_numpy(float) - n0[metric].to_numpy(float)
            trait_rows.append({
                "break_id": bid, "break_name": bname, "anchor_latitude": anchor,
                "metric": metric, "C3_minus_N0_mean": mean_finite(d),
                "C3_minus_N0_q025": q025(d), "C3_minus_N0_q975": q975(d),
                "paired_two_sided_sign_p": p_two_sign(d),
            })
    trait = pd.DataFrame(trait_rows)
    if len(trait):
        trait["BH_q_trait_contrasts"] = bh_adjust(trait["paired_two_sided_sign_p"])
    return turn, rich, trait


def band_analysis(cw: pd.DataFrame, set_maps: dict[str, dict[str, set[str]]], rng: np.random.Generator,
                  iterations: int, equal_n: int):
    bands = [
        ("23-24N", 23.0, 24.0), ("24-26N", 24.0, 26.0), ("26-28N", 26.0, 28.0),
        ("28-30N", 28.0, 30.0), ("30-32N", 30.0, 32.0000001),
    ]
    boundaries = [24.0, 26.0, 28.0, 30.0]
    ids_by_band = {}
    count_rows = []
    for name, lo, hi in bands:
        ids = cw.loc[(cw["centroid_latitude"] >= lo) & (cw["centroid_latitude"] < hi), "cell_id"].to_numpy(str)
        ids_by_band[name] = ids
        count_rows.append({"latitude_band": name, "available_cells": len(ids), "equal_cells": equal_n})
        if len(ids) < equal_n:
            raise RuntimeError(f"Band {name} has only {len(ids)} cells; frozen equal n is {equal_n}")
    rows = []
    for it in range(iterations):
        sampled = {name: rng.choice(ids, size=equal_n, replace=False) for name, ids in ids_by_band.items()}
        for ass, smap in set_maps.items():
            pooled_band = {name: pooled(ids, smap) for name, ids in sampled.items()}
            for idx, boundary in enumerate(boundaries):
                south_name = bands[idx][0]; north_name = bands[idx + 1][0]
                j, sim, a, b, c = beta_metrics(pooled_band[south_name], pooled_band[north_name])
                rows.append({
                    "iteration": it + 1, "assemblage": ass, "boundary_latitude": boundary,
                    "south_band": south_name, "north_band": north_name,
                    "jaccard": j, "simpson": sim, "shared_a": a, "south_unique_b": b, "north_unique_c": c,
                })
    draws = pd.DataFrame(rows)
    sums = []
    for (ass, bound, sb, nb), g in draws.groupby(["assemblage", "boundary_latitude", "south_band", "north_band"], sort=False):
        for metric in ["jaccard", "simpson"]:
            sums.append({
                "assemblage": ass, "boundary_latitude": bound, "south_band": sb, "north_band": nb,
                "metric": metric, "mean": mean_finite(g[metric]), "q025": q025(g[metric]), "q975": q975(g[metric]),
            })
    summaries = pd.DataFrame(sums)
    ranks = []
    for (ass, metric), g in draws.melt(
        id_vars=["iteration", "assemblage", "boundary_latitude"], value_vars=["jaccard", "simpson"],
        var_name="metric", value_name="value"
    ).groupby(["assemblage", "metric"], sort=False):
        pivot = g.pivot(index="iteration", columns="boundary_latitude", values="value")
        # Split ties equally across tied maxima.
        maxv = pivot.max(axis=1)
        tied = pivot.eq(maxv, axis=0)
        weights = tied.div(tied.sum(axis=1), axis=0)
        for b in pivot.columns:
            ranks.append({
                "assemblage": ass, "metric": metric, "boundary_latitude": float(b),
                "frequency_highest_weighted": float(weights[b].mean()),
                "mean_rank": float(pivot.rank(axis=1, ascending=False, method="average")[b].mean()),
            })
    ranks = pd.DataFrame(ranks)
    trait_rows = []
    for bound in boundaries:
        for metric in ["jaccard", "simpson"]:
            c3 = draws[(draws.assemblage == "C3") & (draws.boundary_latitude == bound)].sort_values("iteration")[metric].to_numpy(float)
            n0 = draws[(draws.assemblage == "N0") & (draws.boundary_latitude == bound)].sort_values("iteration")[metric].to_numpy(float)
            d = c3 - n0
            trait_rows.append({
                "boundary_latitude": bound, "metric": metric,
                "C3_minus_N0_mean": mean_finite(d), "C3_minus_N0_q025": q025(d), "C3_minus_N0_q975": q975(d),
                "paired_two_sided_sign_p": p_two_sign(d),
            })
    trait = pd.DataFrame(trait_rows)
    trait["BH_q_band_trait_contrasts"] = bh_adjust(trait["paired_two_sided_sign_p"])
    return draws, summaries, ranks, trait, pd.DataFrame(count_rows)


def make_figures(out_dir: Path, cw: pd.DataFrame, registry: pd.DataFrame, turn: pd.DataFrame,
                 rich: pd.DataFrame, band_sum: pd.DataFrame):
    pub = out_dir / "publication_outputs"
    pub.mkdir(parents=True, exist_ok=True)
    order = registry["break_id"].tolist()
    labels = {r.break_id: f"{int(r.anchor_latitude)}°N\n{r.break_name}" for r in registry.itertuples(index=False)}
    assemblages = ["All", "C3", "N0"]

    # Map-like audit of occupied cell centroids and locked break anchors.
    fig, ax = plt.subplots(figsize=(6.2, 9.0))
    ax.scatter(cw["centroid_longitude"], cw["centroid_latitude"], s=12, alpha=0.65)
    for r in registry.itertuples(index=False):
        ax.axhline(float(r.anchor_latitude), linewidth=1.2)
        ax.text(cw["centroid_longitude"].min() + 0.1, float(r.anchor_latitude) + 0.05,
                f"{r.break_id}: {r.break_name}", fontsize=8, va="bottom")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude (°N)")
    ax.set_title("A priori published-break anchors and occupied 25-km cells")
    fig.tight_layout()
    fig.savefig(pub / "10E_published_break_anchor_map.png", dpi=300)
    fig.savefig(pub / "10E_published_break_anchor_map.pdf")
    plt.close(fig)

    for metric, pretty in [("jaccard", "Jaccard dissimilarity"), ("simpson", "Simpson replacement")]:
        fig, ax = plt.subplots(figsize=(10.2, 5.8))
        sub = turn[turn.metric == metric].copy()
        x0 = np.arange(len(order))
        offsets = {"All": -0.23, "C3": 0.0, "N0": 0.23}
        for ass in assemblages:
            g = sub[sub.assemblage == ass].set_index("break_id").reindex(order)
            x = x0 + offsets[ass]
            y = g["observed_mean"].to_numpy(float)
            lo = y - g["observed_q025"].to_numpy(float)
            hi = g["observed_q975"].to_numpy(float) - y
            ax.errorbar(x, y, yerr=[lo, hi], marker="o", linestyle="none", capsize=3, label=ass)
            ax.scatter(x, g["null_mean"].to_numpy(float), marker="x", s=45)
        ax.set_xticks(x0)
        ax.set_xticklabels([labels[x] for x in order])
        ax.set_ylabel(pretty)
        ax.set_title(f"Published-break concordance: {pretty}\nPoints = observed equal-cell means; x marks = local random-split null means")
        ax.legend(title="Assemblage")
        fig.tight_layout()
        fig.savefig(pub / f"10E_local_{metric}_published_break_concordance.png", dpi=300)
        fig.savefig(pub / f"10E_local_{metric}_published_break_concordance.pdf")
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.2, 5.8))
    x0 = np.arange(len(order))
    offsets = {"All": -0.23, "C3": 0.0, "N0": 0.23}
    for ass in assemblages:
        g = rich[rich.assemblage == ass].set_index("break_id").reindex(order)
        x = x0 + offsets[ass]
        y = g["north_minus_south_mean"].to_numpy(float)
        lo = y - g["north_minus_south_q025"].to_numpy(float)
        hi = g["north_minus_south_q975"].to_numpy(float) - y
        ax.errorbar(x, y, yerr=[lo, hi], marker="o", linestyle="none", capsize=3, label=ass)
    ax.axhline(0, linewidth=0.9)
    ax.set_xticks(x0)
    ax.set_xticklabels([labels[x] for x in order])
    ax.set_ylabel("North minus south pooled genus richness")
    ax.set_title("Equal-cell richness discontinuity across published breaks")
    ax.legend(title="Assemblage")
    fig.tight_layout()
    fig.savefig(pub / "10E_local_richness_discontinuity.png", dpi=300)
    fig.savefig(pub / "10E_local_richness_discontinuity.pdf")
    plt.close(fig)

    for metric, pretty in [("jaccard", "Jaccard dissimilarity"), ("simpson", "Simpson replacement")]:
        fig, ax = plt.subplots(figsize=(8.8, 5.5))
        sub = band_sum[band_sum.metric == metric]
        for ass in assemblages:
            g = sub[sub.assemblage == ass].sort_values("boundary_latitude")
            ax.plot(g["boundary_latitude"], g["mean"], marker="o", label=ass)
            ax.fill_between(g["boundary_latitude"], g["q025"], g["q975"], alpha=0.14)
        ax.set_xticks([24, 26, 28, 30])
        ax.set_xlabel("Boundary latitude (°N)")
        ax.set_ylabel(pretty)
        ax.set_title(f"Five-band equal-cell turnover at a priori published boundaries ({pretty})")
        ax.legend(title="Assemblage")
        fig.tight_layout()
        fig.savefig(pub / f"10E_five_band_{metric}_turnover.png", dpi=300)
        fig.savefig(pub / f"10E_five_band_{metric}_turnover.pdf")
        plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Test a priori Baja biogeographic latitudes using equal-cell C3/N0 turnover."
    )
    parser.add_argument("project_root", type=Path)
    parser.add_argument("seed", nargs="?", type=int, default=20260715)
    parser.add_argument("iterations", nargs="?", type=int, default=DEFAULT_ITERATIONS)
    args = parser.parse_args()
    project_root = args.project_root.expanduser().resolve()
    seed = args.seed
    iterations = args.iterations
    out_dir = project_root / "04_analysis/C3_pipeline_rebuild/09_C3_biogeographic_concordance/10E_published_break_concordance"
    out_dir.mkdir(parents=True, exist_ok=True)
    script_dir = Path(__file__).resolve().parent
    registry_path = script_dir / "10E_PUBLISHED_BREAK_REGISTRY.csv"
    try:
        paths = source_paths(project_root)
        manifest = pd.DataFrame([
            {"role": role, "file": str(path), "md5": md5_file(path), "size_bytes": path.stat().st_size}
            for role, path in paths.items()
        ] + [{"role": "published_break_registry", "file": str(registry_path), "md5": md5_file(registry_path), "size_bytes": registry_path.stat().st_size}])
        manifest.to_csv(out_dir / "10E_input_manifest.csv", index=False)

        cw0 = read_csv_flexible(paths["crosswalk"])
        id_col = identify_column(cw0, ["step10b_cell_id", "grid25km_id", "cell_id"], ["cell", "id"])
        lat_col = identify_column(cw0, ["centroid_latitude", "latitude", "lat"], ["latitude"])
        lon_col = identify_column(cw0, ["centroid_longitude", "longitude", "lon"], ["longitude"])
        c3_count_col = identify_column(cw0, ["C3_positive_genera", "C3"], ["c3", "genera"])
        n0_count_col = identify_column(cw0, ["N0_reference_genera", "N0"], ["n0", "genera"])
        cw = pd.DataFrame({
            "cell_id": cw0[id_col].map(clean_str),
            "centroid_latitude": pd.to_numeric(cw0[lat_col], errors="coerce"),
            "centroid_longitude": pd.to_numeric(cw0[lon_col], errors="coerce"),
            "C3_positive_genera": pd.to_numeric(cw0[c3_count_col], errors="coerce").fillna(0).astype(int),
            "N0_reference_genera": pd.to_numeric(cw0[n0_count_col], errors="coerce").fillna(0).astype(int),
        }).drop_duplicates("cell_id")
        if len(cw) != EXPECTED_CELLS:
            raise RuntimeError(f"Frozen occupied-cell count failed: expected {EXPECTED_CELLS}, found {len(cw)}")
        if cw[["centroid_latitude", "centroid_longitude"]].isna().any().any():
            raise RuntimeError("Missing cell centroid coordinates in Step 10B table.")

        incidence = load_incidence(paths["incidence"], set(cw.cell_id))
        traits = load_traits(paths["trait"])
        inc = incidence.merge(traits, on="genus", how="left")
        all_genera = set(inc.genus.unique())
        if len(all_genera) != EXPECTED_GENERA:
            raise RuntimeError(f"Frozen genus count failed: expected {EXPECTED_GENERA}, found {len(all_genera)}")
        set_maps: dict[str, dict[str, set[str]]] = {}
        for ass in ["All", "C3", "N0"]:
            if ass == "All":
                z = inc
            else:
                z = inc[inc.analysis_class.astype(str).str.upper() == ass.upper()]
            groups = z.groupby("cell_id")["genus"].apply(lambda s: set(s.astype(str))).to_dict()
            set_maps[ass] = {cid: groups.get(cid, set()) for cid in cw.cell_id}
        validation = cw[["cell_id", "C3_positive_genera", "N0_reference_genera"]].copy()
        validation["reconstructed_C3"] = validation.cell_id.map(lambda x: len(set_maps["C3"][x]))
        validation["reconstructed_N0"] = validation.cell_id.map(lambda x: len(set_maps["N0"][x]))
        validation["C3_exact"] = validation.C3_positive_genera == validation.reconstructed_C3
        validation["N0_exact"] = validation.N0_reference_genera == validation.reconstructed_N0
        validation.to_csv(out_dir / "10E_frozen_count_validation.csv", index=False)
        if not validation.C3_exact.all() or not validation.N0_exact.all():
            raise RuntimeError("Frozen per-cell C3/N0 validation failed; refusing to run Step 10E.")

        registry = read_csv_flexible(registry_path)
        expected_anchors = [24.0, 26.0, 28.0, 30.0]
        anchors = pd.to_numeric(registry.anchor_latitude, errors="coerce").tolist()
        if anchors != expected_anchors:
            raise RuntimeError(f"Published-break registry differs from design lock: {anchors}")
        registry.to_csv(out_dir / "10E_break_registry_used.csv", index=False)

        # Derive and freeze local sample size; current data should give 16.
        flank_counts = []
        for a in expected_anchors:
            flank_counts.extend([
                int(((cw.centroid_latitude >= a - 1) & (cw.centroid_latitude < a)).sum()),
                int(((cw.centroid_latitude >= a) & (cw.centroid_latitude < a + 1)).sum()),
            ])
        local_n = min(flank_counts)
        if local_n != EXPECTED_LOCAL_EQUAL_CELLS:
            raise RuntimeError(f"Frozen local equal-cell count failed: expected {EXPECTED_LOCAL_EQUAL_CELLS}, found {local_n}")

        rng_local = np.random.default_rng(seed)
        local_draws, local_nulls, local_counts = local_analysis(cw, set_maps, registry, rng_local, iterations, local_n)
        local_turn, local_rich, local_trait = summarize_local(local_draws, local_nulls)
        local_counts.to_csv(out_dir / "10E_local_window_cell_counts.csv", index=False)
        local_draws.to_csv(out_dir / "10E_local_equal_cell_draws.csv", index=False)
        local_nulls.to_csv(out_dir / "10E_local_random_split_null_draws.csv", index=False)
        local_turn.to_csv(out_dir / "10E_local_turnover_tests.csv", index=False)
        local_rich.to_csv(out_dir / "10E_local_richness_tests.csv", index=False)
        local_trait.to_csv(out_dir / "10E_local_C3_vs_N0_turnover_contrasts.csv", index=False)

        rng_band = np.random.default_rng(seed + 1000003)
        band_draws, band_sum, band_ranks, band_trait, band_counts = band_analysis(
            cw, set_maps, rng_band, iterations, EXPECTED_BAND_EQUAL_CELLS
        )
        band_counts.to_csv(out_dir / "10E_five_band_cell_counts.csv", index=False)
        band_draws.to_csv(out_dir / "10E_five_band_turnover_draws.csv", index=False)
        band_sum.to_csv(out_dir / "10E_five_band_turnover_summaries.csv", index=False)
        band_ranks.to_csv(out_dir / "10E_boundary_rank_frequencies.csv", index=False)
        band_trait.to_csv(out_dir / "10E_five_band_C3_vs_N0_contrasts.csv", index=False)

        make_figures(out_dir, cw, registry, local_turn, local_rich, band_sum)

        summary = {
            "step": STEP,
            "audit_status": "PASS_FROZEN_COUNTS_REPRODUCED_AND_PUBLISHED_BREAK_ANALYSIS_COMPLETE",
            "seed": seed, "iterations": iterations,
            "project_root": str(project_root), "output_dir": str(out_dir),
            "occupied_cells": len(cw), "total_genera": len(all_genera),
            "C3_exact_fraction": float(validation.C3_exact.mean()),
            "N0_exact_fraction": float(validation.N0_exact.mean()),
            "published_break_anchors": expected_anchors,
            "local_equal_cells_per_flank": local_n,
            "five_band_equal_cells": EXPECTED_BAND_EQUAL_CELLS,
            "local_turnover_tests": len(local_turn),
            "local_richness_tests": len(local_rich),
            "next_step": "Interpret Step 10E jointly with Steps 10C-10D, freeze the biogeographic synthesis, and package final scripts/results.",
        }
        (out_dir / "10E_analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        (out_dir / "README_RESULTS_FIRST.txt").write_text(
            "STEP 10E COMPLETE\n\n"
            "Start with:\n"
            "  10E_local_turnover_tests.csv\n"
            "  10E_local_richness_tests.csv\n"
            "  10E_local_C3_vs_N0_turnover_contrasts.csv\n"
            "  10E_five_band_turnover_summaries.csv\n"
            "  10E_boundary_rank_frequencies.csv\n"
            "  publication_outputs/\n\n"
            "Interpretation rule: published breaks are supported only when the observed local equal-cell statistic exceeds its local random-split null after considering BH correction and when the five-band analysis is directionally concordant. Absence of support is a valid negative result.\n",
            encoding="utf-8"
        )
        print("STEP 10E COMPLETE")
        print("AUDIT_STATUS=PASS_FROZEN_COUNTS_REPRODUCED_AND_PUBLISHED_BREAK_ANALYSIS_COMPLETE")
        print(f"OUTPUT_DIR={out_dir}")
        print(f"OCCUPIED_CELLS={len(cw)}")
        print(f"TOTAL_GENERA={len(all_genera)}")
        print(f"LOCAL_EQUAL_CELLS_PER_FLANK={local_n}")
        print(f"FIVE_BAND_EQUAL_CELLS={EXPECTED_BAND_EQUAL_CELLS}")
        print(f"ITERATIONS={iterations}")
        return 0
    except Exception as exc:
        (out_dir / "10E_FAILURE.txt").write_text(
            f"STEP 10E FAILED\n{type(exc).__name__}: {exc}\n\n{traceback.format_exc()}", encoding="utf-8"
        )
        print(f"STEP 10E FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
