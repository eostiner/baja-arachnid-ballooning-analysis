#!/usr/bin/env python3
"""
13D_construct_paired_C3_N0_community_dissimilarities.py

Purpose
-------
Construct paired ballooning (C3 = D1+D2+D3) and non-ballooning (N0)
community dissimilarities on EXACTLY the same Phase 13C eligible cells.

This step is construction/validation only. It does NOT test whether historical
boundaries or contemporary environment explain either assemblage.

Primary retained trait definition
---------------------------------
C3 ballooning-capable = D1 + D2 + D3
N0 non-ballooning     = N0
D4                    = excluded from primary inference

Expected retained genus totals across the 267-genus incidence matrix:
C3 = 87\nD4 excluded = 40\nN0 = 140\n\nTrait provenance:\n  retained 11G_normalized_trait_lookup.csv

Community dissimilarities
-------------------------
For each unordered cell pair and each trait class:
  a = shared genera
  b = genera only in cell i
  c = genera only in cell j

Jaccard total dissimilarity:
  (b + c) / (a + b + c)

Simpson replacement dissimilarity:
  min(b, c) / (a + min(b, c))

If both cells contain zero genera for a trait class, dissimilarity is NA.
If one cell is empty and the other is not, Jaccard = 1 and Simpson = 0.

Outputs
-------
13D_genus_trait_classification_retained.csv
13D_cell_trait_richness.csv
13D_pairwise_C3_N0_dissimilarity_wide.csv
13D_pairwise_C3_N0_dissimilarity_long.csv
13D_dissimilarity_summary.csv
13D_run_manifest.txt
13D_README.txt
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


EXPECTED = {"C3": 87, "D4_excluded": 40, "N0": 140}
EXPECTED_TOTAL = 267


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args():
    p = argparse.ArgumentParser(
        description="Construct paired C3/N0 community dissimilarities for Phase 13."
    )
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--incidence", type=Path, default=None)
    p.add_argument("--trait-lookup", type=Path, default=None)
    p.add_argument("--eligible-cells", type=Path, default=None)
    p.add_argument("--output-dir", type=Path, default=None)
    return p.parse_args()


def compute_pair_metrics(pres_i: np.ndarray, pres_j: np.ndarray):
    a = int(np.logical_and(pres_i, pres_j).sum())
    b = int(np.logical_and(pres_i, ~pres_j).sum())
    c = int(np.logical_and(~pres_i, pres_j).sum())
    ri = int(pres_i.sum())
    rj = int(pres_j.sum())

    denom_jac = a + b + c
    if denom_jac == 0:
        jaccard = np.nan
        simpson = np.nan
    else:
        jaccard = (b + c) / denom_jac
        m = min(b, c)
        denom_sim = a + m
        simpson = 0.0 if denom_sim == 0 else m / denom_sim

    return a, b, c, ri, rj, jaccard, simpson


def main():
    args = parse_args()
    root = args.project_root.expanduser().resolve()

    incidence_path = args.incidence or (
        root
        / "02_data_clean"
        / "08_grid25km_incidence"
        / "10_biodiversity_final_genus_by_grid25km_incidence.csv"
    )
    trait_path = args.trait_lookup or (
        root
        / "04_analysis_USE _THIS"
        / "11G_trait_partitioned_equal_cell"
        / "11G_normalized_trait_lookup.csv"
    )
    eligible_path = args.eligible_cells or (
        root
        / "04_analysis_USE _THIS"
        / "13_historical_vs_contemporary"
        / "13C_contemporary_environment_signal"
        / "13C_cell_environment_standardized.csv"
    )
    outdir = args.output_dir or (
        root
        / "04_analysis_USE _THIS"
        / "13_historical_vs_contemporary"
        / "13D_paired_C3_N0_community_dissimilarity"
    )

    incidence_path = incidence_path.expanduser().resolve()
    trait_path = trait_path.expanduser().resolve()
    eligible_path = eligible_path.expanduser().resolve()
    outdir = outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print("PHASE 13D — PAIRED C3/N0 COMMUNITY DISSIMILARITIES")
    print("=" * 78)
    print(f"PROJECT ROOT   : {root}")
    print(f"INCIDENCE      : {incidence_path}")
    print(f"TRAIT LOOKUP   : {trait_path}")
    print(f"ELIGIBLE CELLS : {eligible_path}")
    print(f"OUTPUT DIR     : {outdir}")

    for p in [incidence_path, trait_path, eligible_path]:
        if not p.exists():
            raise FileNotFoundError(f"Required input not found: {p}")

    inc = pd.read_csv(incidence_path)
    traits = pd.read_csv(trait_path)
    eligible = pd.read_csv(eligible_path)

    if "genus" not in inc.columns:
        raise ValueError("Incidence matrix must contain a 'genus' column.")
    if inc["genus"].duplicated().any():
        dups = inc.loc[inc["genus"].duplicated(), "genus"].tolist()[:20]
        raise ValueError(f"Incidence matrix contains duplicated genera: {dups}")

    if len(inc) != EXPECTED_TOTAL:
        raise ValueError(
            f"Expected {EXPECTED_TOTAL} retained genera; incidence contains {len(inc)}."
        )

    required_trait_cols = {"genus", "evidence_class", "analysis_class"}
    missing_trait_cols = required_trait_cols - set(traits.columns)
    if missing_trait_cols:
        raise ValueError(
            "Authoritative 11G trait lookup is missing required columns: "
            f"{sorted(missing_trait_cols)}"
        )

    print("\nAUTHORITATIVE TRAIT SOURCE")
    print("  genus column    : genus")
    print("  evidence column : evidence_class")
    print("  class column    : analysis_class")

    t = traits[["genus", "evidence_class", "analysis_class"]].copy()
    t["genus"] = t["genus"].astype(str).str.strip()
    t["evidence_class"] = t["evidence_class"].astype(str).str.strip()
    t["analysis_class"] = t["analysis_class"].astype(str).str.strip()

    if len(t) != EXPECTED_TOTAL:
        raise ValueError(
            f"Expected authoritative lookup to contain {EXPECTED_TOTAL} genera; "
            f"found {len(t)}."
        )
    if t["genus"].duplicated().any():
        dups = t.loc[t["genus"].duplicated(), "genus"].tolist()[:25]
        raise ValueError(f"Authoritative trait lookup has duplicate genera: {dups}")

    merged = inc[["genus"]].merge(t, on="genus", how="left", validate="one_to_one")

    missing = merged["analysis_class"].isna()
    if missing.any():
        names = merged.loc[missing, "genus"].tolist()
        raise ValueError(
            f"{len(names)} retained incidence genera are missing from the authoritative "
            f"11G trait lookup. First examples: {names[:25]}"
        )

    allowed_classes = {"C3", "D4_excluded", "N0"}
    observed_classes = set(merged["analysis_class"].dropna().astype(str))
    unexpected = sorted(observed_classes - allowed_classes)
    if unexpected:
        raise ValueError(
            f"Unexpected analysis_class values in authoritative trait lookup: {unexpected}"
        )

    counts = merged["analysis_class"].value_counts().to_dict()
    print("\nRETAINED PRIMARY TRAIT COUNTS")
    for k in ["C3", "D4_excluded", "N0"]:
        print(f"  {k:12s}: {counts.get(k, 0)}")

    for k, expected in EXPECTED.items():
        observed = int(counts.get(k, 0))
        if observed != expected:
            raise ValueError(
                f"Trait count mismatch for {k}: expected {expected}, observed {observed}. "
                "Stop rather than constructing Phase 13D from the wrong trait mapping."
            )

    cell_ids = eligible["grid_cell_id"].astype(str).tolist()
    if len(cell_ids) != 189:
        raise ValueError(
            f"Expected exactly 189 Phase 13C eligible cells; found {len(cell_ids)}."
        )
    if len(set(cell_ids)) != len(cell_ids):
        raise ValueError("Eligible-cell list contains duplicate grid_cell_id values.")

    incidence_cell_cols = {str(c): c for c in inc.columns if c != "genus"}
    missing_cells = [c for c in cell_ids if c not in incidence_cell_cols]
    if missing_cells:
        raise ValueError(
            f"{len(missing_cells)} Phase 13C eligible cells are absent from incidence matrix. "
            f"First examples: {missing_cells[:20]}"
        )

    # Reorder incidence columns to exactly match Phase 13C eligible cells.
    ordered_inc_cols = [incidence_cell_cols[c] for c in cell_ids]
    X = inc[ordered_inc_cols].copy()
    X.columns = cell_ids
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0)
    X = X.gt(0)

    classes = merged.set_index("genus").loc[inc["genus"], "analysis_class"].to_numpy()
    c3_mask = classes == "C3"
    n0_mask = classes == "N0"

    X_c3 = X.loc[c3_mask].to_numpy(dtype=bool).T   # cells x genera
    X_n0 = X.loc[n0_mask].to_numpy(dtype=bool).T

    richness = eligible[
        ["grid_cell_id", "grid_cell_order", "centroid_latitude",
         "centroid_longitude", "latitude_band"]
    ].copy()
    richness["C3_richness"] = X_c3.sum(axis=1)
    richness["N0_richness"] = X_n0.sum(axis=1)
    richness["classified_richness_C3_plus_N0"] = (
        richness["C3_richness"] + richness["N0_richness"]
    )

    print("\nCELL-LEVEL RICHNESS CHECK")
    print(
        richness[["C3_richness", "N0_richness", "classified_richness_C3_plus_N0"]]
        .describe()
        .to_string()
    )

    pair_rows = []
    long_rows = []

    for i, j in combinations(range(len(cell_ids)), 2):
        ci = cell_ids[i]
        cj = cell_ids[j]

        c3 = compute_pair_metrics(X_c3[i], X_c3[j])
        n0 = compute_pair_metrics(X_n0[i], X_n0[j])

        row = {
            "cell_i": ci,
            "cell_j": cj,
            "C3_shared_a": c3[0],
            "C3_unique_i_b": c3[1],
            "C3_unique_j_c": c3[2],
            "C3_richness_i": c3[3],
            "C3_richness_j": c3[4],
            "C3_jaccard": c3[5],
            "C3_simpson_replacement": c3[6],
            "N0_shared_a": n0[0],
            "N0_unique_i_b": n0[1],
            "N0_unique_j_c": n0[2],
            "N0_richness_i": n0[3],
            "N0_richness_j": n0[4],
            "N0_jaccard": n0[5],
            "N0_simpson_replacement": n0[6],
        }
        row["delta_jaccard_C3_minus_N0"] = (
            row["C3_jaccard"] - row["N0_jaccard"]
            if np.isfinite(row["C3_jaccard"]) and np.isfinite(row["N0_jaccard"])
            else np.nan
        )
        row["delta_simpson_C3_minus_N0"] = (
            row["C3_simpson_replacement"] - row["N0_simpson_replacement"]
            if np.isfinite(row["C3_simpson_replacement"])
            and np.isfinite(row["N0_simpson_replacement"])
            else np.nan
        )
        row["paired_valid_jaccard"] = bool(
            np.isfinite(row["C3_jaccard"]) and np.isfinite(row["N0_jaccard"])
        )
        row["paired_valid_simpson"] = bool(
            np.isfinite(row["C3_simpson_replacement"])
            and np.isfinite(row["N0_simpson_replacement"])
        )
        pair_rows.append(row)

        for trait, vals in [("C3", c3), ("N0", n0)]:
            long_rows.append({
                "cell_i": ci,
                "cell_j": cj,
                "trait_class": trait,
                "shared_a": vals[0],
                "unique_i_b": vals[1],
                "unique_j_c": vals[2],
                "richness_i": vals[3],
                "richness_j": vals[4],
                "jaccard": vals[5],
                "simpson_replacement": vals[6],
            })

    wide = pd.DataFrame(pair_rows)
    long = pd.DataFrame(long_rows)

    expected_pairs = len(cell_ids) * (len(cell_ids) - 1) // 2
    if len(wide) != expected_pairs:
        raise ValueError(
            f"Pair-count mismatch: expected {expected_pairs}, found {len(wide)}."
        )

    summary_rows = []
    for trait in ["C3", "N0"]:
        for metric in ["jaccard", "simpson_replacement"]:
            s = long.loc[long["trait_class"] == trait, metric].dropna()
            summary_rows.append({
                "trait_class": trait,
                "metric": metric,
                "n_valid_pairs": len(s),
                "mean": float(s.mean()),
                "sd": float(s.std(ddof=1)),
                "median": float(s.median()),
                "q025": float(s.quantile(0.025)),
                "q975": float(s.quantile(0.975)),
                "min": float(s.min()),
                "max": float(s.max()),
            })

    for metric in ["delta_jaccard_C3_minus_N0", "delta_simpson_C3_minus_N0"]:
        s = wide[metric].dropna()
        summary_rows.append({
            "trait_class": "C3_minus_N0",
            "metric": metric,
            "n_valid_pairs": len(s),
            "mean": float(s.mean()),
            "sd": float(s.std(ddof=1)),
            "median": float(s.median()),
            "q025": float(s.quantile(0.025)),
            "q975": float(s.quantile(0.975)),
            "min": float(s.min()),
            "max": float(s.max()),
        })

    summary = pd.DataFrame(summary_rows)

    trait_out = merged.sort_values(["analysis_class", "genus"]).reset_index(drop=True)

    trait_path_out = outdir / "13D_genus_trait_classification_retained.csv"
    richness_path = outdir / "13D_cell_trait_richness.csv"
    wide_path = outdir / "13D_pairwise_C3_N0_dissimilarity_wide.csv"
    long_path = outdir / "13D_pairwise_C3_N0_dissimilarity_long.csv"
    summary_path = outdir / "13D_dissimilarity_summary.csv"
    manifest_path = outdir / "13D_run_manifest.txt"
    readme_path = outdir / "13D_README.txt"

    trait_out.to_csv(trait_path_out, index=False)
    richness.to_csv(richness_path, index=False)
    wide.to_csv(wide_path, index=False)
    long.to_csv(long_path, index=False)
    summary.to_csv(summary_path, index=False)

    readme = """PHASE 13D — PAIRED C3/N0 COMMUNITY DISSIMILARITIES

This step constructs community-response distances only.

Primary trait definition:
- C3 ballooning-capable = D1 + D2 + D3 (87 retained genera)
- N0 non-ballooning = N0 (140 retained genera)
- D4 excluded from primary inference (40 retained genera)

Cells:
Exactly the 189 cells retained by Phase 13C are used.

Metrics:
- Jaccard total dissimilarity
- Simpson replacement dissimilarity

No historical-boundary or environmental explanatory test is performed here.
No pair is selected based on outcome size.
No boundary location or environmental predictor is altered.

The wide pairwise table contains both C3 and N0 values for the exact same
unordered cell pairs, enabling paired inference in the next step.
"""
    readme_path.write_text(readme, encoding="utf-8")

    manifest = [
        "PHASE 13D RUN MANIFEST",
        f"python={sys.version.replace(chr(10), ' ')}",
        f"project_root={root}",
        f"incidence={incidence_path}",
        f"incidence_sha256={sha256(incidence_path)}",
        f"trait_lookup={trait_path}",
        f"trait_lookup_sha256={sha256(trait_path)}",
        f"eligible_cells={eligible_path}",
        f"eligible_cells_sha256={sha256(eligible_path)}",
        f"n_retained_genera={len(inc)}",
        f"C3_genera={EXPECTED['C3']}",
        f"D4_excluded_genera={EXPECTED['D4_excluded']}",
        f"N0_genera={EXPECTED['N0']}",
        f"n_cells={len(cell_ids)}",
        f"n_unordered_cell_pairs={expected_pairs}",
        "primary_metrics=Jaccard total dissimilarity; Simpson replacement dissimilarity",
        "historical_boundary_testing_performed=NO",
        "environment_testing_performed=NO",
        "outcome_based_pair_selection=NO",
    ]
    manifest_path.write_text("\n".join(manifest) + "\n", encoding="utf-8")

    print("\nDISSIMILARITY SUMMARY")
    print(summary.to_string(index=False))

    print("\nPAIR VALIDITY")
    print(f"  total unordered pairs        : {expected_pairs}")
    print(f"  paired valid Jaccard         : {int(wide['paired_valid_jaccard'].sum())}")
    print(f"  paired valid Simpson         : {int(wide['paired_valid_simpson'].sum())}")

    print("\nFILES WRITTEN")
    for p in [
        trait_path_out, richness_path, wide_path, long_path,
        summary_path, manifest_path, readme_path
    ]:
        print(f"  {p}")

    print("\nSTATUS: PASS")
    print("13D constructed paired C3/N0 community dissimilarities on the same 189 cells.")
    print("No historical-boundary or contemporary-environment inference was performed.")
    print("Next step: 13E join historical, contemporary, geographic, and community pairwise signals.")


if __name__ == "__main__":
    main()
