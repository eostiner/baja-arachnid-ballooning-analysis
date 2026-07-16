#!/usr/bin/env python3
"""Collect available manuscript figures into one publication-output folder."""
from __future__ import annotations
import argparse
import shutil
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project_root", type=Path)
    args = parser.parse_args()
    root = args.project_root.expanduser().resolve()
    if not root.is_dir():
        parser.error(f"Project root not found: {root}")
    out = root / "04_analysis" / "PUBLICATION_FIGURES"
    out.mkdir(parents=True, exist_ok=True)
    patterns = {
        "Figure_1_Equal_Cell_Rarefaction": root / "04_analysis/11F_equal_cell_rarefaction/figures/Figure_1_equal_cell_rarefaction_primary",
        "Figure_2A_Equal_Cell_Jaccard_Map": root / "04_analysis/11E_clean_jaccard_map/figures/Figure_2_Jaccard_map_clean_primary",
        "Figure_2B_Five_Band_Simpson_Replacement": root / "04_analysis/11B_equal_cell_resampling/figures/11B_equal_cell_adjacent_simpson_turnover",
        "Figure_3_Biogeographic_Dispersal_Balance": root / "04_analysis/C3_pipeline_rebuild/09_C3_biogeographic_concordance/10I_final_map_and_ratios/publication_outputs/Figure_3_Biogeographic_Dispersal_Balance",
        "Figure_4_SpatialPlus_Environmental_Effects": root / "04_analysis/12K_spatial_plus_trait_composition/figures/Figure_3_spatial_plus_environmental_effects",
        "Figure_S_SpatialPlus_Map": root / "04_analysis/12K1_combined_spatialplus_map/figures/Figure_4_observed_ballooning_and_SpatialPlus_contribution",
    }
    for label, stem in patterns.items():
        found = False
        for ext in (".png", ".pdf", ".svg", ".tif", ".jpg"):
            source = stem.with_suffix(ext)
            if source.exists():
                shutil.copy2(source, out / (label + ext))
                print("COPIED", source)
                found = True
        if not found:
            print("NOT FOUND", stem)
    print("Publication figures:", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
