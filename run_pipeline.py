#!/usr/bin/env python3
"""Run the retained Baja Ballooning publication workflow in dependency order."""
from __future__ import annotations
import argparse, shutil, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEPS = [
    ("10", "python", "scripts/core/step_10/10_build_25km_genus_grid_incidence.py", ["--project-root", "{root}"]),
    ("10T", "python", "tools/normalize_trait_table.py", ["{root}"]),
    ("11", "R", "scripts/core/step_11/11_latitude_band_diversity_turnover.R", ["{root}"]),
    ("11B", "python", "scripts/core/step_11B/11B_equal_cell_resampling.py", ["{root}", "5000", "20260713"]),
    ("11F", "python", "scripts/core/step_11F/11F_equal_cell_rarefaction.py", ["{root}", "5000", "20260712"]),
    ("11G", "python", "scripts/core/step_11G/11G_equal_cell_trait_partitioned_richness_turnover.py", ["{root}", "5000", "20260713"]),
    ("11E", "python", "scripts/mapping/step_11E/11E_clean_jaccard_map.py", ["{root}"]),
    ("10A", "R", "scripts/biogeography/step_10A/10A_download_and_audit_baja_ecoregions.R", ["{root}", "20260715"]),
    ("10B", "R", "scripts/biogeography/step_10B/10B_build_cell_ecoregion_crosswalk.R", ["{root}", "20260715"]),
    ("10C", "python", "scripts/biogeography/step_10C/10C_equal_cell_ecoregion_richness.py", ["{root}", "20260715", "5000"]),
    ("10D", "python", "scripts/biogeography/step_10D/10D_ecoregion_boundary_turnover.py", ["{root}", "20260715", "5000"]),
    ("10E", "python", "scripts/biogeography/step_10E/10E_published_break_concordance.py", ["{root}", "20260715", "5000"]),
    ("10G", "python", "scripts/biogeography/step_10G/10G_build_map_and_band_ratio.py", ["{root}", "20260715", "5000"]),
    ("10I", "python", "scripts/mapping/step_10I/10I_build_final_map_and_ratios_v7.py", ["{root}", "--dpi", "400"]),
    ("12A", "R", "scripts/environment/step_12A/12A_environmental_glm_input_audit.R", ["{root}"]),
    ("12B", "R", "scripts/environment/step_12B/12B_install_and_audit_GEE_exports.R", ["{root}"]),
    ("12C", "R", "scripts/environment/step_12C/12C_build_cell_environment_model_table.R", ["{root}"]),
    ("12D", "R", "scripts/environment/step_12D/12D_screen_predictors_compare_models.R", ["{root}"]),
    ("12E", "R", "scripts/environment/step_12E/12E_test_coefficient_stability_and_reduced_models.R", ["{root}"]),
    ("12F", "R", "scripts/environment/step_12F/12F_prepare_environment_predictor_scores.R", ["{root}"]),
    ("12K", "R", "scripts/environment/step_12K/12K_spatial_plus_trait_composition.R", ["{root}", "paper", "20260714"]),
    ("12K1", "R", "scripts/mapping/step_12K1/12K1_make_combined_spatialplus_map.R", ["{root}"]),
]

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True, type=Path)
    ap.add_argument("--from-step", default=STEPS[0][0])
    ap.add_argument("--to-step", default=STEPS[-1][0])
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--skip-step", action="append", default=[], help="Step ID to skip; repeat as needed.")
    ap.add_argument("--environment-source", type=Path, help="Optional Step 12B GEE export folder; defaults to ~/Downloads/Baja_Ballooning_12B.")
    ap.add_argument("--continue-on-error", action="store_true")
    ap.add_argument("--include-gdm", action="store_true", help="Run supplementary Step 12J v6 after the retained pipeline.")
    ap.add_argument("--include-omi", action="store_true", help="Run exploratory Step 12N OMI niche analysis after the retained pipeline.")
    args = ap.parse_args()
    root = args.project_root.expanduser().resolve()
    if not root.is_dir(): ap.error(f"Project root not found: {root}")
    ids = [x[0] for x in STEPS]
    if args.from_step not in ids or args.to_step not in ids:
        ap.error(f"Steps must be one of: {', '.join(ids)}")
    start_index = ids.index(args.from_step)
    end_index = ids.index(args.to_step)
    if start_index > end_index:
        ap.error("--from-step must occur before or equal to --to-step in PIPELINE_ORDER.md")
    unknown_skips = sorted(set(args.skip_step) - set(ids))
    if unknown_skips:
        ap.error(f"Unknown --skip-step value(s): {', '.join(unknown_skips)}")
    selected = [row for row in STEPS[start_index:end_index + 1] if row[0] not in set(args.skip_step)]
    needs_r = any(kind == "R" for _, kind, _, _ in selected) or args.include_gdm or args.include_omi
    if not args.dry_run and needs_r and shutil.which("Rscript") is None:
        print("ERROR: Rscript was not found in PATH.", file=sys.stderr)
        return 2
    failures = []
    for step, kind, rel, template_args in selected:
        exe = sys.executable if kind == "python" else "Rscript"
        rendered_args = [x.format(root=str(root)) for x in template_args]
        if step == "12B" and args.environment_source is not None:
            rendered_args.append(str(args.environment_source.expanduser().resolve()))
        cmd = [exe, str(HERE/rel)] + rendered_args
        print("\n" + "="*78); print(f"STEP {step}: {' '.join(cmd)}"); print("="*78)
        if args.dry_run: continue
        result = subprocess.run(cmd)
        if result.returncode:
            failures.append((step, result.returncode))
            if not args.continue_on_error: break
    if args.include_gdm and not failures:
        cmd = ["Rscript", str(HERE/"scripts/supplementary/step_12J/12J_master_trait_stratified_gdm.R"), str(root), "paper", "199", "500", "20260714", "4"]
        print("\nSUPPLEMENTARY 12J:", " ".join(cmd)); result = subprocess.run(cmd)
        if result.returncode: failures.append(("12J", result.returncode))
    if args.include_omi and not failures:
        cmd = [
            "Rscript",
            str(HERE/"scripts/environment/step_12N/12N_omi_environmental_niche_analysis.R"),
            str(root),
            "5000",
            "20260717",
        ]
        print("\nEXPLORATORY 12N:", " ".join(cmd))
        if not args.dry_run:
            result = subprocess.run(cmd)
            if result.returncode:
                failures.append(("12N", result.returncode))

    if not failures:
        if not args.dry_run:
            subprocess.run([sys.executable, str(HERE/"tools/collect_publication_figures.py"), str(root)], check=False)
            snapshot_dir = root / "04_analysis" / "REPRODUCIBILITY_SNAPSHOT"
            snapshot_dir.mkdir(parents=True, exist_ok=True)
            with (snapshot_dir / "python_pip_freeze.txt").open("w", encoding="utf-8") as stream:
                subprocess.run(
                    [sys.executable, "-m", "pip", "freeze"],
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            if shutil.which("Rscript") is not None:
                r_output = str(snapshot_dir / "R_sessionInfo.txt").replace("\\", "/")
                r_expr = f'writeLines(capture.output(sessionInfo()), "{r_output}")'
                subprocess.run(["Rscript", "-e", r_expr], check=False)
        print("\nPIPELINE DRY RUN COMPLETE" if args.dry_run else "\nPIPELINE COMPLETED SUCCESSFULLY"); return 0
    print("\nFAILED STEPS:", failures, file=sys.stderr); return 1
if __name__ == "__main__": raise SystemExit(main())
