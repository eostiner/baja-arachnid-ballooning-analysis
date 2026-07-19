#!/usr/bin/env python3
"""Run Phase 14 temporal H3 workflow in canonical order."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ORDER = ("14A", "14B", "14C0", "14C1", "14C2", "14D", "14E0", "14E1", "14E2", "14F", "14G", "14H", "14I")
FILES = {
    "14A": "14A_audit_temporal_coverage.py",
    "14B": "14B_build_temporal_community_turnover.py",
    "14C0": "14C0_prepare_eligible_cells.py",
    "14C1": "14C1_extract_real_stress_earth_engine.py",
    "14C2": "14C_join_cell_year_stress_anomalies.py",
    "14D": "14D_test_recent_stress_tracking.py",
    "14E0": "14E0_audit_ecostress_coverage.py",
    "14E1": "14E1_extract_extended_stress_earth_engine.py",
    "14E2": "14E2_join_extended_stress_summaries.py",
    "14F": "14F_test_extended_stress_sensitivities.py",
    "14G": "14G_test_sampling_continuity.py",
    "14H": "14H_test_trait_threshold_temporal_sensitivity.py",
    "14I": "14I_synthesize_h3_evidence.py",
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--from-step", choices=ORDER, default="14A")
    parser.add_argument("--to-step", choices=ORDER, default="14B")
    parser.add_argument("--iterations", type=int, default=1000, help="Event-resampling iterations for 14B, 14G and 14H.")
    parser.add_argument("--bootstrap-iterations", type=int, default=5000, help="Cell-cluster bootstrap iterations for 14D, 14F, 14G and 14H.")
    parser.add_argument("--permutation-iterations", type=int, default=10000, help="Wild-cluster Monte Carlo fallback; exact enumeration is used for <=15 cells.")
    parser.add_argument("--ee-project", help="Google Cloud project registered for Earth Engine.")
    parser.add_argument("--authenticate", action="store_true", help="Run interactive Earth Engine authentication for Earth Engine steps.")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    here = Path(__file__).resolve().parent
    start = ORDER.index(args.from_step)
    stop = ORDER.index(args.to_step)
    if start > stop:
        parser.error("--from-step must occur before or equal to --to-step")

    for step in ORDER[start:stop + 1]:
        cmd = [sys.executable, str(here / FILES[step]), "--project-root", str(args.project_root.expanduser().resolve())]
        if step in {"14B", "14G", "14H"}:
            cmd += ["--iterations", str(args.iterations)]
        if step in {"14D", "14F", "14G", "14H"}:
            cmd += [
                "--bootstrap-iterations", str(args.bootstrap_iterations),
                "--permutation-iterations", str(args.permutation_iterations),
            ]
        if step in {"14C1", "14E0", "14E1"}:
            if args.ee_project:
                cmd += ["--ee-project", args.ee_project]
            if args.authenticate:
                cmd.append("--authenticate")
        print("\n" + "=" * 78)
        print(f"STEP {step}: {' '.join(cmd)}")
        print("=" * 78)
        if args.dry_run:
            continue
        result = subprocess.run(cmd)
        if result.returncode:
            return result.returncode
    print(f"\nPHASE 14 {args.from_step}-{args.to_step} COMPLETED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
