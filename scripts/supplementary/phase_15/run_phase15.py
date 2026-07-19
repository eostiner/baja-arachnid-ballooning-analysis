#!/usr/bin/env python3
"""Run Phase 15A–15F in the frozen order."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

STEPS = ["15A", "15B", "15C", "15D", "15E", "15F"]
FILES = {
    "15A": "15A_prepare_bayesian_input.py",
    "15B": "15B_fit_bayesian_h3_model.py",
    "15C": "15C_prior_sensitivity_and_checks.py",
    "15D": "15D_leave_one_cell_out_prediction.py",
    "15E": "15E_design_power_simulation.py",
    "15F": "15F_synthesize_phase15.py",
}


def parse_args():
    p = argparse.ArgumentParser(description="Run Phase 15 Bayesian H3 analysis.")
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--from-step", choices=STEPS, default="15A")
    p.add_argument("--to-step", choices=STEPS, default="15F")
    p.add_argument("--chains", type=int, default=4)
    p.add_argument("--warmup", type=int, default=3000)
    p.add_argument("--draws", type=int, default=5000)
    p.add_argument("--power-simulations", type=int, default=2000)
    p.add_argument("--seed", type=int, default=20260718)
    p.add_argument("--quick", action="store_true", help="Development check only; not manuscript inference.")
    return p.parse_args()


def main():
    args = parse_args()
    here = Path(__file__).resolve().parent
    start, stop = STEPS.index(args.from_step), STEPS.index(args.to_step)
    if start > stop:
        raise ValueError("--from-step must not follow --to-step")
    chains = 2 if args.quick else args.chains
    warmup = 1000 if args.quick else args.warmup
    draws = 1000 if args.quick else args.draws
    power = 250 if args.quick else args.power_simulations
    for step in STEPS[start:stop + 1]:
        cmd = [sys.executable, str(here / FILES[step]), "--project-root", str(args.project_root), "--seed", str(args.seed)] if step in {"15B", "15C", "15D", "15E"} else [sys.executable, str(here / FILES[step]), "--project-root", str(args.project_root)]
        if step == "15B":
            cmd += ["--chains", str(chains), "--warmup", str(warmup), "--draws", str(draws)]
        elif step == "15C":
            cmd += ["--chains", str(chains), "--warmup", str(max(1000, warmup // 2)), "--draws", str(max(1000, draws // 2))]
        elif step == "15D":
            cmd += ["--chains", "2", "--warmup", str(600 if args.quick else 1200), "--draws", str(800 if args.quick else 1800)]
        elif step == "15E":
            cmd += ["--simulations", str(power)]
        print("\n" + "=" * 78)
        print(f"STEP {step}: {' '.join(cmd)}")
        print("=" * 78, flush=True)
        subprocess.run(cmd, check=True)
    print("\nPHASE 15 COMPLETED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
