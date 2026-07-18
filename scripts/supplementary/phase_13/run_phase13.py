#!/usr/bin/env python3
"""Run the retained Phase 13 extension in canonical order."""
import argparse, subprocess, sys
from pathlib import Path
ORDER=("13A","13B","13C","13D","13E","13F","13G","13H1","13H2","13H3")
FILES={
"13A":"13A_audit_inputs.py","13B":"13B_build_historical_boundary_signals.py",
"13C":"13C_build_contemporary_environment_signal.py","13D":"13D_construct_paired_C3_N0_community_dissimilarities.py",
"13E":"13E_join_pairwise_historical_environment_geography_community.py","13F":"13F_test_historical_vs_contemporary_C3_N0.py",
"13G":"13G_robustness_sensitivity_phase13.py","13H1":"13H1_audit_mulege_ecotone_inputs.py",
"13H2":"13H2_build_mulege_ecotone_predictors.py","13H3":"13H3_test_mulege_ecotone_anomaly.py"}
p=argparse.ArgumentParser(); p.add_argument("--project-root",required=True,type=Path); p.add_argument("--from-step",choices=ORDER,default="13A"); p.add_argument("--to-step",choices=ORDER,default="13H3"); p.add_argument("--dry-run",action="store_true"); a=p.parse_args()
here=Path(__file__).resolve().parent; cfg=here/"configs"; start=ORDER.index(a.from_step); stop=ORDER.index(a.to_step)
if stop<start: raise SystemExit("--to-step precedes --from-step")
for step in ORDER[start:stop+1]:
    cmd=[sys.executable,str(here/FILES[step]),"--project-root",str(a.project_root)]
    if step=="13B": cmd += ["--config",str(cfg/"phase_13_boundaries_frozen.csv")]
    if step=="13C": cmd += ["--config",str(cfg/"phase_13_contemporary_predictors_frozen.csv")]
    if step=="13F": cmd += ["--permutations","4999"]
    if step=="13G": cmd += ["--permutations","1999"]
    if step=="13H3": cmd += ["--bootstrap","1999"]
    print("\n"+"="*88+f"\nRUNNING {step}\n"+" ".join(cmd))
    if not a.dry_run: subprocess.run(cmd,check=True)
print("\nPhase 13 run complete.")
