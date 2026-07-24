#!/usr/bin/env python3
"""Minimal formula checks using fixed a,b,c counts; no biological results are generated."""
import importlib.util
import sys
from pathlib import Path

script = Path(__file__).resolve().parents[1] / "scripts" / "supplementary" / "step_11K" / "11K_publication_qc_core.py"
spec = importlib.util.spec_from_file_location("qc", script)
qc = importlib.util.module_from_spec(spec)
assert spec and spec.loader
sys.modules[spec.name] = qc
spec.loader.exec_module(qc)

# a=3 shared, b=2 unique A, c=1 unique B
# Build masks: shared {0,1,2}; A-only {3,4}; B-only {5}
a = sum(1 << i for i in [0,1,2,3,4])
b = sum(1 << i for i in [0,1,2,5])
r = qc.baselga_pair(a, b)
assert abs(r["jaccard_total"] - 3/6) < 1e-12
assert abs(r["jaccard_turnover"] - 2/5) < 1e-12
assert abs(r["jaccard_nestedness"] - (1/2 - 2/5)) < 1e-12
assert abs(r["sorensen_total"] - 3/9) < 1e-12
assert abs(r["simpson_replacement"] - 1/4) < 1e-12
assert abs(r["sorensen_nestedness"] - (1/3 - 1/4)) < 1e-12
assert abs(r["jaccard_total"] - (r["jaccard_turnover"] + r["jaccard_nestedness"])) < 1e-12
assert abs(r["sorensen_total"] - (r["simpson_replacement"] + r["sorensen_nestedness"])) < 1e-12
print("Baselga/Sorensen formula checks passed.")
