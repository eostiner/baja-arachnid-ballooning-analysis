#!/usr/bin/env python3
"""Check publication inputs and validate the authoritative C3/N0/D4 trait schema."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from normalize_trait_table import normalize_project_traits


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project_root", type=Path)
    args = parser.parse_args()
    root = args.project_root.expanduser().resolve()

    required = [
        root / "02_data_clean/05_final_qc_flags/05_biodiversity_final_records.tsv",
        root / "02_data_clean/05_final_qc_flags/05_ballooning_final_records.tsv",
    ]
    missing = [path for path in required if not path.is_file()]
    for path in required:
        print(("FOUND   " if path.is_file() else "MISSING ") + str(path))
    if missing:
        return 1

    try:
        result = normalize_project_traits(root, write=False)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"FOUND TRAIT SOURCE: {result.source}")
    print(f"TRAIT TIER FIELD: {result.tier_field or '<derived from primary group>'}")
    print(f"TRAIT PRIMARY GROUP FIELD: {result.group_field or '<derived from tier>'}")
    print(
        "PRIMARY CLASSES: "
        f"C3={result.counts['C3']}, N0={result.counts['N0']}, "
        f"D4_excluded={result.counts['D4_excluded']}"
    )
    print("EVIDENCE COUNTS: " + ", ".join(f"{k}={v}" for k, v in sorted(result.evidence_counts.items())))
    print("PREFLIGHT PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
