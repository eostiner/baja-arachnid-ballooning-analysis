#!/usr/bin/env python3
"""Verify SHA-256 checksums recorded in SCRIPT_MANIFEST.tsv."""
from __future__ import annotations
import argparse, csv, hashlib
from pathlib import Path

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        help="Optional manifest path; defaults to SCRIPT_MANIFEST.tsv in the repository root.",
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    manifest = args.manifest.expanduser().resolve() if args.manifest else root / "SCRIPT_MANIFEST.tsv"
    if not manifest.is_file():
        parser.error(f"Manifest not found: {manifest}")
    bad: list[tuple[str, str, str]] = []
    with manifest.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            path = root / row["path"]
            actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else "MISSING"
            if actual != row["sha256"]:
                bad.append((row["path"], actual, row["sha256"]))
    if bad:
        print("MANIFEST FAIL")
        for item in bad:
            print(*item, sep="\t")
        return 1
    print("MANIFEST PASS")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
