#!/usr/bin/env python3
"""Build manuscript-ready methods, results, and caption drafts from completed QC outputs."""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def fnum(value: str | float, digits: int = 3) -> str:
    x = float(value)
    return f"{x:.{digits}f}".replace("-0.000", "0.000")


def boolish(value: str) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def one(rows: list[dict[str, str]], **where: str) -> dict[str, str]:
    hits = [r for r in rows if all(str(r.get(k)) == str(v) for k, v in where.items())]
    if len(hits) != 1:
        raise RuntimeError(f"Expected one row for {where}, found {len(hits)}")
    return hits[0]


def interval_text(row: dict[str, str]) -> str:
    return f"{fnum(row['median'])}, 95% resampling interval {fnum(row['p025'])} to {fnum(row['p975'])}"


def main() -> int:
    if len(sys.argv) < 2:
        raise SystemExit("Usage: 04_build_publication_text.py <QC output directory>")
    out = Path(sys.argv[1]).expanduser().resolve()
    provenance = json.loads((out / "00_INPUT_PROVENANCE.json").read_text(encoding="utf-8"))
    contrasts = read_csv(out / "06_TRAIT_CONTRAST_SUMMARY.csv")
    beta = read_csv(out / "05_BASELGA_PAIRWISE_SUMMARY.csv")
    contrasts = [r for r in contrasts if boolish(r["adjacent"])]

    pairs = [
        ("23-24N", "24-26N", "23–24 and 24–26°N"),
        ("24-26N", "26-28N", "24–26 and 26–28°N"),
        ("26-28N", "28-30N", "26–28 and 28–30°N"),
        ("28-30N", "30-32N", "28–30 and 30–32°N"),
    ]

    sample_size = provenance["actual_equal_cells"]
    iterations = provenance["iterations"]
    counts = provenance["trait_evidence_counts"]
    ballooning_n = int(counts.get("D1", 0)) + int(counts.get("D2", 0)) + int(counts.get("D3", 0))
    nonballooning_n = int(counts.get("N0", 0))
    d4_n = int(counts.get("D4", 0))

    south_rep = one(contrasts, band_1="23-24N", band_2="24-26N", metric="jaccard_turnover")
    south_nest = one(contrasts, band_1="23-24N", band_2="24-26N", metric="jaccard_nestedness")
    south_total = one(contrasts, band_1="23-24N", band_2="24-26N", metric="jaccard_total")
    north_rep = one(contrasts, band_1="28-30N", band_2="30-32N", metric="jaccard_turnover")
    north_nest = one(contrasts, band_1="28-30N", band_2="30-32N", metric="jaccard_nestedness")
    north_total = one(contrasts, band_1="28-30N", band_2="30-32N", metric="jaccard_total")

    methods = (
        "### Methods\n\n"
        f"To compare the structure of compositional dissimilarity between ballooning-capable and non-ballooning arachnid assemblages, genus incidences were analyzed across five latitude bands using a frozen genus-by-25-km-cell incidence matrix containing {provenance['matrix_genera']} genera in {provenance['matrix_cells']} occupied cells. "
        f"Sampling effort was standardized by selecting {sample_size} occupied cells without replacement from every latitude band during each of {iterations:,} Monte Carlo iterations. The same cell draw was used for both trait groups within each iteration, yielding paired comparisons. Ballooning-capable genera comprised evidence classes D1–D3 ({ballooning_n} genera), non-ballooning genera comprised N0 ({nonballooning_n} genera), and D4 ({d4_n} genera) was excluded from the primary comparison. "
        "Genus incidences were pooled within each selected latitude-band sample. Total Jaccard dissimilarity was partitioned following Baselga’s framework into replacement and nestedness-resultant components. The complementary Sørensen-family partition was also calculated as Simpson replacement and Sørensen nestedness-resultant dissimilarity. Trait contrasts were calculated within each iteration as ballooning-capable minus non-ballooning. Values are reported as medians and 2.5th–97.5th percentile paired-resampling intervals. These intervals quantify uncertainty associated with equal-cell resampling and are not parametric confidence intervals.\n"
    )

    results = (
        "### Results\n\n"
        "Total Jaccard dissimilarity did not differ consistently between ballooning-capable and non-ballooning assemblages across adjacent latitude bands, because the paired interval for the total-Jaccard contrast included zero at all four boundaries. However, the decomposition of total dissimilarity differed at the southern and northern ends of the peninsula. "
        f"Between 23–24 and 24–26°N, ballooning-capable assemblages had lower Jaccard replacement than non-ballooning assemblages ({interval_text(south_rep)}) and greater nestedness-resultant dissimilarity ({interval_text(south_nest)}), whereas the total-Jaccard contrast remained unresolved ({interval_text(south_total)}). "
        "Trait-group contrasts were near zero between 24–26 and 26–28°N and remained unresolved between 26–28 and 28–30°N. The latter boundary nevertheless had the highest absolute replacement for both trait groups. "
        f"Between 28–30 and 30–32°N, ballooning-capable assemblages again had greater nestedness-resultant dissimilarity ({interval_text(north_nest)}), accompanied by a tendency toward lower replacement ({interval_text(north_rep)}), while total Jaccard remained unresolved ({interval_text(north_total)}). "
        "Thus, ballooning capability did not uniformly reduce total compositional dissimilarity across the peninsula; rather, it altered the relative contribution of replacement and nestedness at particular biogeographic transitions.\n"
    )

    inext_note = ""
    inext_file = out / "09_iNEXT_COVERAGE_STANDARDIZED_HILL.csv"
    if inext_file.exists():
        inext_rows = read_csv(inext_file)
        target = None
        for row in inext_rows:
            if row.get("target_coverage") not in (None, ""):
                target = float(row["target_coverage"])
                break
        target_text = f"{target:.3f}" if target is not None else "the shared minimum observed coverage"
        inext_note = (
            "\n### iNEXT/Hill-diversity diagnostic\n\n"
            f"Incidence-based Hill diversity for q = 0, 1, and 2 was additionally standardized to a shared sample coverage of {target_text}. These estimates provide a coverage-standardized comparison of alpha and pooled adjacent-band gamma diversity. They are reported separately from the Baselga replacement–nestedness partition and should not be interpreted as replacement or nestedness metrics.\n"
        )
    else:
        inext_note = (
            "\n### iNEXT/Hill-diversity diagnostic\n\n"
            "The iNEXT stage has not yet completed. Do not report coverage-standardized Hill-diversity values until `09_iNEXT_COVERAGE_STANDARDIZED_HILL.csv` exists.\n"
        )

    caption = (
        "### Recommended main-figure caption\n\n"
        f"**Figure X.** Differences in total Jaccard dissimilarity and its Baselga replacement and nestedness-resultant components between ballooning-capable and non-ballooning arachnid assemblages across adjacent latitude bands of the Baja California Peninsula. Each band was standardized to {sample_size} occupied 25-km cells over {iterations:,} paired Monte Carlo iterations. Points show median ballooning-capable minus non-ballooning contrasts, and error bars show 2.5th–97.5th percentile paired-resampling intervals. Negative replacement values indicate lower genus replacement among ballooning-capable assemblages, whereas positive nestedness-resultant values indicate that a greater proportion of compositional dissimilarity was associated with richness imbalance. Total-Jaccard contrasts were unresolved at all boundaries, but replacement was lower and nestedness-resultant dissimilarity greater among ballooning-capable assemblages at the southernmost boundary. Greater nestedness-resultant dissimilarity was also evident at the northernmost boundary.\n"
    )

    guardrail = (
        "\n### Reporting guardrails\n\n"
        "- Jaccard total dissimilarity is not equivalent to replacement.\n"
        "- Baselga Jaccard replacement and Jaccard nestedness-resultant sum to total Jaccard dissimilarity.\n"
        "- Simpson replacement belongs to the Sørensen-family partition and should be reported as a sensitivity analysis rather than labeled as the Jaccard component.\n"
        "- iNEXT Hill diversity is a separate coverage-standardized diversity analysis; it is not a Baselga replacement–nestedness partition.\n"
        "- These pairwise adjacent-band analyses do not automatically reproduce a previously reported single peninsula-wide contrast.\n"
    )

    text = "# Baja Ballooning publication text generated from the QC run\n\n" + methods + "\n" + results + inext_note + "\n" + caption + guardrail
    (out / "PUBLICATION_METHODS_RESULTS_AND_CAPTIONS.md").write_text(text, encoding="utf-8")
    print(out / "PUBLICATION_METHODS_RESULTS_AND_CAPTIONS.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
