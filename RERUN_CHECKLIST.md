# Publication rerun checklist

Use this checklist for the final manuscript rerun and again before creating a tagged GitHub release.

## 1. Create the software environment

```bash
cd baja-arachnid-ballooning-analysis
conda env create -f environment.yml
conda activate baja-ballooning
```

## 2. Verify the downloaded release

```bash
python tools/verify_manifest.py
```

## 3. Check the project inputs

```bash
python tools/preflight_check.py ~/Desktop/Baja_Ballooning_Pipeline
```

The preflight must identify the authoritative source, tier field, and primary-group field and report nonzero C3 and fixed N0 counts. It accepts the authoritative CSV or Excel workbook. Blank `exclusive_tier` values are valid for N0 rows when `primary_C3_group` explicitly marks them as non-ballooning. A legacy binary-only table is intentionally rejected.

## 4. Build and inspect the canonical trait table

The full runner performs this as Step 10T. It can also be run independently:

```bash
python tools/normalize_trait_table.py ~/Desktop/Baja_Ballooning_Pipeline
```

Confirm that it writes `07_reviewed_genus_trait_lookup_normalized.csv` and reports nonzero C3 and N0 counts.

## 5. Inspect the planned commands

```bash
python run_pipeline.py \
  --project-root ~/Desktop/Baja_Ballooning_Pipeline \
  --dry-run
```

## 6. Run the retained publication pipeline

```bash
python run_pipeline.py \
  --project-root ~/Desktop/Baja_Ballooning_Pipeline
```

If the Step 12B exports are stored outside the default Downloads folder, add:

```bash
  --environment-source /path/to/Baja_Ballooning_12B
```

If the original export directory no longer exists, Step 12B automatically audits the already installed rasters in `ANALYSIS_READY_INPUTS/05_environmental_rasters_12B/`.

The exploratory GDM is not required. Run it only for the supplement with `--include-gdm`.

## 7. Check the primary outputs

The runner collects the available manuscript figures in:

```text
PROJECT_ROOT/04_analysis/PUBLICATION_FIGURES/
```

The centerpiece vector figure should be present at:

```text
PROJECT_ROOT/04_analysis/C3_pipeline_rebuild/09_C3_biogeographic_concordance/
10I_final_map_and_ratios/publication_outputs/
Figure_3_Biogeographic_Dispersal_Balance.pdf
Figure_3_Biogeographic_Dispersal_Balance.svg
```

Use the PDF or SVG for manuscript production. See `docs/MAP_QA_CHECKLIST.md` before locking the figure.

## 8. Archive the run

Retain the generated input manifests, analysis summaries, and:

```text
PROJECT_ROOT/04_analysis/REPRODUCIBILITY_SNAPSHOT/
├── python_pip_freeze.txt
└── R_sessionInfo.txt
```

Record the Git commit and tag used for the rerun in the manuscript repository or analysis notes.
