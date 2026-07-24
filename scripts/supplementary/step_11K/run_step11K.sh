#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
PROJECT="${1:-}"
ITERATIONS="${2:-5000}"
SEED="${3:-20260723}"
NBOOT="${4:-200}"

if [[ -z "$PROJECT" ]]; then
  for candidate in \
    "$HOME/Desktop/BALLOONING Overflo/Baja_Ballooning_Pipeline" \
    "$HOME/Desktop/OLD BALLOONING/Baja_Ballooning_Pipeline" \
    "$HOME/Desktop/Baja_Ballooning_Pipeline"; do
    if [[ -d "$candidate" ]]; then PROJECT="$candidate"; break; fi
  done
fi
if [[ -z "$PROJECT" || ! -d "$PROJECT" ]]; then
  echo "ERROR: project root not found."
  echo "Usage: bash scripts/supplementary/step_11K/run_step11K.sh /path/to/Baja_Ballooning_Pipeline [iterations] [seed] [iNEXT_nboot]"
  exit 1
fi

PARENT="$(dirname "$PROJECT")"
AUDIT="$PARENT/D1_D4_INPUT_FILE_AUDIT 2"
OUTPUT="$PROJECT/04_analysis/11K_publication_nestedness_replacement_QC"

INCIDENCE=""
for file in \
  "$AUDIT/10_biodiversity_final_genus_by_grid25km_incidence.csv" \
  "$PROJECT/ANALYSIS_READY_INPUTS/08_grid25km_incidence/10_biodiversity_final_genus_by_grid25km_incidence.csv" \
  "$PROJECT/02_data_clean/08_grid25km_incidence/10_biodiversity_final_genus_by_grid25km_incidence.csv"; do
  [[ -f "$file" ]] && INCIDENCE="$file" && break
done

LOOKUP=""
for file in \
  "$PROJECT/ANALYSIS_READY_INPUTS/04_spatial_reference/10_common_grid25km_cell_lookup.csv" \
  "$PROJECT/02_data_clean/08_grid25km_incidence/10_common_grid25km_cell_lookup.csv"; do
  [[ -f "$file" ]] && LOOKUP="$file" && break
done

TRAITS=""
for file in \
  "$PROJECT/ANALYSIS_READY_INPUTS/03_trait_tables/07_reviewed_genus_trait_lookup_final_REAUDITED_COMPLETE.csv" \
  "$PROJECT/ANALYSIS_READY_INPUTS/03_trait_tables/07_reviewed_genus_trait_lookup_normalized.csv" \
  "$AUDIT/07_reviewed_genus_trait_lookup_final_REAUDITED_COMPLETE.csv" \
  "$AUDIT/USE_GOOD_BalloonID_Baja_Arachnid_GenusSpecies_Long_D1_D4_AUTHORITATIVE.xlsx"; do
  [[ -f "$file" ]] && TRAITS="$file" && break
done

[[ -n "$INCIDENCE" ]] || { echo "ERROR: incidence matrix not found"; exit 1; }
[[ -n "$LOOKUP" ]] || { echo "ERROR: cell lookup not found"; exit 1; }
[[ -n "$TRAITS" ]] || { echo "ERROR: complete D1-D4/N0 trait table not found"; exit 1; }
if [[ "$TRAITS" == *"/07_reviewed_genus_trait_lookup_final.csv" ]]; then
  echo "ERROR: refusing older incomplete trait lookup: $TRAITS"
  exit 1
fi
command -v python3 >/dev/null || { echo "ERROR: python3 not found"; exit 1; }
command -v Rscript >/dev/null || { echo "ERROR: Rscript not found"; exit 1; }

mkdir -p "$OUTPUT"
python3 "$REPO_ROOT/tests/test_step11K_baselga_formulas.py"
python3 -u "$SCRIPT_DIR/11K_publication_qc_core.py" \
  "$PROJECT" \
  --iterations "$ITERATIONS" \
  --seed "$SEED" \
  --expected-equal-cells 22 \
  --incidence "$INCIDENCE" \
  --cell-lookup "$LOOKUP" \
  --traits "$TRAITS" \
  --output-dir "$OUTPUT"

Rscript "$SCRIPT_DIR/11K_make_publication_qc_figures.R" "$OUTPUT"
Rscript "$SCRIPT_DIR/11K_run_iNEXT_hill.R" "$OUTPUT" "$NBOOT" "$SEED"
python3 "$SCRIPT_DIR/11K_build_publication_text.py" "$OUTPUT"

echo "COMPLETE" > "$OUTPUT/COMPLETE_WORKFLOW_STATUS.txt"
echo "core_beta_diversity=present" >> "$OUTPUT/COMPLETE_WORKFLOW_STATUS.txt"
echo "iNEXT_hill_diversity=completed" >> "$OUTPUT/COMPLETE_WORKFLOW_STATUS.txt"
echo "nboot=$NBOOT" >> "$OUTPUT/COMPLETE_WORKFLOW_STATUS.txt"
echo "seed=$SEED" >> "$OUTPUT/COMPLETE_WORKFLOW_STATUS.txt"

echo "STEP 11K COMPLETE: $OUTPUT"
