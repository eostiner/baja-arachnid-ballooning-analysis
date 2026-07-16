#!/usr/bin/env Rscript

# STEP 12J v6 — Landscape-wide trait-stratified generalized dissimilarity analysis
#
# Purpose
# -------
# Replace the original genus-level ballooning GLM with a community-turnover
# analysis that directly compares three explanations for genus replacement:
#   1. geographic separation,
#   2. contemporary environmental differences, and
#   3. crossing an independently published Baja transition zone.
#
# Models are fitted separately for ballooning-capable and non-ballooning
# assemblages. Simpson turnover is the primary response because it isolates
# genus replacement from richness imbalance. Each transition is tested in a
# separate model. A matched across-boundary versus same-side analysis provides
# an easier-to-interpret corroborating test.
#
# Default run:
# Rscript 12J_master_trait_stratified_gdm.R \
#   ~/Desktop/Baja_Ballooning_Pipeline paper 199 500 20260714 4
#
# Positional arguments:
#   1 project_root       default ~/Desktop/Baja_Ballooning_Pipeline
#   2 run_mode           audit | quick | paper; default paper
#   3 permutations       gdm.varImp permutations; default 199
#   4 match_iterations   repeated greedy matching iterations; default 500
#   5 seed               default 20260714
#   6 cores              default min(4, detected cores)
#
# Notes
# -----
# * Model fitting and validation occur at the 25-km cell grain.
# * Published transitions are independent hypotheses and are tested one at a
#   time to limit collinearity with latitude and environment.
# * Oasis proximity, soil water, and surface-water access are optional one-at-a-time sensitivity controls; none enters the core model automatically.
# * The script does not infer vicariance, speciation, gene flow, or direct
#   ballooning events.
# * v6 validates every gdm return object and records non-estimable candidate
#   models as missing rather than halting the complete analysis.

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(path.expand(args[1]), mustWork = TRUE) else normalizePath(path.expand("~/Desktop/Baja_Ballooning_Pipeline"), mustWork = TRUE)
run_mode <- if (length(args) >= 2) tolower(args[2]) else "paper"
n_perm <- if (length(args) >= 3) as.integer(args[3]) else 199L
n_match_iterations <- if (length(args) >= 4) as.integer(args[4]) else 500L
seed <- if (length(args) >= 5) as.integer(args[5]) else 20260714L
detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (!is.finite(detected_cores) || detected_cores < 1) detected_cores <- 1L
cores <- if (length(args) >= 6) as.integer(args[6]) else min(4L, detected_cores)
cores <- max(1L, min(cores, detected_cores))

if (!run_mode %in% c("audit", "quick", "paper")) stop("run_mode must be audit, quick, or paper")
if (!is.finite(n_perm) || n_perm < 0) stop("permutations must be >= 0")
if (!is.finite(n_match_iterations) || n_match_iterations < 1) stop("match_iterations must be >= 1")
set.seed(seed)

SCRIPT_VERSION <- "12J_v6_2026-07-14"
GDM_MIN_RICHNESS <- 3L
MATCH_THRESHOLDS <- c(1L, 2L, 3L)
MATCH_WINDOWS_DEGREES <- c(0.75, 1.00, 1.25)
MAX_LOCAL_PAIR_KM <- 300
MIN_MATCHED_PAIRS <- 20L
MIN_CELLS_PER_SIDE <- 8L
MIN_UNIQUE_MATCHED_CELLS <- 15L
MIN_MEDIAN_RICHNESS <- 2
MAX_MATCHED_ABS_SMD <- 0.35
MIN_TEST_PAIRS <- 10L
MATCH_PROFILES <- data.frame(
  profile_name = c("standard", "relaxed"),
  geo_caliper = c(0.75, 1.00),
  env_caliper = c(0.75, 1.00),
  records_caliper = c(1.25, 1.50),
  midpoint_caliper = c(1.25, 1.50),
  stringsAsFactors = FALSE
)

required_packages <- c(
  "gdm", "vegan", "dplyr", "tidyr", "purrr", "readr", "tibble",
  "ggplot2", "patchwork", "scales", "svglite"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    "\nRun 12J_install_packages.R, then rerun Step 12J."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

utc_now <- function() format(Sys.time(), tz = "UTC", usetz = TRUE)

first_existing <- function(paths, label, required = TRUE) {
  paths <- path.expand(paths)
  found <- paths[file.exists(paths)]
  if (length(found)) return(normalizePath(found[[1]], mustWork = TRUE))
  if (required) stop("Could not find ", label, ". Tried:\n", paste(paths, collapse = "\n"))
  NULL
}

find_col <- function(fields, candidates, label, required = TRUE) {
  lookup <- setNames(fields, tolower(fields))
  for (candidate in candidates) {
    key <- tolower(candidate)
    if (key %in% names(lookup)) return(unname(lookup[[key]]))
  }
  if (required) stop("Could not identify ", label, ". Available fields: ", paste(fields, collapse = ", "))
  NULL
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

clamp01 <- function(x) pmin(1, pmax(0, x))

archive_existing <- function(output_dir, archive_root) {
  if (!dir.exists(output_dir) || length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) == 0) return(NULL)
  stamp <- format(Sys.time(), "%Y%m%dT%H%M%S")
  dir.create(archive_root, recursive = TRUE, showWarnings = FALSE)
  destination <- file.path(archive_root, paste0("12J_trait_stratified_gdm_", stamp))
  suffix <- 1L
  while (file.exists(destination)) {
    destination <- file.path(archive_root, paste0("12J_trait_stratified_gdm_", stamp, "_", suffix))
    suffix <- suffix + 1L
  }

  # Moving the complete directory is the most reliable archive operation when
  # the project and archive are on the same filesystem.
  moved <- file.rename(output_dir, destination)
  if (!isTRUE(moved)) {
    # Cross-filesystem fallback: create the target directory, copy its contents,
    # verify every copy operation, and remove the source only after success.
    dir.create(destination, recursive = TRUE, showWarnings = FALSE)
    entries <- list.files(output_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE)
    ok <- if (length(entries)) {
      file.copy(entries, destination, recursive = TRUE, copy.mode = TRUE, copy.date = TRUE)
    } else {
      logical(0)
    }
    if (length(entries) && !all(ok)) {
      unlink(destination, recursive = TRUE, force = TRUE)
      stop("Failed to archive prior Step 12J output")
    }
    unlink(output_dir, recursive = TRUE, force = TRUE)
  }
  destination
}

write_csv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, path, na = "")
}

capture_table <- function(x, prefix, out_dir) {
  if (is.null(x)) return(invisible(NULL))
  if (is.data.frame(x) || is.matrix(x)) {
    tab <- as.data.frame(x)
    if (!is.null(rownames(tab)) && any(nzchar(rownames(tab)))) tab <- tibble::rownames_to_column(tab, "row_name")
    write_csv_safe(tab, file.path(out_dir, paste0(prefix, ".csv")))
  }
  invisible(NULL)
}

output_dir <- file.path(project_root, "04_analysis", "12J_trait_stratified_gdm")
figure_dir <- file.path(output_dir, "figures")
model_dir <- file.path(output_dir, "models")
archive_root <- file.path(project_root, "08_archive")
archived <- archive_existing(output_dir, archive_root)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(output_dir, "12J_analysis_log.txt")
log_con <- file(log_path, open = "wt")
on.exit(close(log_con), add = TRUE)
log_msg <- function(...) {
  text <- paste0(...)
  cat(text, "\n")
  writeLines(text, log_con)
  flush(log_con)
}

log_msg("STEP 12J MASTER TRAIT-STRATIFIED GDM STARTED")
log_msg("Version: ", SCRIPT_VERSION)
log_msg("Started UTC: ", utc_now())
log_msg("Project root: ", project_root)
log_msg("Run mode: ", run_mode)
log_msg("Permutations: ", n_perm)
log_msg("Matching iterations: ", n_match_iterations)
log_msg("Seed: ", seed)
log_msg("Cores: ", cores)
if (!is.null(archived)) log_msg("Archived prior output: ", archived)
log_msg("gdm package version: ", as.character(utils::packageVersion("gdm")))

analysis_ready <- file.path(project_root, "ANALYSIS_READY_INPUTS")
grid_fallback <- file.path(project_root, "02_data_clean", "08_grid25km_incidence")
trait_fallback <- file.path(project_root, "02_data_clean", "07_final_trait_merge")
step12c_dir <- file.path(project_root, "04_analysis", "12C_cell_environment_model_table")
step12i_dir <- file.path(project_root, "04_analysis", "12I_integrated_oasis_satellite_model")

primary_matrix_path <- first_existing(c(
  file.path(analysis_ready, "02_incidence_matrices_25km", "10_biodiversity_final_genus_by_grid25km_incidence.csv"),
  file.path(grid_fallback, "10_biodiversity_final_genus_by_grid25km_incidence.csv")
), "primary genus-by-cell incidence matrix")
strict_matrix_path <- first_existing(c(
  file.path(analysis_ready, "02_incidence_matrices_25km", "10_biodiversity_taxonomy_strict_genus_by_grid25km_incidence.csv"),
  file.path(grid_fallback, "10_biodiversity_taxonomy_strict_genus_by_grid25km_incidence.csv")
), "taxonomy-strict incidence matrix", required = FALSE)
trait_path <- first_existing(c(
  file.path(analysis_ready, "03_trait_tables", "07_reviewed_genus_trait_lookup_normalized.csv"),
  file.path(analysis_ready, "03_trait_tables", "07_reviewed_genus_trait_lookup_final.csv"),
  file.path(trait_fallback, "07_reviewed_genus_trait_lookup_final.csv")
), "reviewed genus trait table")
cell_lookup_path <- first_existing(c(
  file.path(analysis_ready, "04_spatial_reference", "10_common_grid25km_cell_lookup.csv"),
  file.path(grid_fallback, "10_common_grid25km_cell_lookup.csv")
), "25-km cell lookup")
env_path <- first_existing(c(
  file.path(step12c_dir, "12C_primary_glm_candidate_table.csv"),
  file.path(analysis_ready, "06_environmental_cell_tables_12C", "12C_primary_glm_candidate_table.csv")
), "Step 12C primary environmental table")
oasis_metrics_path <- first_existing(c(
  file.path(step12i_dir, "12I_oasis_metrics_by_training_cell.csv")
), "Step 12I oasis metrics", required = FALSE)

input_manifest <- tibble::tibble(
  input = c("primary_incidence", "taxonomy_strict_incidence", "trait_lookup", "cell_lookup", "environment_table", "oasis_metrics_optional"),
  path = c(primary_matrix_path, strict_matrix_path %||% "", trait_path, cell_lookup_path, env_path, oasis_metrics_path %||% ""),
  exists = c(TRUE, !is.null(strict_matrix_path), TRUE, TRUE, TRUE, !is.null(oasis_metrics_path))
)
write_csv_safe(input_manifest, file.path(output_dir, "12J_input_manifest.csv"))
walk2(input_manifest$input, input_manifest$path, ~log_msg(.x, ": ", ifelse(nzchar(.y), .y, "not available")))

read_incidence <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(x) < 2) stop("Incidence matrix has fewer than two columns: ", path)
  genus_col <- names(x)[1]
  genera <- trimws(as.character(x[[genus_col]]))
  if (any(!nzchar(genera))) stop("Blank genus names in incidence matrix")
  if (anyDuplicated(tolower(genera))) stop("Duplicate genus names in incidence matrix")
  mat <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  if (any(!mat %in% c(0, 1), na.rm = TRUE)) stop("Incidence matrix contains values other than 0/1")
  mat[is.na(mat)] <- 0
  rownames(mat) <- genera
  list(matrix = mat, genera = genera, cells = colnames(mat))
}

align_incidence <- function(inc, genera, cells) {
  out <- matrix(0, nrow = length(genera), ncol = length(cells), dimnames = list(genera, cells))
  common_g <- intersect(tolower(genera), tolower(inc$genera))
  g_target <- match(common_g, tolower(genera))
  g_source <- match(common_g, tolower(inc$genera))
  common_c <- intersect(cells, inc$cells)
  out[g_target, common_c] <- inc$matrix[g_source, common_c, drop = FALSE]
  out
}

parse_evidence_class <- function(value) {
  text <- toupper(trimws(as.character(value)))
  if (is.na(text) || !nzchar(text)) return(NA_character_)
  hits <- regmatches(
    text,
    gregexpr(
      "(?<![A-Z0-9])(D1|D2|D3|D4|N0|C3)(?![A-Z0-9])",
      text,
      perl = TRUE
    )
  )[[1]]
  hits <- unique(hits[hits != ""])
  if (length(hits) == 1L && hits %in% c("D1", "D2", "D3", "D4", "N0", "C3")) return(hits)
  normalized <- tolower(gsub("[^a-z0-9]+", "", text))
  if (normalized %in% c(
    "nonballooning", "fixednonballooning", "referencenonballooning",
    "noballooning", "nonballooningreference"
  )) return("N0")
  if (normalized %in% c("c3", "primaryc3", "d1d2d3", "d1tod3")) return("C3")
  if (normalized %in% c("d4excluded", "excludedd4")) return("D4")
  NA_character_
}

choose_evidence_field <- function(table) {
  preferred <- c(
    "evidence_class", "final_evidence_class", "final_evidence_category",
    "evidence_category", "evidence_level", "trait_evidence_level",
    "d_level", "dlevel", "trait_class", "primary_class", "analysis_class",
    "ballooning_evidence_tier", "ballooning_evidence_category",
    "final_designation", "designation"
  )
  scores <- vapply(names(table), function(field) {
    parsed <- vapply(table[[field]], parse_evidence_class, character(1))
    n_ok <- sum(!is.na(parsed))
    if (n_ok == 0L) return(-Inf)
    fraction <- n_ok / max(1L, nrow(table))
    n_classes <- length(unique(parsed[!is.na(parsed)]))
    clean <- tolower(gsub("[^a-z0-9]+", "", field))
    bonus <- 0
    if (tolower(field) %in% tolower(preferred)) bonus <- bonus + 100
    if (grepl("evidence|tier|class|designation|decision", clean)) bonus <- bonus + 20
    if (n_classes < 2L || fraction < 0.25) return(-Inf)
    bonus + 100 * fraction + 5 * n_classes
  }, numeric(1))
  if (!any(is.finite(scores))) {
    stop(
      "Trait table lacks an explicit D1/D2/D3/D4/N0 or C3/N0 evidence field. ",
      "Legacy binary fields are rejected because they cannot distinguish D4 from fixed N0."
    )
  }
  names(which.max(scores))
}

primary_inc <- read_incidence(primary_matrix_path)
if ("fesa" %in% tolower(primary_inc$genera)) stop("Fesa remains in the primary incidence matrix")
strict_inc <- if (!is.null(strict_matrix_path)) read_incidence(strict_matrix_path) else NULL

trait_raw <- read.csv(trait_path, stringsAsFactors = FALSE, check.names = FALSE)
genus_field <- find_col(names(trait_raw), c("genus", "analysis_genus"), "trait genus")
evidence_field <- choose_evidence_field(trait_raw)
confidence_field <- find_col(names(trait_raw), c("final_confidence", "trait_final_confidence", "trait_confidence", "trait_ballooning_confidence"), "trait confidence", required = FALSE)
order_field <- find_col(names(trait_raw), c("order", "trait_order", "analysis_order"), "taxonomic order", required = FALSE)

trait_raw$genus_key <- tolower(trimws(as.character(trait_raw[[genus_field]])))
trait_raw$evidence_class_resolved <- vapply(trait_raw[[evidence_field]], parse_evidence_class, character(1))
trait_raw$analysis_class_resolved <- ifelse(
  trait_raw$evidence_class_resolved %in% c("D1", "D2", "D3", "C3"),
  "C3",
  ifelse(trait_raw$evidence_class_resolved == "N0", "N0",
         ifelse(trait_raw$evidence_class_resolved == "D4", "D4_excluded", NA_character_))
)
trait_raw$ballooning_binary_resolved <- ifelse(
  trait_raw$analysis_class_resolved == "C3", 1L,
  ifelse(trait_raw$analysis_class_resolved == "N0", 0L, -1L)
)
if (anyDuplicated(trait_raw$genus_key[nzchar(trait_raw$genus_key)])) stop("Duplicate genera in trait table")
trait_lookup <- trait_raw[match(tolower(primary_inc$genera), trait_raw$genus_key), , drop = FALSE]
if (any(is.na(trait_lookup$genus_key))) stop("Some incidence genera are absent from trait table")
if (any(is.na(trait_lookup$analysis_class_resolved))) stop("Some incidence genera lack explicit D1/D2/D3/D4/N0 or C3/N0 traits")
if (any(!trait_lookup$ballooning_binary_resolved %in% c(-1L, 0L, 1L))) stop("Invalid C3/N0/D4 trait classes")
trait_lookup$confidence_resolved <- if (!is.null(confidence_field)) toupper(trimws(as.character(trait_lookup[[confidence_field]]))) else "UNSPECIFIED"
trait_lookup$order_resolved <- if (!is.null(order_field)) trimws(as.character(trait_lookup[[order_field]])) else "UNRESOLVED"
trait_lookup$genus <- primary_inc$genera

primary_mat <- primary_inc$matrix
strict_mat <- if (!is.null(strict_inc)) align_incidence(strict_inc, primary_inc$genera, primary_inc$cells) else NULL

cell_lookup <- read.csv(cell_lookup_path, stringsAsFactors = FALSE, check.names = FALSE)
cell_id_field <- find_col(names(cell_lookup), c("grid_cell_id", "cell_id"), "cell ID")
cell_lookup$grid_cell_id <- as.character(cell_lookup[[cell_id_field]])

env <- read.csv(env_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!"grid_cell_id" %in% names(env)) stop("Environmental table lacks grid_cell_id")
env$grid_cell_id <- as.character(env$grid_cell_id)
if (anyDuplicated(env$grid_cell_id)) stop("Duplicate grid_cell_id in environmental table")

join_cols <- intersect(c("grid_cell_id", "centroid_latitude", "centroid_longitude", "centroid_x_m", "centroid_y_m", "easting_km", "northing_km", "latitude_band"), names(cell_lookup))
if (length(join_cols) > 1) {
  lookup_small <- cell_lookup[, join_cols, drop = FALSE]
  for (nm in setdiff(join_cols, "grid_cell_id")) {
    if (!nm %in% names(env)) env[[nm]] <- lookup_small[[nm]][match(env$grid_cell_id, lookup_small$grid_cell_id)]
  }
}

if (!"easting_km" %in% names(env) && "centroid_x_m" %in% names(env)) env$easting_km <- safe_numeric(env$centroid_x_m) / 1000
if (!"northing_km" %in% names(env) && "centroid_y_m" %in% names(env)) env$northing_km <- safe_numeric(env$centroid_y_m) / 1000
required_spatial <- c("grid_cell_id", "centroid_latitude", "centroid_longitude", "easting_km", "northing_km", "latitude_band")
missing_spatial <- setdiff(required_spatial, names(env))
if (length(missing_spatial)) stop("Environmental table lacks spatial fields: ", paste(missing_spatial, collapse = ", "))
for (nm in c("centroid_latitude", "centroid_longitude", "easting_km", "northing_km")) env[[nm]] <- safe_numeric(env[[nm]])

common_cells <- intersect(primary_inc$cells, env$grid_cell_id)
if (length(common_cells) < 50) stop("Fewer than 50 incidence cells match the environmental table")
env <- env[match(common_cells, env$grid_cell_id), , drop = FALSE]
primary_mat <- primary_mat[, common_cells, drop = FALSE]
if (!is.null(strict_mat)) strict_mat <- strict_mat[, common_cells, drop = FALSE]

if (!is.null(oasis_metrics_path)) {
  oasis_metrics <- read.csv(oasis_metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
  oasis_metrics$grid_cell_id <- as.character(oasis_metrics$grid_cell_id)
  keep_oasis <- intersect(c("grid_cell_id", "nearest_oasis_km", "oasis_proximity", "oasis_count_50km", "oasis_kernel_50km"), names(oasis_metrics))
  env <- left_join(env, oasis_metrics[, keep_oasis, drop = FALSE], by = "grid_cell_id")
  if ("oasis_proximity" %in% names(env)) {
    env$oasis_proximity <- safe_numeric(env$oasis_proximity)
    if (any(!is.finite(env$oasis_proximity))) env$oasis_proximity[!is.finite(env$oasis_proximity)] <- median(env$oasis_proximity, na.rm = TRUE)
  }
}

# ---------------------------------------------------------------------------
# Environmental-axis construction
# ---------------------------------------------------------------------------
transform_environment <- function(df) {
  out <- df
  numeric_candidates <- c(
    "vpd_mean_kpa", "wind_monthly_sd_ms", "soil_water_mean_frac", "precip_annual_mean_mm",
    "surface_water_occurrence_frac", "distance_to_modis_water_km", "evi_mean",
    "evi_monthly_amplitude", "evi_interannual_cv_pct", "lc_shrub_savanna_prop",
    "lc_grassland_prop", "lc_barren_sparse_prop", "relief_5km_m",
    "elevation_sd_5km_m", "slope_deg", "nearest_oasis_km", "oasis_proximity"
  )
  for (nm in intersect(numeric_candidates, names(out))) out[[nm]] <- safe_numeric(out[[nm]])
  if ("precip_annual_mean_mm" %in% names(out)) out$log_precip <- log1p(pmax(0, out$precip_annual_mean_mm))
  if ("surface_water_occurrence_frac" %in% names(out)) out$log_surface_water <- log1p(100 * pmax(0, out$surface_water_occurrence_frac))
  if ("distance_to_modis_water_km" %in% names(out)) out$neg_log_distance_water <- -log1p(pmax(0, out$distance_to_modis_water_km))
  if ("lc_shrub_savanna_prop" %in% names(out)) out$asin_shrub_savanna <- asin(sqrt(clamp01(out$lc_shrub_savanna_prop)))
  if ("lc_grassland_prop" %in% names(out)) out$asin_grassland <- asin(sqrt(clamp01(out$lc_grassland_prop)))
  if ("lc_barren_sparse_prop" %in% names(out)) out$asin_barren_sparse <- asin(sqrt(clamp01(out$lc_barren_sparse_prop)))
  if ("evi_interannual_cv_pct" %in% names(out)) out$log_evi_interannual_cv <- log1p(pmax(0, out$evi_interannual_cv_pct))
  if ("relief_5km_m" %in% names(out)) out$log_relief <- log1p(pmax(0, out$relief_5km_m))
  if ("elevation_sd_5km_m" %in% names(out)) out$log_elevation_sd <- log1p(pmax(0, out$elevation_sd_5km_m))
  if ("slope_deg" %in% names(out)) out$log_slope <- log1p(pmax(0, out$slope_deg))
  if ("nearest_oasis_km" %in% names(out) && !"oasis_proximity" %in% names(out)) out$oasis_proximity <- -log1p(pmax(0, out$nearest_oasis_km))
  out
}

env <- transform_environment(env)

median_impute_z <- function(df, source_name, output_name) {
  if (!source_name %in% names(df)) stop("Required environmental predictor missing: ", source_name)
  x <- safe_numeric(df[[source_name]])
  med <- median(x, na.rm = TRUE)
  if (!is.finite(med)) stop("No finite values for environmental predictor: ", source_name)
  missing_n <- sum(!is.finite(x))
  x[!is.finite(x)] <- med
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) stop("Environmental predictor has zero variance: ", source_name)
  df[[output_name]] <- as.numeric((x - mean(x)) / s)
  attr(df[[output_name]], "source_variable") <- source_name
  attr(df[[output_name]], "median_imputed") <- missing_n
  df
}

env <- median_impute_z(env, "vpd_mean_kpa", "vpd_z")
env <- median_impute_z(env, "wind_monthly_sd_ms", "wind_seasonality_z")

axis_definitions <- list(
  vegetation_axis = c("evi_mean", "evi_monthly_amplitude", "log_evi_interannual_cv", "asin_shrub_savanna", "asin_grassland", "asin_barren_sparse"),
  topography_axis = c("log_relief", "log_elevation_sd", "log_slope")
)
axis_anchor_definitions <- list(
  vegetation_axis = c(evi_mean = 1, evi_monthly_amplitude = 0.25, log_evi_interannual_cv = -0.25, asin_shrub_savanna = 1, asin_grassland = 0.5, asin_barren_sparse = -1),
  topography_axis = c(log_relief = 1, log_elevation_sd = 1, log_slope = 1)
)

build_axis <- function(df, axis_name, candidates, anchor_weights) {
  available <- candidates[candidates %in% names(df)]
  available <- available[vapply(df[available], function(x) sum(is.finite(x)) >= 0.90 * nrow(df) && stats::sd(x, na.rm = TRUE) > 0, logical(1))]
  if (length(available) < 2) stop("Axis ", axis_name, " has fewer than two usable variables: ", paste(available, collapse = ", "))
  x <- as.data.frame(df[, available, drop = FALSE])
  medians <- vapply(x, median, numeric(1), na.rm = TRUE)
  missing_counts <- vapply(x, function(z) sum(!is.finite(z)), integer(1))
  for (nm in names(x)) x[[nm]][!is.finite(x[[nm]])] <- medians[[nm]]
  pca <- prcomp(x, center = TRUE, scale. = TRUE)
  score <- as.numeric(pca$x[, 1])
  loadings <- pca$rotation[, 1]
  anchor_vars <- intersect(names(anchor_weights), available)
  if (length(anchor_vars)) {
    z <- scale(x[, anchor_vars, drop = FALSE])
    anchor <- as.numeric(z %*% anchor_weights[anchor_vars])
    if (is.finite(cor(score, anchor)) && cor(score, anchor) < 0) {
      score <- -score
      loadings <- -loadings
    }
  }
  explained <- 100 * pca$sdev[1]^2 / sum(pca$sdev^2)
  loading_table <- tibble(
    axis = axis_name,
    source_variable = available,
    loading_pc1 = as.numeric(loadings[available]),
    missing_values_median_imputed = as.integer(missing_counts[available]),
    pc1_variance_explained_pct = explained
  )
  list(score = score, pca = pca, loadings = loading_table, variables = available, explained = explained)
}

axis_objects <- imap(axis_definitions, ~build_axis(env, .y, .x, axis_anchor_definitions[[.y]]))
for (axis_name in names(axis_objects)) env[[axis_name]] <- axis_objects[[axis_name]]$score
axis_loadings <- bind_rows(map(axis_objects, "loadings"))
write_csv_safe(axis_loadings, file.path(output_dir, "12J_environment_axis_loadings.csv"))

# Core predictors deliberately keep the two previously validated atmospheric
# variables separate from the two interpretable landscape PCA axes.
env_predictors <- c("vpd_z", "wind_seasonality_z", "vegetation_axis", "topography_axis")

# Optional localized-moisture controls are tested one at a time and never
# enter the core model automatically.
sensitivity_predictor_sets <- list()
if ("soil_water_mean_frac" %in% names(env) && sum(is.finite(env$soil_water_mean_frac)) >= 0.9 * nrow(env)) {
  env <- median_impute_z(env, "soil_water_mean_frac", "soil_water_z")
  sensitivity_predictor_sets$soil_water <- "soil_water_z"
}
water_candidates <- intersect(c("log_surface_water", "neg_log_distance_water"), names(env))
if (length(water_candidates) >= 2) {
  surface_obj <- build_axis(env, "surface_water_axis", water_candidates, c(log_surface_water = 1, neg_log_distance_water = 1))
  env$surface_water_axis <- surface_obj$score
  axis_loadings <- bind_rows(axis_loadings, surface_obj$loadings)
  write_csv_safe(axis_loadings, file.path(output_dir, "12J_environment_axis_loadings.csv"))
  sensitivity_predictor_sets$surface_water <- "surface_water_axis"
}
if ("oasis_proximity" %in% names(env) && sum(is.finite(env$oasis_proximity)) >= 0.9 * nrow(env)) {
  env <- median_impute_z(env, "oasis_proximity", "oasis_proximity_z")
  sensitivity_predictor_sets$oasis <- "oasis_proximity_z"
}

axis_manifest <- bind_rows(
  tibble(axis = c("vpd_z", "wind_seasonality_z"), n_source_variables = 1L,
         source_variables = c("vpd_mean_kpa", "wind_monthly_sd_ms"),
         pc1_variance_explained_pct = NA_real_, role = "core direct predictor"),
  bind_rows(imap(axis_objects, function(obj, nm) tibble(
    axis = nm, n_source_variables = length(obj$variables),
    source_variables = paste(obj$variables, collapse = ";"),
    pc1_variance_explained_pct = obj$explained, role = "core PCA axis"
  ))),
  if (length(sensitivity_predictor_sets)) tibble(
    axis = unname(unlist(sensitivity_predictor_sets)), n_source_variables = NA_integer_,
    source_variables = names(sensitivity_predictor_sets), pc1_variance_explained_pct = NA_real_,
    role = "optional localized-moisture control"
  ) else tibble()
)
write_csv_safe(axis_manifest, file.path(output_dir, "12J_environment_predictor_manifest.csv"))
score_cols <- unique(c("grid_cell_id", "centroid_latitude", "centroid_longitude", "latitude_band", "easting_km", "northing_km", env_predictors, unname(unlist(sensitivity_predictor_sets))))
write_csv_safe(env %>% select(all_of(intersect(score_cols, names(env)))), file.path(output_dir, "12J_environment_predictor_scores_by_cell.csv"))
log_msg("Core environmental predictors: ", paste(env_predictors, collapse = "; "))
walk2(names(axis_objects), map_dbl(axis_objects, "explained"), ~log_msg(.x, " PC1 variance explained: ", sprintf("%.1f%%", .y)))
if (length(sensitivity_predictor_sets)) log_msg("Optional localized-moisture controls: ", paste(names(sensitivity_predictor_sets), collapse = "; "))

# ---------------------------------------------------------------------------
# Transition definitions
# ---------------------------------------------------------------------------
transition_zones <- tibble::tribble(
  ~break_id, ~break_number, ~break_label, ~anchor_latitude, ~display_min, ~display_max, ~source,
  "la_paz", 1L, "Isthmus of La Paz", 24.20, 23.90, 24.50, "Dolby et al. 2015",
  "loreto", 2L, "Loreto region", 26.00, 25.70, 26.30, "Harrington et al. 2018; Dolby et al. 2015",
  "mid_peninsula", 3L, "Mid-peninsular discontinuity", 27 + 25 / 60, 27 + 20 / 60, 27 + 30 / 60, "Dolby et al. 2015",
  "north_30N", 4L, "30°N climatic transition", 30.00, 29.75, 30.25, "Dolby et al. 2015"
)
write_csv_safe(transition_zones, file.path(output_dir, "12J_transition_zone_definitions.csv"))

# ---------------------------------------------------------------------------
# Assemblage preparation and dissimilarity functions
# ---------------------------------------------------------------------------
analysis_matrices <- list(primary = primary_mat)
if (!is.null(strict_mat)) analysis_matrices$taxonomy_strict <- strict_mat
low_keep <- trait_lookup$confidence_resolved != "LOW"
analysis_matrices$low_confidence_exclusion <- primary_mat[low_keep, , drop = FALSE]
analysis_trait_lookup <- list(
  primary = trait_lookup,
  taxonomy_strict = trait_lookup,
  low_confidence_exclusion = trait_lookup[low_keep, , drop = FALSE]
)
if (is.null(strict_mat)) {
  analysis_matrices$taxonomy_strict <- NULL
  analysis_trait_lookup$taxonomy_strict <- NULL
}

trait_masks_for <- function(traits) {
  list(
    ballooning = which(traits$ballooning_binary_resolved == 1L),
    non_ballooning = which(traits$ballooning_binary_resolved == 0L)
  )
}

simpson_matrix <- function(site_species) {
  site_species <- (site_species > 0) * 1
  richness <- rowSums(site_species)
  shared <- tcrossprod(site_species)
  b <- outer(richness, rep(1, length(richness))) - shared
  c <- t(b)
  numerator <- pmin(b, c)
  denominator <- shared + numerator
  out <- matrix(NA_real_, nrow = nrow(site_species), ncol = nrow(site_species), dimnames = list(rownames(site_species), rownames(site_species)))
  ok <- denominator > 0
  out[ok] <- numerator[ok] / denominator[ok]
  diag(out) <- 0
  out
}

jaccard_matrix <- function(site_species) {
  site_species <- (site_species > 0) * 1
  richness <- rowSums(site_species)
  shared <- tcrossprod(site_species)
  b <- outer(richness, rep(1, length(richness))) - shared
  c <- t(b)
  denominator <- shared + b + c
  out <- matrix(NA_real_, nrow = nrow(site_species), ncol = nrow(site_species), dimnames = list(rownames(site_species), rownames(site_species)))
  ok <- denominator > 0
  out[ok] <- (b[ok] + c[ok]) / denominator[ok]
  diag(out) <- 0
  out
}

make_site_data <- function(matrix_genus_cell, traits, trait_class, min_richness, common_both = FALSE) {
  masks <- trait_masks_for(traits)
  if (common_both) {
    bmat <- t(matrix_genus_cell[masks$ballooning, , drop = FALSE])
    nmat <- t(matrix_genus_cell[masks$non_ballooning, , drop = FALSE])
    eligible <- rowSums(bmat) >= min_richness & rowSums(nmat) >= min_richness
    selected <- if (trait_class == "ballooning") bmat else nmat
  } else {
    selected <- t(matrix_genus_cell[masks[[trait_class]], , drop = FALSE])
    eligible <- rowSums(selected) >= min_richness
  }
  selected <- selected[eligible, , drop = FALSE]
  selected
}

build_pair_table <- function(dissim, site_df, predictors) {
  predictors <- clean_predictors(predictors)
  missing_site_predictors <- setdiff(predictors, names(site_df))
  if (length(missing_site_predictors)) {
    stop("Site table lacks requested GDM predictors: ", paste(missing_site_predictors, collapse = ", "))
  }
  sites <- intersect(rownames(dissim), site_df$grid_cell_id)
  if (length(sites) < 5) stop("Too few sites to build site-pair table")
  dissim <- dissim[sites, sites, drop = FALSE]
  site_df <- site_df[match(sites, site_df$grid_cell_id), , drop = FALSE]
  cmb <- utils::combn(seq_along(sites), 2)
  i <- cmb[1, ]
  j <- cmb[2, ]
  response <- dissim[cbind(i, j)]
  keep <- is.finite(response)
  i <- i[keep]
  j <- j[keep]
  response <- response[keep]
  sp <- data.frame(
    distance = response,
    weights = 1,
    s1.xCoord = site_df$easting_km[i],
    s1.yCoord = site_df$northing_km[i],
    s2.xCoord = site_df$easting_km[j],
    s2.yCoord = site_df$northing_km[j],
    check.names = FALSE
  )
  for (p in predictors) sp[[paste0("s1.", p)]] <- site_df[[p]][i]
  for (p in predictors) sp[[paste0("s2.", p)]] <- site_df[[p]][j]
  class(sp) <- c("gdmData", "data.frame")
  meta <- tibble(
    pair_id = seq_along(i),
    site1 = sites[i],
    site2 = sites[j],
    band1 = as.character(site_df$latitude_band[i]),
    band2 = as.character(site_df$latitude_band[j]),
    latitude1 = site_df$centroid_latitude[i],
    latitude2 = site_df$centroid_latitude[j]
  )
  list(sp = sp, meta = meta, site_df = site_df)
}

clean_predictors <- function(predictors) {
  predictors <- as.character(predictors)
  predictors <- predictors[!is.na(predictors) & nzchar(predictors)]
  unique(predictors)
}

subset_sp_predictors <- function(sp, predictors) {
  predictors <- clean_predictors(predictors)
  required <- c("distance", "weights", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord")
  predictor_cols <- if (length(predictors)) {
    c(paste0("s1.", predictors), paste0("s2.", predictors))
  } else {
    character(0)
  }
  cols <- c(required, predictor_cols)
  missing <- setdiff(cols, names(sp))
  if (length(missing)) stop("Missing site-pair columns: ", paste(missing, collapse = ", "))
  out <- as.data.frame(sp[, cols, drop = FALSE], check.names = FALSE)
  class(out) <- c("gdmData", "data.frame")
  out
}

is_valid_gdm_fit <- function(x) {
  if (is.null(x) || !is.list(x) || !inherits(x, "gdm")) return(FALSE)
  required <- c("nulldeviance", "gdmdeviance", "explained", "intercept")
  all(required %in% names(x))
}

fit_gdm_safe <- function(sp, predictors, label) {
  predictors <- clean_predictors(predictors)
  dat <- subset_sp_predictors(sp, predictors)
  warning_messages <- character()
  fit <- withCallingHandlers(
    tryCatch(gdm::gdm(dat, geo = TRUE), error = function(e) e),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (length(warning_messages)) {
    for (msg in unique(warning_messages)) log_msg("GDM warning [", label, "]: ", msg)
  }
  if (inherits(fit, "error")) {
    log_msg("GDM failed [", label, "]: ", conditionMessage(fit))
    return(NULL)
  }
  # gdm() normally returns a list-class gdm object, but it can return NULL when
  # no spline coefficient is estimable. Protect against any other non-model
  # return value so one failed candidate does not halt the full model set.
  if (!is_valid_gdm_fit(fit)) {
    log_msg(
      "GDM produced no usable model [", label, "]. Returned type/class: ",
      typeof(fit), " / ", paste(class(fit), collapse = ",")
    )
    return(NULL)
  }
  fit
}

safe_model_scalar <- function(model, component) {
  if (!is_valid_gdm_fit(model) || !component %in% names(model)) return(NA_real_)
  value <- suppressWarnings(as.numeric(model[[component]]))
  if (!length(value) || !is.finite(value[[1]])) return(NA_real_)
  value[[1]]
}

extract_spline_effect <- function(model, predictor) {
  if (!is_valid_gdm_fit(model)) return(NA_real_)
  ext <- tryCatch(gdm::isplineExtract(model), error = function(e) NULL)
  if (is.null(ext) || length(ext) < 2 || is.null(colnames(ext[[2]])) || !predictor %in% colnames(ext[[2]])) return(0)
  y <- ext[[2]][, predictor]
  if (!any(is.finite(y))) return(0)
  max(y, na.rm = TRUE) - min(y, na.rm = TRUE)
}

model_row <- function(model, model_name, n_sites, n_pairs) {
  tibble(
    model = model_name,
    n_sites = n_sites,
    n_pairs = n_pairs,
    model_fitted = is_valid_gdm_fit(model),
    null_deviance = safe_model_scalar(model, "nulldeviance"),
    fitted_deviance = safe_model_scalar(model, "gdmdeviance"),
    deviance_explained_pct = safe_model_scalar(model, "explained"),
    intercept = safe_model_scalar(model, "intercept")
  )
}

cv_metrics <- function(obs, pred, baseline) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  if (length(obs) < MIN_TEST_PAIRS) return(tibble(n_test_pairs = length(obs), rmse = NA_real_, mae = NA_real_, correlation = NA_real_, test_deviance_explained_pct = NA_real_))
  rmse <- sqrt(mean((obs - pred)^2))
  mae <- mean(abs(obs - pred))
  correlation <- suppressWarnings(cor(obs, pred))
  null_pred <- rep(baseline, length(obs))
  dev_model <- tryCatch(gdm::calculate.gdm.deviance(pred, obs), error = function(e) NA_real_)
  dev_null <- tryCatch(gdm::calculate.gdm.deviance(null_pred, obs), error = function(e) NA_real_)
  dev_exp <- if (is.finite(dev_model) && is.finite(dev_null) && dev_null > 0) 100 * (1 - dev_model / dev_null) else NA_real_
  tibble(n_test_pairs = length(obs), rmse = rmse, mae = mae, correlation = correlation, test_deviance_explained_pct = dev_exp)
}

blocked_cv <- function(pair_object, predictor_sets, folds) {
  sp <- pair_object$sp
  meta <- pair_object$meta
  rows <- list()
  predictions <- list()
  for (fold in folds) {
    train_idx <- meta$band1 != fold & meta$band2 != fold
    test_idx <- meta$band1 == fold & meta$band2 == fold
    if (sum(train_idx) < 30 || sum(test_idx) < MIN_TEST_PAIRS) next
    for (model_name in names(predictor_sets)) {
      predictors <- predictor_sets[[model_name]]
      train_sp <- sp[train_idx, , drop = FALSE]
      test_sp <- sp[test_idx, , drop = FALSE]
      fit <- fit_gdm_safe(train_sp, predictors, paste0("CV ", fold, " ", model_name))
      if (!is_valid_gdm_fit(fit)) next
      predictors <- clean_predictors(predictors)
      train_dat <- subset_sp_predictors(train_sp, predictors)
      test_dat <- subset_sp_predictors(test_sp, predictors)
      pred <- tryCatch(as.numeric(stats::predict(fit, test_dat)), error = function(e) rep(NA_real_, nrow(test_dat)))
      metrics <- cv_metrics(test_sp$distance, pred, mean(train_sp$distance, na.rm = TRUE))
      rows[[length(rows) + 1]] <- bind_cols(tibble(fold = fold, model = model_name, n_train_pairs = sum(train_idx)), metrics)
      predictions[[length(predictions) + 1]] <- tibble(
        fold = fold, model = model_name,
        pair_id = meta$pair_id[test_idx], site1 = meta$site1[test_idx], site2 = meta$site2[test_idx],
        observed = test_sp$distance, predicted = pred
      )
    }
  }
  list(by_fold = bind_rows(rows), predictions = bind_rows(predictions))
}

summarize_cv <- function(cv_rows) {
  if (!nrow(cv_rows)) return(tibble())
  cv_rows %>% group_by(model) %>% summarise(
    folds = sum(is.finite(rmse)),
    total_test_pairs = sum(n_test_pairs, na.rm = TRUE),
    mean_rmse = weighted.mean(rmse, pmax(1, n_test_pairs), na.rm = TRUE),
    mean_mae = weighted.mean(mae, pmax(1, n_test_pairs), na.rm = TRUE),
    mean_correlation = weighted.mean(correlation, pmax(1, n_test_pairs), na.rm = TRUE),
    mean_test_deviance_explained_pct = weighted.mean(test_deviance_explained_pct, pmax(1, n_test_pairs), na.rm = TRUE),
    .groups = "drop"
  )
}

# ---------------------------------------------------------------------------
# Matched boundary contrasts
# ---------------------------------------------------------------------------
make_all_pair_covariates <- function(cells, sim_balloon, sim_non, anchor) {
  n <- nrow(cells)
  cmb <- utils::combn(seq_len(n), 2)
  i <- cmb[1, ]; j <- cmb[2, ]
  dx <- cells$easting_km[i] - cells$easting_km[j]
  dy <- cells$northing_km[i] - cells$northing_km[j]
  axes <- env_predictors
  env_dist <- sqrt(rowSums((as.matrix(cells[i, axes, drop = FALSE]) - as.matrix(cells[j, axes, drop = FALSE]))^2))
  side1 <- ifelse(cells$centroid_latitude[i] < anchor, "south", "north")
  side2 <- ifelse(cells$centroid_latitude[j] < anchor, "south", "north")
  log_records <- if ("biodiversity_record_count" %in% names(cells)) log1p(pmax(0, safe_numeric(cells$biodiversity_record_count))) else rep(0, n)
  tibble(
    pair_id = seq_along(i),
    site1 = cells$grid_cell_id[i], site2 = cells$grid_cell_id[j],
    side1 = side1, side2 = side2,
    across_boundary = side1 != side2,
    geo_km = sqrt(dx^2 + dy^2),
    env_distance = env_dist,
    midpoint_offset_deg = abs((cells$centroid_latitude[i] + cells$centroid_latitude[j]) / 2 - anchor),
    mean_log_records = (log_records[i] + log_records[j]) / 2,
    simpson_ballooning = sim_balloon[cbind(i, j)],
    simpson_non_ballooning = sim_non[cbind(i, j)]
  )
}

standardize_match_covariates <- function(across, controls) {
  vars <- c("geo_km", "env_distance", "midpoint_offset_deg", "mean_log_records")
  combined <- bind_rows(across[, vars], controls[, vars])
  centers <- vapply(combined, mean, numeric(1), na.rm = TRUE)
  scales <- vapply(combined, sd, numeric(1), na.rm = TRUE)
  scales[!is.finite(scales) | scales == 0] <- 1
  for (v in vars) {
    across[[paste0("z_", v)]] <- (across[[v]] - centers[[v]]) / scales[[v]]
    controls[[paste0("z_", v)]] <- (controls[[v]] - centers[[v]]) / scales[[v]]
  }
  list(across = across, controls = controls, centers = centers, scales = scales)
}

greedy_match_once <- function(across, controls, seed_value, profile_name = "standard") {
  profile <- MATCH_PROFILES[MATCH_PROFILES$profile_name == profile_name, , drop = FALSE]
  if (nrow(profile) != 1) stop("Unknown matching profile: ", profile_name)
  set.seed(seed_value)
  if (!nrow(across) || !nrow(controls)) return(tibble())
  vars <- paste0("z_", c("geo_km", "env_distance", "midpoint_offset_deg", "mean_log_records"))
  across$order_random <- runif(nrow(across))
  # Harder pairs first: farthest from the center of candidate covariate space.
  across$difficulty <- rowSums(abs(as.matrix(across[, vars, drop = FALSE])))
  order_idx <- order(-across$difficulty, across$order_random)
  used <- rep(FALSE, nrow(controls))
  matches <- list()
  for (idx in order_idx) {
    a <- across[idx, , drop = FALSE]
    available <- which(!used)
    if (!length(available)) break
    csub <- controls[available, , drop = FALSE]
    no_shared_site <- !(csub$site1 %in% c(a$site1, a$site2) | csub$site2 %in% c(a$site1, a$site2))
    caliper <- abs(csub$z_geo_km - a$z_geo_km) <= profile$geo_caliper &
      abs(csub$z_env_distance - a$z_env_distance) <= profile$env_caliper &
      abs(csub$z_mean_log_records - a$z_mean_log_records) <= profile$records_caliper &
      abs(csub$z_midpoint_offset_deg - a$z_midpoint_offset_deg) <= profile$midpoint_caliper &
      no_shared_site
    eligible <- available[caliper]
    if (!length(eligible)) next
    diffs <- sweep(as.matrix(controls[eligible, vars, drop = FALSE]), 2, as.numeric(unlist(a[1, vars], use.names = FALSE)), "-")
    dist2 <- rowSums(diffs^2)
    best <- eligible[which.min(dist2)]
    used[best] <- TRUE
    matches[[length(matches) + 1]] <- tibble(
      across_pair_id = a$pair_id,
      control_pair_id = controls$pair_id[best],
      across_site1 = a$site1, across_site2 = a$site2,
      control_site1 = controls$site1[best], control_site2 = controls$site2[best],
      match_distance = sqrt(min(dist2)),
      across_geo_km = a$geo_km, control_geo_km = controls$geo_km[best],
      across_env_distance = a$env_distance, control_env_distance = controls$env_distance[best],
      across_midpoint_offset_deg = a$midpoint_offset_deg, control_midpoint_offset_deg = controls$midpoint_offset_deg[best],
      across_mean_log_records = a$mean_log_records, control_mean_log_records = controls$mean_log_records[best],
      across_simpson_ballooning = a$simpson_ballooning,
      control_simpson_ballooning = controls$simpson_ballooning[best],
      across_simpson_non_ballooning = a$simpson_non_ballooning,
      control_simpson_non_ballooning = controls$simpson_non_ballooning[best]
    )
  }
  bind_rows(matches)
}

match_score <- function(matches) {
  if (!nrow(matches)) return(Inf)
  mean(matches$match_distance, na.rm = TRUE) + 5 / sqrt(nrow(matches))
}

contrast_from_matches <- function(matches) {
  if (!nrow(matches)) return(c(ballooning = NA_real_, non_ballooning = NA_real_, differential = NA_real_))
  excess_b <- matches$across_simpson_ballooning - matches$control_simpson_ballooning
  excess_n <- matches$across_simpson_non_ballooning - matches$control_simpson_non_ballooning
  c(ballooning = mean(excess_b, na.rm = TRUE), non_ballooning = mean(excess_n, na.rm = TRUE), differential = mean(excess_n - excess_b, na.rm = TRUE))
}

cell_jackknife <- function(matches) {
  full_est <- contrast_from_matches(matches)
  cells <- unique(c(matches$across_site1, matches$across_site2, matches$control_site1, matches$control_site2))
  if (length(cells) < 8) return(tibble(trait = names(full_est), estimate = as.numeric(full_est), se = NA_real_, lower95 = NA_real_, upper95 = NA_real_, n_jackknife_cells = length(cells)))
  leave_out <- map_dfr(cells, function(cell) {
    keep <- !(matches$across_site1 == cell | matches$across_site2 == cell | matches$control_site1 == cell | matches$control_site2 == cell)
    est <- contrast_from_matches(matches[keep, , drop = FALSE])
    tibble(cell = cell, trait = names(est), value = as.numeric(est))
  })
  leave_out %>% group_by(trait) %>% summarise(
    estimate = full_est[match(first(trait), names(full_est))],
    jackknife_mean = mean(value, na.rm = TRUE),
    se = sqrt((n() - 1) / n() * sum((value - mean(value, na.rm = TRUE))^2, na.rm = TRUE)),
    n_jackknife_cells = n(),
    .groups = "drop"
  ) %>% mutate(lower95 = estimate - 1.96 * se, upper95 = estimate + 1.96 * se)
}

match_balance <- function(across, controls, matches) {
  if (!nrow(matches)) return(tibble())
  vars <- c("geo_km", "env_distance", "midpoint_offset_deg", "mean_log_records")
  rows <- map_dfr(vars, function(v) {
    before_smd <- (mean(across[[v]], na.rm = TRUE) - mean(controls[[v]], na.rm = TRUE)) / sqrt((var(across[[v]], na.rm = TRUE) + var(controls[[v]], na.rm = TRUE)) / 2)
    av <- matches[[paste0("across_", v)]]
    cv <- matches[[paste0("control_", v)]]
    after_smd <- (mean(av, na.rm = TRUE) - mean(cv, na.rm = TRUE)) / sqrt((var(av, na.rm = TRUE) + var(cv, na.rm = TRUE)) / 2)
    tibble(covariate = v, smd_before = before_smd, smd_after = after_smd)
  })
  rows
}

# ---------------------------------------------------------------------------
# Primary GDM workflow
# ---------------------------------------------------------------------------
# env_predictors was defined during environmental preparation.
fold_order <- c("23-24N", "24-26N", "26-28N", "28-30N", "30-32N")
primary_results <- list()
cv_results <- list()
cv_predictions_all <- list()
spline_results <- list()
partition_results <- list()
varimp_index <- list()
model_objects <- list()
feasibility_rows <- list()
matched_rows <- list()
matched_pair_rows <- list()
matched_balance_rows <- list()
match_iteration_rows <- list()

run_primary_boundary <- function(trait_class, zone) {
  traits <- analysis_trait_lookup$primary
  trait_site_species <- make_site_data(primary_mat, traits, trait_class, GDM_MIN_RICHNESS, common_both = FALSE)
  sites <- intersect(rownames(trait_site_species), env$grid_cell_id)
  if (length(sites) < 10) {
    log_msg("Insufficient sites for primary GDM [", trait_class, " / ", zone$break_id, "]: ", length(sites))
    primary_results[[paste(trait_class, zone$break_id, sep = "__")]] <<- tibble(
      trait_class = trait_class, break_id = zone$break_id, break_label = zone$break_label,
      model = c("geography", "geography_environment", "geography_environment_boundary"),
      n_sites = length(sites), n_pairs = NA_integer_, null_deviance = NA_real_, fitted_deviance = NA_real_,
      deviance_explained_pct = NA_real_, intercept = NA_real_
    )
    return(invisible(NULL))
  }
  trait_site_species <- trait_site_species[sites, , drop = FALSE]
  site_df <- env[match(sites, env$grid_cell_id), , drop = FALSE]
  site_df$boundary_side <- as.numeric(site_df$centroid_latitude >= zone$anchor_latitude)
  diss <- simpson_matrix(trait_site_species)
  # Build the site-pair table with every predictor that may be requested by
  # either the core models or one-at-a-time localized-moisture sensitivities.
  # Earlier versions built only the core columns, causing soil-water,
  # surface-water, and oasis sensitivity models to halt when their standardized
  # columns were requested later.
  pair_predictors <- unique(c(
    env_predictors,
    "boundary_side",
    unname(unlist(sensitivity_predictor_sets, use.names = FALSE))
  ))
  pair_predictors <- pair_predictors[pair_predictors %in% names(site_df)]
  pair <- build_pair_table(diss, site_df, pair_predictors)
  predictor_sets <- list(
    geography = character(),
    geography_environment = env_predictors,
    geography_environment_boundary = c(env_predictors, "boundary_side")
  )
  if (length(sensitivity_predictor_sets)) {
    for (control_name in names(sensitivity_predictor_sets)) {
      control_predictor <- sensitivity_predictor_sets[[control_name]]
      if (control_predictor %in% names(site_df) && all(is.finite(site_df[[control_predictor]]))) {
        predictor_sets[[paste0("geography_environment_", control_name)]] <- c(env_predictors, control_predictor)
        predictor_sets[[paste0("geography_environment_", control_name, "_boundary")]] <- c(env_predictors, control_predictor, "boundary_side")
      }
    }
  }
  fits <- map2(predictor_sets, names(predictor_sets), ~fit_gdm_safe(pair$sp, .x, paste(trait_class, zone$break_id, .y)))
  names(fits) <- names(predictor_sets)
  model_tab <- imap_dfr(fits, ~model_row(.x, .y, nrow(site_df), nrow(pair$sp))) %>%
    mutate(trait_class = trait_class, break_id = zone$break_id, break_label = zone$break_label, .before = 1)
  primary_results[[paste(trait_class, zone$break_id, sep = "__")]] <<- model_tab
  model_objects[[paste(trait_class, zone$break_id, sep = "__")]] <<- fits

  base_explained <- model_tab$deviance_explained_pct[model_tab$model == "geography_environment"]
  boundary_explained <- model_tab$deviance_explained_pct[model_tab$model == "geography_environment_boundary"]
  spline_results[[paste(trait_class, zone$break_id, sep = "__")]] <<- tibble(
    trait_class = trait_class, break_id = zone$break_id, break_label = zone$break_label,
    boundary_increment_deviance_pct = boundary_explained - base_explained,
    boundary_spline_effect = extract_spline_effect(fits$geography_environment_boundary, "boundary_side")
  )

  cv <- blocked_cv(pair, predictor_sets[c("geography_environment", "geography_environment_boundary")], fold_order)
  cv_results[[paste(trait_class, zone$break_id, sep = "__")]] <<- cv$by_fold %>% mutate(trait_class = trait_class, break_id = zone$break_id, break_label = zone$break_label, .before = 1)
  cv_predictions_all[[paste(trait_class, zone$break_id, sep = "__")]] <<- cv$predictions %>% mutate(trait_class = trait_class, break_id = zone$break_id, break_label = zone$break_label, .before = 1)

  full_boundary_sp <- subset_sp_predictors(pair$sp, c(env_predictors, "boundary_side"))
  part <- tryCatch(gdm::gdm.partition.deviance(full_boundary_sp, varSets = list(environment = env_predictors, boundary = "boundary_side"), partSpace = TRUE), error = function(e) NULL)
  if (!is.null(part)) {
    part_df <- as.data.frame(part)
    if (!is.null(rownames(part_df))) part_df <- tibble::rownames_to_column(part_df, "partition")
    partition_results[[paste(trait_class, zone$break_id, sep = "__")]] <<- part_df %>% mutate(trait_class = trait_class, break_id = zone$break_id, break_label = zone$break_label, .before = 1)
  }

  if (run_mode != "audit" && n_perm > 0) {
    vi <- tryCatch(
      gdm::gdm.varImp(
        full_boundary_sp, geo = TRUE, predSelect = FALSE, nPerm = n_perm,
        parallel = cores > 1, cores = cores,
        sampleSites = 1,
        sampleSitePairs = if (run_mode == "quick") 0.5 else 1
      ),
      error = function(e) {
        log_msg("gdm.varImp failed [", trait_class, " / ", zone$break_id, "]: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(vi)) {
      vi_prefix <- paste0("12J_varImp_", trait_class, "_", zone$break_id)
      saveRDS(vi, file.path(model_dir, paste0(vi_prefix, ".rds")))
      for (k in seq_along(vi)) capture_table(vi[[k]], paste0(vi_prefix, "_component", k), output_dir)
      varimp_index[[paste(trait_class, zone$break_id, sep = "__")]] <<- tibble(trait_class = trait_class, break_id = zone$break_id, component_count = length(vi), rds_path = file.path(model_dir, paste0(vi_prefix, ".rds")))
    }
  }
  invisible(NULL)
}

if (run_mode != "audit") {
  for (trait_class in c("ballooning", "non_ballooning")) {
    for (z in seq_len(nrow(transition_zones))) run_primary_boundary(trait_class, transition_zones[z, , drop = FALSE])
  }
}

# ---------------------------------------------------------------------------
# Feasibility audit and matched contrasts using a common cell set
# ---------------------------------------------------------------------------
run_matched_for_dataset <- function(analysis_name, matrix_genus_cell, traits, min_richness, zone,
                                    window_degrees, profile_name, n_iterations = n_match_iterations,
                                    save_pairs = FALSE) {
  bmat <- make_site_data(matrix_genus_cell, traits, "ballooning", min_richness, common_both = TRUE)
  nmat <- make_site_data(matrix_genus_cell, traits, "non_ballooning", min_richness, common_both = TRUE)
  sites <- Reduce(intersect, list(rownames(bmat), rownames(nmat), env$grid_cell_id))
  bmat <- bmat[sites, , drop = FALSE]
  nmat <- nmat[sites, , drop = FALSE]
  cells <- env[match(sites, env$grid_cell_id), , drop = FALSE]
  local <- abs(cells$centroid_latitude - zone$anchor_latitude) <= window_degrees
  cells <- cells[local, , drop = FALSE]
  bmat <- bmat[cells$grid_cell_id, , drop = FALSE]
  nmat <- nmat[cells$grid_cell_id, , drop = FALSE]
  south_n <- sum(cells$centroid_latitude < zone$anchor_latitude)
  north_n <- sum(cells$centroid_latitude >= zone$anchor_latitude)
  med_b <- if (nrow(bmat)) median(rowSums(bmat)) else NA_real_
  med_n <- if (nrow(nmat)) median(rowSums(nmat)) else NA_real_
  base_cols <- list(analysis = analysis_name, min_richness = min_richness,
                    window_degrees = window_degrees, match_profile = profile_name,
                    break_id = zone$break_id, break_label = zone$break_label)
  if (nrow(cells) < 4 || south_n < 2 || north_n < 2) {
    audit <- as_tibble(c(base_cols, list(
      cells_south = south_n, cells_north = north_n, candidate_across_pairs = 0L,
      candidate_control_pairs = 0L, matched_pairs = 0L, unique_cells_matched = 0L,
      median_ballooning_richness = med_b, median_non_ballooning_richness = med_n,
      max_abs_smd_after = NA_real_, feasible = FALSE, reason = "insufficient cells")))
    return(list(audit = audit, contrasts = tibble(), pairs = tibble(), balance = tibble(), iterations = tibble()))
  }
  sim_b <- simpson_matrix(bmat)
  sim_n <- simpson_matrix(nmat)
  all_pairs <- make_all_pair_covariates(cells, sim_b, sim_n, zone$anchor_latitude) %>%
    filter(is.finite(simpson_ballooning), is.finite(simpson_non_ballooning), geo_km <= MAX_LOCAL_PAIR_KM)
  across <- all_pairs %>% filter(across_boundary)
  controls <- all_pairs %>% filter(!across_boundary)
  if (!nrow(across) || !nrow(controls)) {
    audit <- as_tibble(c(base_cols, list(
      cells_south = south_n, cells_north = north_n, candidate_across_pairs = nrow(across),
      candidate_control_pairs = nrow(controls), matched_pairs = 0L, unique_cells_matched = 0L,
      median_ballooning_richness = med_b, median_non_ballooning_richness = med_n,
      max_abs_smd_after = NA_real_, feasible = FALSE, reason = "no usable across-boundary or control pairs")))
    return(list(audit = audit, contrasts = tibble(), pairs = tibble(), balance = tibble(), iterations = tibble()))
  }
  std <- standardize_match_covariates(across, controls)
  across <- std$across; controls <- std$controls
  best <- tibble(); best_score <- Inf
  iteration_summary <- vector("list", n_iterations)
  for (iter in seq_len(n_iterations)) {
    matches <- greedy_match_once(across, controls,
      seed + iter + zone$break_number * 10000 + min_richness * 100000 + round(window_degrees * 1000),
      profile_name = profile_name)
    score <- match_score(matches)
    est <- contrast_from_matches(matches)
    iteration_summary[[iter]] <- tibble(
      analysis = analysis_name, min_richness = min_richness, window_degrees = window_degrees,
      match_profile = profile_name, break_id = zone$break_id, iteration = iter,
      matched_pairs = nrow(matches), match_score = score,
      excess_ballooning = est[["ballooning"]], excess_non_ballooning = est[["non_ballooning"]],
      differential_non_minus_balloon = est[["differential"]]
    )
    if (is.finite(score) && score < best_score) { best <- matches; best_score <- score }
  }
  unique_matched <- if (nrow(best)) unique(c(best$across_site1, best$across_site2, best$control_site1, best$control_site2)) else character()
  balance <- if (nrow(best)) match_balance(across, controls, best) else tibble()
  max_abs_smd <- if (nrow(balance) && any(is.finite(balance$smd_after))) max(abs(balance$smd_after), na.rm = TRUE) else NA_real_
  feasible <- south_n >= MIN_CELLS_PER_SIDE && north_n >= MIN_CELLS_PER_SIDE &&
    nrow(best) >= MIN_MATCHED_PAIRS && length(unique_matched) >= MIN_UNIQUE_MATCHED_CELLS &&
    med_b >= MIN_MEDIAN_RICHNESS && med_n >= MIN_MEDIAN_RICHNESS &&
    is.finite(max_abs_smd) && max_abs_smd <= MAX_MATCHED_ABS_SMD
  reason_bits <- c()
  if (south_n < MIN_CELLS_PER_SIDE || north_n < MIN_CELLS_PER_SIDE) reason_bits <- c(reason_bits, "too few cells on one side")
  if (nrow(best) < MIN_MATCHED_PAIRS) reason_bits <- c(reason_bits, "too few matched pairs")
  if (length(unique_matched) < MIN_UNIQUE_MATCHED_CELLS) reason_bits <- c(reason_bits, "too few unique matched cells")
  if (!is.finite(med_b) || !is.finite(med_n) || med_b < MIN_MEDIAN_RICHNESS || med_n < MIN_MEDIAN_RICHNESS) reason_bits <- c(reason_bits, "median trait richness below safeguard")
  if (!is.finite(max_abs_smd) || max_abs_smd > MAX_MATCHED_ABS_SMD) reason_bits <- c(reason_bits, "post-match balance above safeguard")
  audit <- as_tibble(c(base_cols, list(
    cells_south = south_n, cells_north = north_n,
    candidate_across_pairs = nrow(across), candidate_control_pairs = nrow(controls),
    matched_pairs = nrow(best), unique_cells_matched = length(unique_matched),
    median_ballooning_richness = med_b, median_non_ballooning_richness = med_n,
    max_abs_smd_after = max_abs_smd, feasible = feasible,
    reason = ifelse(feasible, "feasible", paste(reason_bits, collapse = "; ")))))
  contrasts <- if (nrow(best)) cell_jackknife(best) %>% mutate(
    analysis = analysis_name, min_richness = min_richness, window_degrees = window_degrees,
    match_profile = profile_name, break_id = zone$break_id, break_label = zone$break_label,
    matched_pairs = nrow(best), unique_cells_matched = length(unique_matched),
    max_abs_smd_after = max_abs_smd, feasible = feasible, .before = 1
  ) else tibble()
  balance <- if (nrow(balance)) balance %>% mutate(
    analysis = analysis_name, min_richness = min_richness, window_degrees = window_degrees,
    match_profile = profile_name, break_id = zone$break_id, break_label = zone$break_label, .before = 1
  ) else tibble()
  if (save_pairs && nrow(best)) best <- best %>% mutate(
    analysis = analysis_name, min_richness = min_richness, window_degrees = window_degrees,
    match_profile = profile_name, break_id = zone$break_id, break_label = zone$break_label, .before = 1
  )
  list(audit = audit, contrasts = contrasts, pairs = if (save_pairs) best else tibble(),
       balance = balance, iterations = bind_rows(iteration_summary))
}

# Expanded primary feasibility grid. This is a supporting analysis; the
# landscape-wide GDM remains the primary inferential framework.
audit_iterations <- if (run_mode == "audit") min(n_match_iterations, 50L) else n_match_iterations
for (threshold in MATCH_THRESHOLDS) {
  for (window_degrees in MATCH_WINDOWS_DEGREES) {
    for (profile_name in MATCH_PROFILES$profile_name) {
      for (z in seq_len(nrow(transition_zones))) {
        zone <- transition_zones[z, , drop = FALSE]
        matched <- run_matched_for_dataset(
          "primary", analysis_matrices$primary, analysis_trait_lookup$primary,
          threshold, zone, window_degrees, profile_name,
          n_iterations = audit_iterations, save_pairs = FALSE
        )
        feasibility_rows[[length(feasibility_rows) + 1]] <- matched$audit
        if (nrow(matched$contrasts)) matched_rows[[length(matched_rows) + 1]] <- matched$contrasts
        if (nrow(matched$balance)) matched_balance_rows[[length(matched_balance_rows) + 1]] <- matched$balance
        if (nrow(matched$iterations)) match_iteration_rows[[length(match_iteration_rows) + 1]] <- matched$iterations
      }
    }
  }
}

# Taxonomy and LOW-confidence checks use the moderate predeclared setting and
# are run only outside audit mode to keep the audit fast.
if (run_mode != "audit") {
  for (analysis_name in setdiff(names(analysis_matrices), "primary")) {
    matrix_use <- analysis_matrices[[analysis_name]]
    traits_use <- analysis_trait_lookup[[analysis_name]]
    for (z in seq_len(nrow(transition_zones))) {
      zone <- transition_zones[z, , drop = FALSE]
      matched <- run_matched_for_dataset(
        analysis_name, matrix_use, traits_use, 2L, zone, 1.00, "standard",
        n_iterations = n_match_iterations, save_pairs = FALSE
      )
      feasibility_rows[[length(feasibility_rows) + 1]] <- matched$audit
      if (nrow(matched$contrasts)) matched_rows[[length(matched_rows) + 1]] <- matched$contrasts
      if (nrow(matched$balance)) matched_balance_rows[[length(matched_balance_rows) + 1]] <- matched$balance
      if (nrow(matched$iterations)) match_iteration_rows[[length(match_iteration_rows) + 1]] <- matched$iterations
    }
  }
}

feasibility <- bind_rows(feasibility_rows)
matched_contrasts <- bind_rows(matched_rows)
matched_pairs <- bind_rows(matched_pair_rows)
matched_balance <- bind_rows(matched_balance_rows)
match_iterations <- bind_rows(match_iteration_rows)

# Select one configuration per transition using an explicit hierarchy:
# feasible first, then the highest richness threshold, narrowest window,
# standard matching, more matched pairs, and better balance.
selected_match_configs <- feasibility %>%
  filter(analysis == "primary") %>%
  mutate(profile_rank = ifelse(match_profile == "standard", 0L, 1L),
         median_trait_richness = pmin(median_ballooning_richness, median_non_ballooning_richness)) %>%
  group_by(break_id, break_label) %>%
  group_modify(~{
    supported <- .x %>% filter(feasible)
    if (nrow(supported)) {
      supported %>% arrange(desc(min_richness), window_degrees, profile_rank,
                            desc(matched_pairs), max_abs_smd_after) %>% slice(1)
    } else {
      .x %>% arrange(desc(matched_pairs), desc(unique_cells_matched),
                     desc(median_trait_richness), desc(min_richness),
                     window_degrees, profile_rank, max_abs_smd_after) %>% slice(1)
    }
  }) %>%
  ungroup() %>%
  mutate(selection_status = ifelse(feasible, "supported matched configuration", "best available exploratory configuration")) %>%
  select(-profile_rank, -median_trait_richness)

selected_matched_contrasts <- if (nrow(matched_contrasts)) {
  matched_contrasts %>% inner_join(
    selected_match_configs %>% select(break_id, min_richness, window_degrees, match_profile, selection_status),
    by = c("break_id", "min_richness", "window_degrees", "match_profile")
  )
} else tibble()

write_csv_safe(feasibility, file.path(output_dir, "12J_feasibility_audit.csv"))
write_csv_safe(selected_match_configs, file.path(output_dir, "12J_selected_match_configuration.csv"))
write_csv_safe(matched_contrasts, file.path(output_dir, "12J_matched_boundary_contrasts_all_configurations.csv"))
write_csv_safe(selected_matched_contrasts, file.path(output_dir, "12J_selected_matched_boundary_contrasts.csv"))
write_csv_safe(matched_pairs, file.path(output_dir, "12J_matched_pairs_primary.csv"))
write_csv_safe(matched_balance, file.path(output_dir, "12J_matched_pair_balance.csv"))
write_csv_safe(match_iterations, file.path(output_dir, "12J_matching_iteration_sensitivity.csv"))

if (run_mode == "audit") {
  log_msg("Audit mode requested; GDM fitting skipped after feasibility and environmental-predictor construction.")
}

primary_model_comparison <- bind_rows(primary_results)
primary_spline_summary <- bind_rows(spline_results)
primary_cv_by_fold <- bind_rows(cv_results)
primary_cv_predictions <- bind_rows(cv_predictions_all)
primary_cv_summary <- if (nrow(primary_cv_by_fold)) primary_cv_by_fold %>% group_by(trait_class, break_id, break_label) %>% group_modify(~summarize_cv(.x)) %>% ungroup() else tibble()

if (nrow(primary_model_comparison)) write_csv_safe(primary_model_comparison, file.path(output_dir, "12J_primary_gdm_model_comparison.csv"))
if (nrow(primary_model_comparison)) {
  landscape_model_comparison <- primary_model_comparison %>%
    filter(model %in% c("geography", "geography_environment")) %>%
    group_by(trait_class, model) %>%
    summarise(across_boundary_models = n(), n_sites = first(n_sites), n_pairs = first(n_pairs),
              deviance_explained_pct = median(deviance_explained_pct, na.rm = TRUE), .groups = "drop")
  write_csv_safe(landscape_model_comparison, file.path(output_dir, "12J_landscape_gdm_model_comparison.csv"))
}
if (nrow(primary_spline_summary)) write_csv_safe(primary_spline_summary, file.path(output_dir, "12J_primary_gdm_boundary_increment.csv"))
if (nrow(primary_cv_by_fold)) write_csv_safe(primary_cv_by_fold, file.path(output_dir, "12J_primary_gdm_cross_validation_by_fold.csv"))
if (nrow(primary_cv_predictions)) write_csv_safe(primary_cv_predictions, file.path(output_dir, "12J_primary_gdm_cross_validation_predictions.csv"))
if (nrow(primary_cv_summary)) write_csv_safe(primary_cv_summary, file.path(output_dir, "12J_primary_gdm_cross_validation_summary.csv"))
if (length(partition_results)) write_csv_safe(bind_rows(partition_results), file.path(output_dir, "12J_primary_gdm_deviance_partition.csv"))
if (length(varimp_index)) write_csv_safe(bind_rows(varimp_index), file.path(output_dir, "12J_gdm_variable_importance_index.csv"))
if (length(model_objects)) saveRDS(model_objects, file.path(model_dir, "12J_primary_gdm_models.rds"))

# CV increment table.
cv_increment <- tibble()
if (nrow(primary_cv_summary)) {
  cv_wide <- primary_cv_summary %>% select(trait_class, break_id, break_label, model, mean_rmse, mean_mae, mean_test_deviance_explained_pct) %>%
    pivot_wider(names_from = model, values_from = c(mean_rmse, mean_mae, mean_test_deviance_explained_pct))
  cv_increment <- cv_wide %>% transmute(
    trait_class, break_id, break_label,
    rmse_improvement_boundary = mean_rmse_geography_environment - mean_rmse_geography_environment_boundary,
    mae_improvement_boundary = mean_mae_geography_environment - mean_mae_geography_environment_boundary,
    test_deviance_gain_boundary_pct = mean_test_deviance_explained_pct_geography_environment_boundary - mean_test_deviance_explained_pct_geography_environment
  )
  write_csv_safe(cv_increment, file.path(output_dir, "12J_primary_gdm_cross_validation_boundary_increment.csv"))
}

# Optional localized-moisture control sensitivity. Each control is added one
# at a time so it cannot dominate or obscure the core environmental model.
control_sensitivity <- tibble()
if (nrow(primary_model_comparison) && length(sensitivity_predictor_sets)) {
  rows <- list()
  for (control_name in names(sensitivity_predictor_sets)) {
    base_name <- paste0("geography_environment_", control_name)
    boundary_name <- paste0(base_name, "_boundary")
    if (all(c(base_name, boundary_name) %in% primary_model_comparison$model)) {
      tmp <- primary_model_comparison %>%
        select(trait_class, break_id, break_label, model, deviance_explained_pct) %>%
        filter(model %in% c("geography_environment", "geography_environment_boundary", base_name, boundary_name)) %>%
        mutate(model = case_when(
          model == base_name ~ "control_model",
          model == boundary_name ~ "control_boundary_model",
          TRUE ~ model
        )) %>%
        pivot_wider(names_from = model, values_from = deviance_explained_pct) %>%
        mutate(
          control = control_name,
          boundary_increment_core = geography_environment_boundary - geography_environment,
          boundary_increment_with_control = control_boundary_model - control_model,
          change_in_boundary_increment = boundary_increment_with_control - boundary_increment_core,
          control_increment_without_boundary = control_model - geography_environment
        )
      rows[[length(rows) + 1]] <- tmp
    }
  }
  control_sensitivity <- bind_rows(rows)
  if (nrow(control_sensitivity)) write_csv_safe(control_sensitivity, file.path(output_dir, "12J_local_moisture_control_sensitivity.csv"))
}

# Spider-only feasibility sensitivity.
spider_only_status <- tibble()
if (any(tolower(trait_lookup$order_resolved) == "araneae")) {
  spider_traits <- trait_lookup[tolower(trait_lookup$order_resolved) == "araneae", , drop = FALSE]
  spider_mat <- primary_mat[tolower(trait_lookup$order_resolved) == "araneae", , drop = FALSE]
  counts <- table(spider_traits$ballooning_binary_resolved)
  spider_only_status <- tibble(
    ballooning_spider_genera = unname(counts["1"] %||% 0),
    non_ballooning_spider_genera = unname(counts["0"] %||% 0),
    status = ifelse((counts["0"] %||% 0) >= 10, "potentially feasible as a supplementary sensitivity", "insufficient non-ballooning spider genera for a stable full GDM"),
    note = "The master analysis remains taxonomically confounded because ballooning capability is concentrated within Araneae."
  )
} else {
  spider_only_status <- tibble(ballooning_spider_genera = NA_integer_, non_ballooning_spider_genera = NA_integer_, status = "taxonomic order unavailable", note = "No spider-only audit could be performed.")
}
write_csv_safe(spider_only_status, file.path(output_dir, "12J_spider_only_feasibility.csv"))

# ---------------------------------------------------------------------------
# Publication figure
# ---------------------------------------------------------------------------
make_publication_figure <- function() {
  if (!nrow(primary_spline_summary) || !nrow(selected_matched_contrasts)) return(NULL)
  boundary_order <- transition_zones$break_label
  trait_labels <- c(ballooning = "Ballooning-capable", non_ballooning = "Non-ballooning")

  model_plot_data <- primary_model_comparison %>%
    filter(model %in% c("geography", "geography_environment", "geography_environment_boundary")) %>%
    select(trait_class, break_id, break_label, model, deviance_explained_pct) %>%
    pivot_wider(names_from = model, values_from = deviance_explained_pct) %>%
    mutate(
      geography_increment = geography,
      environment_increment = geography_environment - geography,
      boundary_increment = geography_environment_boundary - geography_environment
    ) %>%
    select(trait_class, break_id, break_label, geography_increment, environment_increment, boundary_increment) %>%
    pivot_longer(cols = ends_with("increment"), names_to = "component", values_to = "deviance_increment") %>%
    mutate(
      break_label = factor(break_label, levels = rev(boundary_order)),
      trait_label = factor(trait_labels[trait_class], levels = c("Ballooning-capable", "Non-ballooning")),
      component = factor(component, levels = c("geography_increment", "environment_increment", "boundary_increment"), labels = c("Geographic distance", "Environmental axes", "Boundary crossing"))
    )

  pA <- ggplot(model_plot_data, aes(x = deviance_increment, y = interaction(break_label, trait_label, sep = " — "), fill = component)) +
    geom_col(width = 0.72) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "A  Sequential contributions to Simpson replacement", x = "Additional deviance explained (%)", y = NULL, fill = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold"))

  matched_primary <- selected_matched_contrasts %>%
    filter(trait %in% c("ballooning", "non_ballooning")) %>%
    mutate(
      break_label = factor(break_label, levels = rev(boundary_order)),
      trait_label = factor(trait_labels[trait], levels = c("Ballooning-capable", "Non-ballooning"))
    )
  pB <- ggplot(matched_primary, aes(x = estimate, y = break_label, shape = trait_label)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
    geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.16, position = position_dodge(width = 0.45)) +
    geom_point(size = 2.6, position = position_dodge(width = 0.45)) +
    scale_shape_manual(values = c(16, 17)) +
    labs(title = "B  Excess replacement across matched boundaries", x = "Across-boundary minus matched same-side Simpson turnover", y = NULL, shape = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold"))

  differential_primary <- selected_matched_contrasts %>%
    filter(trait == "differential") %>%
    mutate(break_label = factor(break_label, levels = rev(boundary_order)))
  pC <- ggplot(differential_primary, aes(x = estimate, y = break_label)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
    geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.16) +
    geom_point(size = 2.8) +
    labs(title = "C  Trait difference in residual boundary effect", subtitle = "Positive values indicate a stronger effect among non-ballooning genera", x = "Non-ballooning minus ballooning excess replacement", y = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 8.5))

  combined <- pA / (pB | pC) + plot_layout(heights = c(1.05, 1))
  paths <- c(
    png = file.path(figure_dir, "Figure_3_environment_boundary_trait_turnover.png"),
    pdf = file.path(figure_dir, "Figure_3_environment_boundary_trait_turnover.pdf"),
    svg = file.path(figure_dir, "Figure_3_environment_boundary_trait_turnover.svg"),
    tif = file.path(figure_dir, "Figure_3_environment_boundary_trait_turnover.tif")
  )
  ggsave(paths[["png"]], combined, width = 12.5, height = 9.2, dpi = 400, bg = "white")
  ggsave(paths[["pdf"]], combined, width = 12.5, height = 9.2, device = grDevices::pdf, bg = "white")
  ggsave(paths[["svg"]], combined, width = 12.5, height = 9.2, bg = "white")
  ggsave(paths[["tif"]], combined, width = 12.5, height = 9.2, dpi = 600, compression = "lzw", bg = "white")
  paths
}

figure_paths <- tryCatch(make_publication_figure(), error = function(e) {
  log_msg("Publication figure failed: ", conditionMessage(e))
  NULL
})
if (!is.null(figure_paths)) walk(figure_paths, ~log_msg("Figure written: ", .x))

# Axis loadings figure.
p_axis <- axis_loadings %>% filter(is.finite(loading_pc1)) %>% mutate(axis = factor(axis, levels = rev(unique(axis)))) %>%
  ggplot(aes(x = loading_pc1, y = reorder(source_variable, loading_pc1))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_point() +
  facet_wrap(~axis, scales = "free_y", ncol = 1) +
  labs(title = "Environmental-axis loadings", x = "PC1 loading", y = NULL) +
  theme_bw(base_size = 10) + theme(panel.grid.major.y = element_blank())
ggsave(file.path(figure_dir, "Figure_S12J_environment_axis_loadings.png"), p_axis, width = 8, height = 8, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Recommendation, manuscript text, caption, and validation
# ---------------------------------------------------------------------------
recommendation_lines <- c(
  "STEP 12J MODEL RECOMMENDATION",
  "=============================",
  "",
  "Primary response: Simpson genus replacement.",
  paste0("Landscape-wide GDM minimum richness: ", GDM_MIN_RICHNESS, " genera within the modeled trait class per cell."),
  "Matched boundary configurations were selected from thresholds 1–3, windows 0.75–1.25 degrees, and standard/relaxed calipers using predeclared safeguards.",
  "Published transitions were tested one at a time.",
  ""
)

if (nrow(primary_spline_summary)) {
  recommendation_table <- primary_spline_summary %>%
    left_join(cv_increment, by = c("trait_class", "break_id", "break_label")) %>%
    left_join(selected_matched_contrasts %>% filter(trait == "differential") %>%
      select(break_id, matched_differential = estimate, matched_lower95 = lower95,
             matched_upper95 = upper95, selection_status), by = "break_id")
  write_csv_safe(recommendation_table, file.path(output_dir, "12J_model_decision_table.csv"))
  for (z in seq_len(nrow(transition_zones))) {
    zone <- transition_zones[z, ]
    sub <- recommendation_table %>% filter(break_id == zone$break_id)
    non <- sub %>% filter(trait_class == "non_ballooning")
    bal <- sub %>% filter(trait_class == "ballooning")
    recommendation_lines <- c(recommendation_lines,
      paste0(zone$break_label, ":"),
      paste0("  Boundary deviance increment — ballooning: ", sprintf("%.2f", bal$boundary_increment_deviance_pct %||% NA_real_), "; non-ballooning: ", sprintf("%.2f", non$boundary_increment_deviance_pct %||% NA_real_), " percentage points."),
      paste0("  Cross-validated RMSE improvement — ballooning: ", sprintf("%.4f", bal$rmse_improvement_boundary %||% NA_real_), "; non-ballooning: ", sprintf("%.4f", non$rmse_improvement_boundary %||% NA_real_), "."),
      paste0("  Matched non-minus-balloon differential: ", sprintf("%.3f", unique(sub$matched_differential) %||% NA_real_), " (95% cell-jackknife CI ", sprintf("%.3f", unique(sub$matched_lower95) %||% NA_real_), " to ", sprintf("%.3f", unique(sub$matched_upper95) %||% NA_real_), ")."),
      ""
    )
  }
  supported <- recommendation_table %>% filter(is.finite(boundary_increment_deviance_pct), boundary_increment_deviance_pct > 0, is.finite(rmse_improvement_boundary), rmse_improvement_boundary > 0)
  trait_supported <- recommendation_table %>% distinct(break_id, break_label, matched_differential, matched_lower95, matched_upper95) %>% filter(is.finite(matched_lower95), matched_lower95 > 0)
  recommendation_lines <- c(recommendation_lines,
    if (nrow(supported)) paste0("Boundaries with positive in-sample and cross-validated additions for at least one trait class: ", paste(unique(supported$break_label), collapse = "; "), ".") else "No transition added stable predictive information beyond geography and environment.",
    if (nrow(trait_supported)) paste0("Transitions with a clearly stronger matched residual effect among non-ballooning genera: ", paste(unique(trait_supported$break_label), collapse = "; "), ".") else "No transition showed a clearly stronger matched residual effect among non-ballooning genera.",
    "",
    "Interpretation rule:",
    "  A historical-boundary interpretation is strongest only when boundary crossing adds predictive value beyond geography and environmental axes and the matched across-boundary contrast is positive. A positive non-minus-ballooning differential indicates that dispersal class may moderate boundary permeability. These results remain community-level associations rather than proof of vicariance or speciation."
  )
} else {
  recommendation_lines <- c(recommendation_lines, "Audit mode: no GDM recommendation was generated.")
}
writeLines(recommendation_lines, file.path(output_dir, "12J_model_recommendation.txt"))

caption <- c(
  "Figure 3. Geographic, environmental, and transition-zone contributions to trait-stratified arachnid genus replacement across the Baja California Peninsula.",
  "(A) Sequential deviance explained by geographic distance, four core environmental predictors (vapor-pressure deficit, wind seasonality, vegetation structure, and topographic heterogeneity), and crossing each independently defined transition zone in generalized dissimilarity models fitted separately to ballooning-capable and non-ballooning assemblages. Simpson turnover was used to isolate genus replacement from richness imbalance. Sequential components are descriptive increments and should not be interpreted as fully independent variance fractions. (B) Excess Simpson replacement across each transition relative to geographically and environmentally matched same-side cell pairs. Points show cell-jackknife estimates with 95% intervals. (C) Difference in matched excess replacement between non-ballooning and ballooning assemblages; positive values indicate a stronger residual transition-zone effect among non-ballooning genera. Models were fitted at the 25-km cell grain, and published transitions were tested one at a time."
)
writeLines(caption, file.path(output_dir, "Figure_3_publication_caption.txt"))

manuscript_template <- c(
  "STEP 12J MANUSCRIPT REPLACEMENT GUIDE",
  "======================================",
  "",
  "KEEP:",
  "  Figure 1 — equal-cell genus accumulation/rarefaction.",
  "  Figure 2 — equal-cell adjacent-band Jaccard dissimilarity map.",
  "",
  "REPLACE:",
  "  Delete the original latitude-band ballooning GLM, old prediction heat map, Fisher exact test, seven-fold increase claim, and 26–28°N dispersal-hotspot interpretation.",
  "  Replace them with the Step 12J trait-stratified GDM and matched-boundary analysis and the Step 12J v2 landscape-wide trait-stratified GDM outputs and Figure 3.",
  "",
  "METHODS CORE SENTENCE:",
  "  Generalized dissimilarity models were fitted separately to ballooning-capable and non-ballooning genus assemblages to quantify Simpson genus replacement as a nonlinear function of geographic distance, four core environmental predictors (vapor-pressure deficit, wind seasonality, vegetation structure, and topographic heterogeneity), and crossing independently defined biogeographic transition zones. Each transition was tested in a separate model, and model transferability was assessed using leave-one-latitude-band-out validation. Matched across-boundary and same-side cell pairs provided a complementary estimate of excess boundary-associated replacement after controlling for geographic distance, environmental difference, and sampling intensity.",
  "",
  "RESULTS ORDER:",
  "  1. Dataset and equal-cell richness.",
  "  2. Adjacent-band Jaccard turnover and the strongest interior transition.",
  "  3. GDM contributions of geography and environmental axes.",
  "  4. Incremental and cross-validated effects of individual transition zones.",
  "  5. Matched trait-dependent excess replacement and sensitivities.",
  "",
  "CLAIM LIMITS:",
  "  Do not describe the results as proof of vicariance, hybrid zones, gene flow, ecological speciation, or direct ballooning events. Ballooning classification is concentrated within Araneae and remains partly confounded with higher taxonomy."
)
writeLines(manuscript_template, file.path(output_dir, "12J_manuscript_replacement_guide.txt"))

validation <- tibble(
  check = c(
    "primary_incidence_has_more_than_200_genera",
    "environment_table_has_at_least_150_cells",
    "no_fesa_in_primary_matrix",
    "both_trait_classes_present",
    "four_core_environmental_predictors_constructed",
    "environment_axes_have_finite_scores",
    "four_transition_zones_defined",
    "feasibility_audit_written",
    "primary_matched_results_written",
    "primary_gdm_models_completed_or_audit_mode",
    "cross_validation_completed_or_audit_mode",
    "publication_figure_written_or_audit_mode",
    "model_recommendation_written"
  ),
  passed = c(
    nrow(primary_mat) > 200,
    nrow(env) >= 150,
    !"fesa" %in% tolower(rownames(primary_mat)),
    all(c(0L, 1L) %in% unique(trait_lookup$ballooning_binary_resolved)),
    all(env_predictors %in% names(env)) && length(env_predictors) == 4,
    all(vapply(env[env_predictors], function(x) all(is.finite(x)), logical(1))),
    nrow(transition_zones) == 4,
    file.exists(file.path(output_dir, "12J_feasibility_audit.csv")),
    file.exists(file.path(output_dir, "12J_selected_matched_boundary_contrasts.csv")),
    run_mode == "audit" || sum(is.finite(primary_model_comparison$deviance_explained_pct)) >= 24,
    run_mode == "audit" || n_distinct(paste(primary_cv_by_fold$trait_class, primary_cv_by_fold$break_id, primary_cv_by_fold$model)) >= 16,
    run_mode == "audit" || (!is.null(figure_paths) && all(file.exists(figure_paths))),
    file.exists(file.path(output_dir, "12J_model_recommendation.txt"))
  ),
  detail = c(
    as.character(nrow(primary_mat)),
    as.character(nrow(env)),
    ifelse("fesa" %in% tolower(rownames(primary_mat)), "Fesa present", "Fesa absent"),
    paste(table(trait_lookup$ballooning_binary_resolved), collapse = "/"),
    paste(env_predictors, collapse = ";"),
    paste(vapply(env[env_predictors], function(x) sum(is.finite(x)), integer(1)), collapse = "/"),
    paste(transition_zones$break_id, collapse = ";"),
    paste0("rows=", nrow(feasibility)),
    paste0("selected_rows=", nrow(selected_matched_contrasts)),
    paste0("rows=", nrow(primary_model_comparison)),
    paste0("rows=", nrow(primary_cv_by_fold)),
    ifelse(is.null(figure_paths), "not written", paste(basename(figure_paths), collapse = ";")),
    file.path(output_dir, "12J_model_recommendation.txt")
  )
)
write_csv_safe(validation, file.path(output_dir, "12J_validation.csv"))

log_msg("Primary genera: ", nrow(primary_mat))
log_msg("Environmental cells: ", nrow(env))
log_msg("Matched feasibility grid rows: ", nrow(feasibility))
log_msg("Selected matched configurations: ", nrow(selected_match_configs), "; feasible: ", sum(selected_match_configs$feasible))
log_msg("Primary GDM comparison rows: ", nrow(primary_model_comparison))
log_msg("Validation checks passed: ", sum(validation$passed), "/", nrow(validation))
if (!all(validation$passed)) {
  failed <- validation$check[!validation$passed]
  stop("Step 12J validation failed: ", paste(failed, collapse = "; "))
}
log_msg("STEP 12J MASTER TRAIT-STRATIFIED GDM COMPLETED SUCCESSFULLY")
