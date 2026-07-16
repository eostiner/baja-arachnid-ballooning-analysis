#!/usr/bin/env Rscript

# STEP 12F — Prepare environmental predictor scores for Spatial+
# This script extracts the score-construction component formerly embedded in
# exploratory Step 12J. It makes the retained Spatial+ analysis independent of
# the supplementary GDM workflow and does not change the score definitions.

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(path.expand(args[1]), mustWork = TRUE) else stop(
  "Usage: Rscript 12F_prepare_environment_predictor_scores.R <project_root>"
)

version <- "12F_v1_2026-07-16"
step12c_dir <- file.path(project_root, "04_analysis", "12C_cell_environment_model_table")
out_dir <- file.path(project_root, "04_analysis", "12F_environment_predictor_scores")
archive_root <- file.path(project_root, "08_archive")
input_path <- file.path(step12c_dir, "12C_primary_glm_candidate_table.csv")
if (!file.exists(input_path)) stop("Missing Step 12C primary table: ", input_path)

if (dir.exists(out_dir)) {
  dir.create(archive_root, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(archive_root, paste0("12F_environment_predictor_scores_", format(Sys.time(), "%Y%m%dT%H%M%S")))
  if (!file.rename(out_dir, dest)) stop("Could not archive previous Step 12F output.")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_con <- file(file.path(out_dir, "12F_analysis_log.txt"), open = "wt")
on.exit(close(log_con), add = TRUE)
log_msg <- function(...) { x <- paste0(...); cat(x, "\n"); writeLines(x, log_con); flush(log_con) }
log_msg("STEP 12F STARTED")
log_msg("Version: ", version)
log_msg("Input: ", input_path)

d <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)
required_spatial <- c("grid_cell_id", "centroid_latitude", "centroid_longitude", "latitude_band", "easting_km", "northing_km")
missing_spatial <- setdiff(required_spatial, names(d))
if (length(missing_spatial)) stop("Step 12C table lacks spatial fields: ", paste(missing_spatial, collapse = "; "))
if (anyDuplicated(d$grid_cell_id)) stop("Duplicate grid_cell_id in Step 12C primary table.")

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))
clamp01 <- function(x) pmin(1, pmax(0, x))
for (nm in setdiff(required_spatial, c("grid_cell_id", "latitude_band"))) d[[nm]] <- safe_numeric(d[[nm]])

numeric_candidates <- c(
  "vpd_mean_kpa", "wind_monthly_sd_ms", "evi_mean", "evi_monthly_amplitude",
  "evi_interannual_cv_pct", "lc_shrub_savanna_prop", "lc_grassland_prop",
  "lc_barren_sparse_prop", "relief_5km_m", "elevation_sd_5km_m", "slope_deg"
)
for (nm in intersect(numeric_candidates, names(d))) d[[nm]] <- safe_numeric(d[[nm]])
if (!"vpd_mean_kpa" %in% names(d)) stop("Required predictor missing: vpd_mean_kpa")
if (!"wind_monthly_sd_ms" %in% names(d)) stop("Required predictor missing: wind_monthly_sd_ms")

if ("lc_shrub_savanna_prop" %in% names(d)) d$asin_shrub_savanna <- asin(sqrt(clamp01(d$lc_shrub_savanna_prop)))
if ("lc_grassland_prop" %in% names(d)) d$asin_grassland <- asin(sqrt(clamp01(d$lc_grassland_prop)))
if ("lc_barren_sparse_prop" %in% names(d)) d$asin_barren_sparse <- asin(sqrt(clamp01(d$lc_barren_sparse_prop)))
if ("evi_interannual_cv_pct" %in% names(d)) d$log_evi_interannual_cv <- log1p(pmax(0, d$evi_interannual_cv_pct))
if ("relief_5km_m" %in% names(d)) d$log_relief <- log1p(pmax(0, d$relief_5km_m))
if ("elevation_sd_5km_m" %in% names(d)) d$log_elevation_sd <- log1p(pmax(0, d$elevation_sd_5km_m))
if ("slope_deg" %in% names(d)) d$log_slope <- log1p(pmax(0, d$slope_deg))

median_impute_z <- function(x, label) {
  x <- safe_numeric(x); med <- median(x, na.rm = TRUE)
  if (!is.finite(med)) stop("No finite values for ", label)
  missing_n <- sum(!is.finite(x)); x[!is.finite(x)] <- med
  sx <- sd(x); if (!is.finite(sx) || sx == 0) stop("Zero variance for ", label)
  list(score = as.numeric((x - mean(x)) / sx), median = med, missing = missing_n)
}
vpd <- median_impute_z(d$vpd_mean_kpa, "vpd_mean_kpa")
wind <- median_impute_z(d$wind_monthly_sd_ms, "wind_monthly_sd_ms")
d$vpd_z <- vpd$score
d$wind_seasonality_z <- wind$score

build_axis <- function(df, axis_name, candidates, anchor_weights) {
  available <- candidates[candidates %in% names(df)]
  available <- available[vapply(df[available], function(x) sum(is.finite(x)) >= 0.90 * nrow(df) && sd(x, na.rm = TRUE) > 0, logical(1))]
  if (length(available) < 2L) stop("Axis ", axis_name, " has fewer than two usable variables: ", paste(available, collapse = ", "))
  x <- as.data.frame(df[, available, drop = FALSE])
  medians <- vapply(x, median, numeric(1), na.rm = TRUE)
  missing_counts <- vapply(x, function(z) sum(!is.finite(z)), integer(1))
  for (nm in names(x)) x[[nm]][!is.finite(x[[nm]])] <- medians[[nm]]
  pca <- prcomp(x, center = TRUE, scale. = TRUE)
  score <- as.numeric(pca$x[, 1]); loadings <- pca$rotation[, 1]
  anchor_vars <- intersect(names(anchor_weights), available)
  if (length(anchor_vars)) {
    z <- scale(x[, anchor_vars, drop = FALSE])
    anchor <- as.numeric(z %*% anchor_weights[anchor_vars])
    relation <- suppressWarnings(cor(score, anchor))
    if (is.finite(relation) && relation < 0) { score <- -score; loadings <- -loadings }
  }
  explained <- 100 * pca$sdev[1]^2 / sum(pca$sdev^2)
  loading_table <- data.frame(
    axis = axis_name, source_variable = available,
    loading_pc1 = as.numeric(loadings[available]),
    missing_values_median_imputed = as.integer(missing_counts[available]),
    pc1_variance_explained_pct = explained, stringsAsFactors = FALSE
  )
  list(score = score, loadings = loading_table, variables = available, explained = explained)
}

vegetation <- build_axis(
  d, "vegetation_axis",
  c("evi_mean", "evi_monthly_amplitude", "log_evi_interannual_cv", "asin_shrub_savanna", "asin_grassland", "asin_barren_sparse"),
  c(evi_mean = 1, evi_monthly_amplitude = 0.25, log_evi_interannual_cv = -0.25, asin_shrub_savanna = 1, asin_grassland = 0.5, asin_barren_sparse = -1)
)
topography <- build_axis(
  d, "topography_axis",
  c("log_relief", "log_elevation_sd", "log_slope"),
  c(log_relief = 1, log_elevation_sd = 1, log_slope = 1)
)
d$vegetation_axis <- vegetation$score
d$topography_axis <- topography$score

score_cols <- c(required_spatial, "vpd_z", "wind_seasonality_z", "vegetation_axis", "topography_axis")
scores <- d[, score_cols, drop = FALSE]
if (any(!complete.cases(scores))) stop("Non-finite values remain in the retained Step 12F score table.")
write.csv(scores, file.path(out_dir, "12F_environment_predictor_scores_by_cell.csv"), row.names = FALSE, na = "")
write.csv(rbind(vegetation$loadings, topography$loadings), file.path(out_dir, "12F_environment_axis_loadings.csv"), row.names = FALSE, na = "")
manifest <- data.frame(
  predictor = c("vpd_z", "wind_seasonality_z", "vegetation_axis", "topography_axis"),
  source = c("vpd_mean_kpa", "wind_monthly_sd_ms", paste(vegetation$variables, collapse = ";"), paste(topography$variables, collapse = ";")),
  variance_explained_pct = c(NA, NA, vegetation$explained, topography$explained),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(out_dir, "12F_environment_predictor_manifest.csv"), row.names = FALSE, na = "")
writeLines(c(
  "STEP 12F COMPLETE", paste0("Version: ", version), paste0("Cells: ", nrow(scores)),
  "Spatial+ input: 12F_environment_predictor_scores_by_cell.csv",
  "These score definitions match the retained environmental score construction formerly embedded in exploratory Step 12J."
), file.path(out_dir, "README_12F_OUTPUTS.txt"))
log_msg("Cells written: ", nrow(scores))
log_msg("Vegetation PC1 variance explained: ", sprintf("%.1f%%", vegetation$explained))
log_msg("Topography PC1 variance explained: ", sprintf("%.1f%%", topography$explained))
log_msg("STEP 12F COMPLETE")
