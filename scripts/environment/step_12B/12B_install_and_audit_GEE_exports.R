#!/usr/bin/env Rscript

# ============================================================================
# STEP 12B — INSTALL AND AUDIT GOOGLE EARTH ENGINE EXPORTS
# Version: 12B_v2_2026-07-16
#
# Usage:
#   Rscript 12B_install_and_audit_GEE_exports.R \
#     ~/Desktop/Baja_Ballooning_Pipeline \
#     ~/Downloads/Baja_Ballooning_12B
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) path.expand(args[[1]]) else
  path.expand("~/Desktop/Baja_Ballooning_Pipeline")
source_dir <- if (length(args) >= 2) path.expand(args[[2]]) else
  path.expand("~/Downloads/Baja_Ballooning_12B")

if (!requireNamespace("terra", quietly = TRUE)) {
  stop("Package 'terra' is required. Install it with install.packages('terra').")
}

version <- "12B_v2_2026-07-16"
out_dir <- file.path(project_root, "04_analysis", "12B_fresh_environmental_rasters")
raster_dir <- file.path(
  project_root,
  "ANALYSIS_READY_INPUTS",
  "05_environmental_rasters_12B"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "12B_analysis_log.txt")
log_con <- file(log_file, open = "wt")
on.exit(close(log_con), add = TRUE)
log_msg <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n")
  writeLines(txt, log_con)
  flush(log_con)
}

expected <- list(
  "12B_ERA5Land_climate_2001_2024" = c(
    "tmean_c", "tseason_monthly_sd_c", "tmax_warmest_month_c",
    "precip_annual_mean_mm", "precip_monthly_cv_pct",
    "precip_interannual_cv_pct", "u_wind_mean_ms", "v_wind_mean_ms",
    "wind_speed_mean_ms", "wind_monthly_sd_ms", "vpd_mean_kpa",
    "vpd_driest_month_kpa", "soil_water_mean_frac",
    "soil_water_monthly_sd", "upward_sensible_heat_mean_wm2",
    "upward_sensible_heat_monthly_sd_wm2"
  ),
  "12B_MODIS_EVI_2001_2024" = c(
    "evi_mean", "evi_monthly_amplitude", "evi_interannual_sd",
    "evi_interannual_cv_pct", "evi_valid_observation_count"
  ),
  "12B_MODIS_phenology_2001_2024" = c(
    "greenup_mean_doy", "greenup_interannual_sd_days",
    "phenology_evi_amplitude_mean", "phenology_evi_area_mean",
    "phenology_valid_year_count"
  ),
  "12B_MODIS_landcover_IGBP_mode_2001_2024" = c(
    "lc_igbp_mode"
  ),
  "12B_MODIS_landcover_broad_binary_mode_2001_2024" = c(
    "lc_forest_binary", "lc_shrub_savanna_binary", "lc_grassland_binary",
    "lc_wetland_binary", "lc_cropland_binary", "lc_urban_binary",
    "lc_snow_ice_binary", "lc_barren_sparse_binary", "lc_water_binary"
  ),
  "12B_Copernicus_topography_2024_1" = c(
    "elevation_m", "slope_deg", "elevation_sd_5km_m", "relief_5km_m"
  ),
  "12B_JRC_surface_water_1984_2021" = c(
    "surface_water_occurrence_frac", "surface_water_seasonality_frac"
  )
)

log_msg("STEP 12B INSTALL/AUDIT STARTED")
log_msg("Version: ", version)
log_msg("Project root: ", project_root)
log_msg("Source directory: ", source_dir)
log_msg("Destination raster directory: ", raster_dir)

if (!dir.exists(source_dir)) {
  existing_tifs <- list.files(
    raster_dir,
    pattern = "\\.(tif|tiff)$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  if (length(existing_tifs)) {
    log_msg("Source directory unavailable; auditing existing installed rasters instead.")
    source_dir <- raster_dir
  } else {
    stop("Neither source exports nor installed rasters were found. Missing source: ", source_dir)
  }
}

all_tifs <- list.files(
  source_dir,
  pattern = "\\.(tif|tiff)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

manifest_rows <- list()
validation_rows <- list()

for (prefix in names(expected)) {
  hits <- all_tifs[grepl(prefix, basename(all_tifs), fixed = TRUE)]
  one_file <- length(hits) == 1

  validation_rows[[length(validation_rows) + 1]] <- data.frame(
    check = paste0("one_export_file_", prefix),
    passed = one_file,
    detail = if (length(hits) == 0) "not found" else paste(basename(hits), collapse = "; "),
    stringsAsFactors = FALSE
  )

  if (!one_file) next

  src <- hits[[1]]
  dest <- file.path(raster_dir, paste0(prefix, ".tif"))
  same_file <- identical(
    normalizePath(src, winslash = "/", mustWork = TRUE),
    normalizePath(dest, winslash = "/", mustWork = FALSE)
  )
  if (!same_file) {
    copied <- file.copy(src, dest, overwrite = TRUE)
    if (!copied) stop("Failed to copy: ", src)
  }

  r <- terra::rast(dest)
  actual_names <- names(r)
  expected_names <- expected[[prefix]]
  band_match <- identical(actual_names, expected_names)

  validation_rows[[length(validation_rows) + 1]] <- data.frame(
    check = paste0("band_names_", prefix),
    passed = band_match,
    detail = paste(actual_names, collapse = "; "),
    stringsAsFactors = FALSE
  )

  manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
    prefix = prefix,
    source_path = normalizePath(src, winslash = "/", mustWork = TRUE),
    installed_path = normalizePath(dest, winslash = "/", mustWork = TRUE),
    file_size_bytes = file.info(dest)$size,
    md5 = unname(tools::md5sum(dest)),
    n_layers = terra::nlyr(r),
    layer_names = paste(actual_names, collapse = ";"),
    crs = terra::crs(r, proj = TRUE),
    resolution_x = terra::res(r)[1],
    resolution_y = terra::res(r)[2],
    xmin = terra::xmin(r),
    xmax = terra::xmax(r),
    ymin = terra::ymin(r),
    ymax = terra::ymax(r),
    stringsAsFactors = FALSE
  )

  log_msg("Installed: ", prefix, " — ", terra::nlyr(r), " layers")
}

manifest <- if (length(manifest_rows)) do.call(rbind, manifest_rows) else data.frame()
validation <- do.call(rbind, validation_rows)

write.csv(
  manifest,
  file.path(out_dir, "12B_environmental_raster_manifest.csv"),
  row.names = FALSE
)
write.csv(
  validation,
  file.path(out_dir, "12B_validation.csv"),
  row.names = FALSE
)

n_expected <- length(expected)
n_installed <- if (nrow(manifest)) nrow(manifest) else 0
n_passed <- sum(validation$passed)
n_checks <- nrow(validation)

log_msg("Expected exports: ", n_expected)
log_msg("Installed exports: ", n_installed)
log_msg("Validation checks passed: ", n_passed, "/", n_checks)

if (n_installed != n_expected || !all(validation$passed)) {
  log_msg("STEP 12B INCOMPLETE — inspect 12B_validation.csv")
  quit(save = "no", status = 1)
}

writeLines(
  c(
    "STEP 12B COMPLETED SUCCESSFULLY",
    paste0("Version: ", version),
    paste0("Installed raster directory: ", raster_dir),
    "Next step: Step 12C cell-polygon extraction and verified cell-level model table."
  ),
  file.path(out_dir, "12B_completion.txt")
)

log_msg("STEP 12B COMPLETED SUCCESSFULLY")
log_msg("Next: build Step 12C cell-level environmental table.")
