#!/usr/bin/env Rscript

# =============================================================================
# STEP 12C — BUILD VERIFIED 25-KM CELL ENVIRONMENTAL MODEL TABLE
# Baja Ballooning Publication
# Version: 12C_C3_publication_v6_2026-07-16
#
# Usage:
#   Rscript 12C_build_cell_environment_model_table.R \
#     ~/Desktop/Baja_Ballooning_Pipeline
#
# Purpose:
#   1. Read the frozen 25-km genus × cell incidence matrices and reviewed
#      genus-level ballooning traits.
#   2. Calculate ballooning and non-ballooning genus counts for every occupied
#      25-km cell (primary, taxonomy-strict, and explicit LOW-confidence
#      exclusion sensitivity responses).
#   3. Extract cell-polygon means from the seven freshly installed Step 12B
#      raster products. Broad binary land-cover layers are converted to cell
#      proportions; modal IGBP class is retained as a diagnostic field.
#   4. Join response, sampling, geometry, and environmental information by the
#      frozen grid_cell_id — never by row position.
#   5. Write audited analysis tables and freeze the main model table in
#      ANALYSIS_READY_INPUTS.
#
# Required R package: terra
# Optional R package: jsonlite (for JSON provenance; a text fallback is used)
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

SCRIPT_VERSION <- "12C_C3_publication_v6_2026-07-16"
EXPECTED_CELL_COUNT <- 205L
EXPECTED_GENUS_COUNT <- 267L

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) {
  normalizePath(path.expand(args[[1]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(
    path.expand("~/Desktop/Baja_Ballooning_Pipeline"),
    winslash = "/",
    mustWork = TRUE
  )
}

if (!requireNamespace("terra", quietly = TRUE)) {
  stop("Package 'terra' is required. Install it with install.packages('terra').")
}

# ------------------------------- helpers -----------------------------------

first_existing <- function(paths, label, required = TRUE) {
  paths <- path.expand(paths)
  hit <- paths[file.exists(paths)]
  if (length(hit)) return(normalizePath(hit[[1]], winslash = "/", mustWork = TRUE))
  if (required) {
    stop("Could not find ", label, ". Tried:\n", paste(paths, collapse = "\n"))
  }
  NA_character_
}

find_field <- function(fields, candidates, label, required = TRUE) {
  lower <- tolower(trimws(fields))
  for (candidate in candidates) {
    idx <- match(tolower(candidate), lower)
    if (!is.na(idx)) return(fields[[idx]])
  }
  if (required) {
    stop(
      "Could not identify ", label, ". Tried: ",
      paste(candidates, collapse = ", "),
      ". Available fields: ", paste(fields, collapse = ", ")
    )
  }
  NA_character_
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))
trim_character <- function(x) trimws(as.character(x))

parse_evidence_one <- function(value) {
  s <- toupper(trimws(as.character(value)))
  if (!nzchar(s) || is.na(s)) return(NA_character_)
  hit <- regmatches(s, gregexpr(
    "(?<![A-Z0-9])(D1|D2|D3|D4|N0|C3)(?![A-Z0-9])",
    s, perl = TRUE
  ))[[1]]
  hit <- unique(toupper(hit[nzchar(hit)]))
  if (length(hit) == 1L && hit %in% c("D1", "D2", "D3", "D4", "N0", "C3")) {
    return(hit)
  }
  n <- gsub("[^a-z0-9]+", "", tolower(s))
  if (n %in% c("nonballooning", "fixednonballooning", "referencenonballooning", "noballooning")) return("N0")
  if (n %in% c("c3", "primaryc3", "d1d2d3", "d1tod3")) return("C3")
  if (n %in% c("d4excluded", "excludedd4")) return("D4")
  NA_character_
}

normalize_evidence <- function(x) {
  vapply(x, parse_evidence_one, character(1), USE.NAMES = FALSE)
}

find_evidence_field <- function(df) {
  preferred <- c(
    "evidence_class", "final_evidence_class", "final_evidence_category",
    "evidence_category", "d_level", "dlevel", "trait_class",
    "primary_class", "analysis_class", "ballooning_evidence_tier",
    "ballooning_evidence_category", "final_designation", "designation"
  )
  scores <- rep(-Inf, ncol(df)); names(scores) <- names(df)
  for (nm in names(df)) {
    parsed <- normalize_evidence(df[[nm]])
    frac <- mean(!is.na(parsed))
    classes <- unique(parsed[!is.na(parsed)])
    if (frac >= 0.25 && length(classes) >= 2L) {
      clean <- gsub("[^a-z0-9]+", "", tolower(nm))
      bonus <- if (tolower(nm) %in% tolower(preferred)) 100 else 0
      if (grepl("evidence|tier|class|designation|decision", clean)) bonus <- bonus + 20
      scores[[nm]] <- bonus + 100 * frac + 5 * length(classes)
    }
  }
  if (!any(is.finite(scores))) {
    stop(
      "Trait table lacks an explicit D1/D2/D3/D4/N0 or C3/N0 evidence field. ",
      "Legacy binary fields are intentionally not accepted because they cannot distinguish D4 from N0."
    )
  }
  names(which.max(scores))
}

md5_file <- function(path) unname(tools::md5sum(path))

archive_directory <- function(path, archive_root, prefix) {
  if (!dir.exists(path) || length(list.files(path, all.files = TRUE, no.. = TRUE)) == 0L) {
    return(NA_character_)
  }
  stamp <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid())
  dest <- file.path(archive_root, paste0(prefix, "_", stamp))
  dir.create(archive_root, recursive = TRUE, showWarnings = FALSE)

  # A rename is reliable and fast when source and archive are on the same disk.
  if (isTRUE(file.rename(path, dest))) {
    return(normalizePath(dest, winslash = "/", mustWork = TRUE))
  }

  # Cross-filesystem fallback: copy the directory contents into an existing
  # destination. Copying the source directory itself to a nonexistent target
  # is version-dependent in base R and caused the prior Step 12C failure.
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  items <- list.files(path, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  ok <- if (length(items)) file.copy(items, dest, recursive = TRUE, copy.date = TRUE) else TRUE
  if (!all(ok)) stop("Failed to archive existing output directory: ", path)
  unlink(path, recursive = TRUE, force = TRUE)
  normalizePath(dest, winslash = "/", mustWork = TRUE)
}

read_incidence_matrix <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(x) < 2L || nrow(x) < 1L) stop("Invalid incidence matrix: ", path)
  genera <- trim_character(x[[1]])
  cells <- trim_character(names(x)[-1])
  if (any(!nzchar(genera))) stop("Blank genus in incidence matrix: ", path)
  if (anyDuplicated(tolower(genera))) stop("Duplicate genus in incidence matrix: ", path)
  if (any(!nzchar(cells)) || anyDuplicated(cells)) stop("Invalid cell columns: ", path)
  mat <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  if (any(is.na(mat)) || any(!mat %in% c(0, 1))) {
    stop("Incidence matrix contains values other than 0/1: ", path)
  }
  rownames(mat) <- genera
  colnames(mat) <- cells
  list(path = path, genera = genera, cells = cells, matrix = mat)
}

align_incidence <- function(obj, target_genera, target_cells) {
  genus_idx <- match(tolower(target_genera), tolower(obj$genera))
  cell_idx <- match(target_cells, obj$cells)
  if (anyNA(genus_idx)) {
    stop("Incidence matrix lacks required genera: ",
         paste(target_genera[is.na(genus_idx)][1:min(10, sum(is.na(genus_idx)))], collapse = ", "))
  }
  if (anyNA(cell_idx)) {
    stop("Incidence matrix lacks required cells: ",
         paste(target_cells[is.na(cell_idx)][1:min(10, sum(is.na(cell_idx)))], collapse = ", "))
  }
  mat <- obj$matrix[genus_idx, cell_idx, drop = FALSE]
  rownames(mat) <- target_genera
  colnames(mat) <- target_cells
  mat
}

summarize_numeric <- function(name, values) {
  values <- safe_numeric(values)
  good <- is.finite(values)
  q <- if (any(good)) stats::quantile(values[good], c(0, .25, .5, .75, 1), na.rm = TRUE) else rep(NA_real_, 5)
  data.frame(
    predictor = name,
    n_cells = length(values),
    n_nonmissing = sum(good),
    n_missing = sum(!good),
    missing_percent = 100 * sum(!good) / length(values),
    minimum = unname(q[[1]]),
    q25 = unname(q[[2]]),
    median = unname(q[[3]]),
    mean = if (any(good)) mean(values[good]) else NA_real_,
    standard_deviation = if (sum(good) > 1L) stats::sd(values[good]) else NA_real_,
    q75 = unname(q[[4]]),
    maximum = unname(q[[5]]),
    stringsAsFactors = FALSE
  )
}

log_con <- NULL
log_msg <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n")
  if (!is.null(log_con)) {
    writeLines(txt, log_con)
    flush(log_con)
  }
}

# ------------------------------ directories --------------------------------

analysis_ready <- file.path(project_root, "ANALYSIS_READY_INPUTS")
grid_fallback <- file.path(project_root, "02_data_clean", "08_grid25km_incidence")
trait_fallback <- file.path(project_root, "02_data_clean", "07_final_trait_merge")
raster_dir <- file.path(analysis_ready, "05_environmental_rasters_12B")
out_dir <- file.path(project_root, "04_analysis", "12C_cell_environment_model_table")
frozen_dir <- file.path(analysis_ready, "06_environmental_cell_tables_12C")
archive_root <- file.path(project_root, "08_archive")

archived_output <- archive_directory(out_dir, archive_root, "12C_cell_environment_model_table")
if (dir.exists(frozen_dir) && length(list.files(frozen_dir, all.files = TRUE, no.. = TRUE))) {
  archive_directory(frozen_dir, archive_root, "12C_frozen_environmental_cell_tables")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(frozen_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "12C_analysis_log.txt")
log_con <- file(log_file, open = "wt")
on.exit(if (!is.null(log_con)) close(log_con), add = TRUE)

log_msg("STEP 12C STARTED")
log_msg("Version: ", SCRIPT_VERSION)
log_msg("Project root: ", project_root)
if (!is.na(archived_output)) log_msg("Archived prior Step 12C output: ", archived_output)

# -------------------------------- inputs -----------------------------------

primary_path <- first_existing(c(
  file.path(analysis_ready, "02_incidence_matrices_25km", "10_biodiversity_final_genus_by_grid25km_incidence.csv"),
  file.path(grid_fallback, "10_biodiversity_final_genus_by_grid25km_incidence.csv")
), "primary biodiversity incidence matrix")

strict_path <- first_existing(c(
  file.path(analysis_ready, "02_incidence_matrices_25km", "10_biodiversity_taxonomy_strict_genus_by_grid25km_incidence.csv"),
  file.path(grid_fallback, "10_biodiversity_taxonomy_strict_genus_by_grid25km_incidence.csv")
), "taxonomy-strict incidence matrix")

cell_lookup_path <- first_existing(c(
  file.path(analysis_ready, "04_spatial_reference", "10_common_grid25km_cell_lookup.csv"),
  file.path(grid_fallback, "10_common_grid25km_cell_lookup.csv")
), "25-km cell lookup")

cell_geojson_path <- first_existing(c(
  file.path(analysis_ready, "04_spatial_reference", "10_common_grid25km_cells.geojson"),
  file.path(grid_fallback, "10_common_grid25km_cells.geojson")
), "25-km cell polygon GeoJSON")

trait_path <- first_existing(c(
  file.path(analysis_ready, "03_trait_tables", "07_reviewed_genus_trait_lookup_normalized.csv"),
  file.path(analysis_ready, "03_trait_tables", "07_reviewed_genus_trait_lookup_final.csv"),
  file.path(trait_fallback, "07_reviewed_genus_trait_lookup_final.csv")
), "reviewed genus trait lookup")

expected_rasters <- c(
  era5 = "12B_ERA5Land_climate_2001_2024.tif",
  evi = "12B_MODIS_EVI_2001_2024.tif",
  phenology = "12B_MODIS_phenology_2001_2024.tif",
  landcover_igbp = "12B_MODIS_landcover_IGBP_mode_2001_2024.tif",
  landcover_broad = "12B_MODIS_landcover_broad_binary_mode_2001_2024.tif",
  topography = "12B_Copernicus_topography_2024_1.tif",
  surface_water = "12B_JRC_surface_water_1984_2021.tif"
)
raster_paths <- file.path(raster_dir, expected_rasters)
names(raster_paths) <- names(expected_rasters)
if (any(!file.exists(raster_paths))) {
  stop("Missing Step 12B raster(s):\n", paste(raster_paths[!file.exists(raster_paths)], collapse = "\n"))
}

log_msg("Primary matrix: ", primary_path)
log_msg("Taxonomy-strict matrix: ", strict_path)
log_msg("Cell lookup: ", cell_lookup_path)
log_msg("Cell polygons: ", cell_geojson_path)
log_msg("Trait lookup: ", trait_path)
log_msg("Raster products found: ", length(raster_paths), "/", length(expected_rasters))

# -------------------------- response construction --------------------------

primary_obj <- read_incidence_matrix(primary_path)
strict_obj <- read_incidence_matrix(strict_path)
genera <- primary_obj$genera
cells <- primary_obj$cells

if (any(tolower(genera) == "fesa")) stop("Fesa remains in the primary incidence matrix.")
if (length(genera) != EXPECTED_GENUS_COUNT) {
  log_msg("NOTICE: Primary genus count is ", length(genera), ", not documented ", EXPECTED_GENUS_COUNT, ".")
}
if (length(cells) != EXPECTED_CELL_COUNT) {
  log_msg("NOTICE: Primary occupied-cell count is ", length(cells), ", not documented ", EXPECTED_CELL_COUNT, ".")
}

primary_mat <- primary_obj$matrix
strict_mat <- align_incidence(strict_obj, genera, cells)

traits <- read.csv(trait_path, check.names = FALSE, stringsAsFactors = FALSE)
genus_field <- find_field(names(traits), c("genus", "analysis_genus"), "trait genus")
evidence_field <- find_evidence_field(traits)
confidence_field <- find_field(names(traits), c(
  "final_confidence", "trait_final_confidence", "trait_confidence", "trait_ballooning_confidence"
), "trait confidence field", required = FALSE)

traits$genus_key <- tolower(trim_character(traits[[genus_field]]))
if (anyDuplicated(traits$genus_key[traits$genus_key != ""])) stop("Duplicate genus in trait table.")
trait_idx <- match(tolower(genera), traits$genus_key)
if (anyNA(trait_idx)) stop("Genera missing from trait table: ", paste(genera[is.na(trait_idx)], collapse = ", "))

evidence_class <- normalize_evidence(traits[[evidence_field]][trait_idx])
if (anyNA(evidence_class)) {
  stop(
    "Unresolved D1/D2/D3/D4/N0 or C3/N0 traits for: ",
    paste(genera[is.na(evidence_class)], collapse = ", ")
  )
}
analysis_class <- ifelse(
  evidence_class %in% c("D1", "D2", "D3", "C3"), "C3",
  ifelse(evidence_class == "N0", "N0", "D4_excluded")
)
trait_confidence <- if (!is.na(confidence_field)) {
  toupper(trim_character(traits[[confidence_field]][trait_idx]))
} else {
  rep("UNSPECIFIED", length(genera))
}
trait_confidence[!nzchar(trait_confidence)] <- "UNSPECIFIED"

c3_rows <- analysis_class == "C3"
n0_rows <- analysis_class == "N0"
d4_rows <- analysis_class == "D4_excluded"
low_rows <- trait_confidence == "LOW"
if (any(c3_rows & n0_rows) || any(c3_rows & d4_rows) || any(n0_rows & d4_rows)) {
  stop("Trait classes are not disjoint.")
}
if (!all(c3_rows | n0_rows | d4_rows)) stop("Trait classes do not cover all genera.")

count_response <- function(mat, keep_rows = rep(TRUE, nrow(mat))) {
  keep_c3 <- keep_rows & c3_rows
  keep_n0 <- keep_rows & n0_rows
  keep_d4 <- keep_rows & d4_rows
  balloon <- if (any(keep_c3)) colSums(mat[keep_c3, , drop = FALSE]) else rep(0, ncol(mat))
  nonballoon <- if (any(keep_n0)) colSums(mat[keep_n0, , drop = FALSE]) else rep(0, ncol(mat))
  d4_excluded <- if (any(keep_d4)) colSums(mat[keep_d4, , drop = FALSE]) else rep(0, ncol(mat))
  classified <- balloon + nonballoon
  all_genera <- classified + d4_excluded
  data.frame(
    grid_cell_id = colnames(mat),
    ballooning_genera = as.integer(balloon),
    non_ballooning_genera = as.integer(nonballoon),
    classified_genera = as.integer(classified),
    excluded_D4_genera = as.integer(d4_excluded),
    all_genera = as.integer(all_genera),
    ballooning_proportion = ifelse(classified > 0, balloon / classified, NA_real_),
    stringsAsFactors = FALSE
  )
}

resp_primary <- count_response(primary_mat)
resp_strict <- count_response(strict_mat)
resp_low_excl <- count_response(primary_mat, keep_rows = !low_rows)

names(resp_primary)[-1] <- paste0(names(resp_primary)[-1], "_primary")
names(resp_strict)[-1] <- paste0(names(resp_strict)[-1], "_taxonomy_strict")
names(resp_low_excl)[-1] <- paste0(names(resp_low_excl)[-1], "_low_conf_exclusion")

response <- Reduce(function(x, y) merge(x, y, by = "grid_cell_id", all = TRUE, sort = FALSE),
                   list(resp_primary, resp_strict, resp_low_excl))
response <- response[match(cells, response$grid_cell_id), , drop = FALSE]

# ----------------------------- cell geometry -------------------------------

lookup <- read.csv(cell_lookup_path, check.names = FALSE, stringsAsFactors = FALSE)
lookup_id_field <- find_field(names(lookup), c("grid_cell_id", "cell_id"), "cell lookup ID")
lookup$grid_cell_id <- trim_character(lookup[[lookup_id_field]])
if (anyDuplicated(lookup$grid_cell_id)) stop("Duplicate grid_cell_id in cell lookup.")
lookup <- lookup[match(cells, lookup$grid_cell_id), , drop = FALSE]
if (anyNA(lookup$grid_cell_id)) stop("Cell lookup does not contain all matrix cells.")

lat_field <- find_field(names(lookup), c("centroid_latitude", "latitude"), "centroid latitude")
lon_field <- find_field(names(lookup), c("centroid_longitude", "longitude"), "centroid longitude")
band_field <- find_field(names(lookup), c("centroid_latitude_band", "latitude_band"), "latitude band")
record_field <- find_field(names(lookup), c("biodiversity_record_count", "occurrence_record_count"), "record count", required = FALSE)
richness_field <- find_field(names(lookup), c("biodiversity_genus_richness", "genus_richness"), "lookup richness", required = FALSE)

cell_info <- data.frame(
  grid_cell_id = lookup$grid_cell_id,
  grid_cell_order = if ("grid_cell_order" %in% names(lookup)) safe_numeric(lookup$grid_cell_order) else seq_len(nrow(lookup)),
  grid_row = if ("grid_row" %in% names(lookup)) safe_numeric(lookup$grid_row) else NA_real_,
  grid_column = if ("grid_column" %in% names(lookup)) safe_numeric(lookup$grid_column) else NA_real_,
  centroid_latitude = safe_numeric(lookup[[lat_field]]),
  centroid_longitude = safe_numeric(lookup[[lon_field]]),
  latitude_band = trim_character(lookup[[band_field]]),
  centroid_x_m = if ("centroid_x_m" %in% names(lookup)) safe_numeric(lookup$centroid_x_m) else NA_real_,
  centroid_y_m = if ("centroid_y_m" %in% names(lookup)) safe_numeric(lookup$centroid_y_m) else NA_real_,
  biodiversity_record_count = if (!is.na(record_field)) safe_numeric(lookup[[record_field]]) else NA_real_,
  lookup_biodiversity_genus_richness = if (!is.na(richness_field)) safe_numeric(lookup[[richness_field]]) else NA_real_,
  stringsAsFactors = FALSE
)

lat_center <- mean(cell_info$centroid_latitude, na.rm = TRUE)
cell_info$latitude_centered <- cell_info$centroid_latitude - lat_center
cell_info$latitude_centered_sq <- cell_info$latitude_centered^2
cell_info$latitude_raw_sq <- cell_info$centroid_latitude^2
cell_info$easting_km <- cell_info$centroid_x_m / 1000
cell_info$northing_km <- cell_info$centroid_y_m / 1000

cell_polygons_all <- terra::vect(cell_geojson_path)
geom_id_field <- find_field(names(cell_polygons_all), c("grid_cell_id", "cell_id"), "polygon grid-cell ID")
geom_ids_raw <- trim_character(cell_polygons_all[[geom_id_field]])

# First attempt the explicit frozen-ID join. Some GeoJSON writers can alter or
# omit plus signs in identifiers such as BJA25K_C+0013_R-0020. If any required
# IDs fail to match, fall back to a spatial join using the verified lookup
# centroids. This preserves the exact frozen polygon geometry while avoiding
# a fragile dependence on GeoJSON attribute-string formatting.
geom_idx <- match(cells, geom_ids_raw)
polygon_match_method <- rep("exact_grid_cell_id", length(cells))

if (anyNA(geom_idx)) {
  log_msg(
    "NOTICE: Exact polygon ID matching failed for ", sum(is.na(geom_idx)),
    " of ", length(cells), " occupied cells; attempting centroid-to-polygon spatial match."
  )

  if (any(!is.finite(cell_info$centroid_longitude)) || any(!is.finite(cell_info$centroid_latitude))) {
    stop("Cannot perform centroid-to-polygon fallback because lookup centroids are missing.")
  }

  cell_polygons_all$polygon_index_12C <- seq_len(nrow(cell_polygons_all))
  lookup_points <- terra::vect(
    data.frame(
      point_order_12C = seq_along(cells),
      longitude = cell_info$centroid_longitude,
      latitude = cell_info$centroid_latitude
    ),
    geom = c("longitude", "latitude"),
    crs = "EPSG:4326"
  )
  if (!terra::same.crs(lookup_points, cell_polygons_all)) {
    lookup_points <- terra::project(lookup_points, terra::crs(cell_polygons_all, proj = TRUE))
  }

  spatial_hits <- terra::intersect(lookup_points, cell_polygons_all)
  # terra::as.data.frame.SpatVector accepts geom = NULL, "WKT", "HEX", or
  # "XY"; logical FALSE is invalid in current terra versions. Omitting the
  # geometry argument returns the intersected feature attributes, which are
  # exactly the point and polygon indices needed for this join.
  spatial_df <- as.data.frame(spatial_hits)
  required_hit_fields <- c("point_order_12C", "polygon_index_12C")
  if (!all(required_hit_fields %in% names(spatial_df))) {
    stop(
      "Centroid-to-polygon fallback did not retain expected index fields. Available fields: ",
      paste(names(spatial_df), collapse = ", ")
    )
  }

  hit_counts <- table(spatial_df$point_order_12C)
  bad_points <- setdiff(seq_along(cells), as.integer(names(hit_counts)[hit_counts == 1L]))
  if (length(bad_points)) {
    stop(
      "Centroid-to-polygon fallback did not find exactly one polygon for required cells: ",
      paste(cells[bad_points], collapse = ", ")
    )
  }

  spatial_df <- spatial_df[order(spatial_df$point_order_12C), , drop = FALSE]
  geom_idx <- as.integer(spatial_df$polygon_index_12C)
  polygon_match_method[] <- "centroid_spatial_join"
}

cell_polygons <- cell_polygons_all[geom_idx, ]
cell_polygons$grid_cell_id <- cells
cell_polygons <- terra::makeValid(cell_polygons)
cell_info$cell_area_km2 <- as.numeric(terra::expanse(cell_polygons, unit = "km"))

polygon_id_audit <- data.frame(
  grid_cell_id = cells,
  original_polygon_id = geom_ids_raw[geom_idx],
  match_method = polygon_match_method,
  centroid_longitude = cell_info$centroid_longitude,
  centroid_latitude = cell_info$centroid_latitude,
  cell_area_km2 = cell_info$cell_area_km2,
  stringsAsFactors = FALSE
)
write.csv(
  polygon_id_audit,
  file.path(out_dir, "12C_polygon_geometry_match_audit.csv"),
  row.names = FALSE,
  na = ""
)
log_msg(
  "Cell polygons matched: ", nrow(cell_polygons), "/", length(cells),
  " using ", paste(unique(polygon_match_method), collapse = "; ")
)

# -------------------------- raster polygon extraction ----------------------

extract_stack_mean <- function(path, cells_vect, product_name) {
  r <- terra::rast(path)
  cv <- terra::crs(r, proj = TRUE)
  cells_r <- if (terra::same.crs(cells_vect, r)) cells_vect else terra::project(cells_vect, cv)
  log_msg("Extracting ", product_name, " — ", terra::nlyr(r), " layer(s)")
  values <- terra::extract(r, cells_r, fun = mean, na.rm = TRUE, exact = TRUE, ID = TRUE)
  values <- values[match(seq_len(nrow(cells_r)), values$ID), , drop = FALSE]
  values$ID <- NULL
  valid_r <- terra::ifel(is.na(r), 0, 1)
  valid <- terra::extract(valid_r, cells_r, fun = mean, na.rm = TRUE, exact = TRUE, ID = TRUE)
  valid <- valid[match(seq_len(nrow(cells_r)), valid$ID), , drop = FALSE]
  valid$ID <- NULL
  names(valid) <- paste0(names(valid), "__valid_fraction")
  list(values = values, valid = valid, raster = r, cells = cells_r)
}

extract_modal <- function(path, cells_vect, product_name) {
  r <- terra::rast(path)
  cv <- terra::crs(r, proj = TRUE)
  cells_r <- if (terra::same.crs(cells_vect, r)) cells_vect else terra::project(cells_vect, cv)
  log_msg("Extracting exact area-weighted modal class for ", product_name)

  # terra::extract() does not provide a built-in modal polygon summary, and
  # exact=TRUE only supports mean, sum, min, max, and table. Extract the
  # categorical pixels with their exact polygon-overlap fractions, sum those
  # fractions by class within each polygon, and select the class with the
  # largest covered area. Ties are resolved deterministically to the smaller
  # numeric IGBP code.
  raw_values <- terra::extract(r, cells_r, exact = TRUE, ID = TRUE, raw = TRUE)
  raw_values <- as.data.frame(raw_values, stringsAsFactors = FALSE)

  value_col <- names(r)[1]
  if (!value_col %in% names(raw_values)) {
    candidates <- setdiff(names(raw_values), c("ID", "cell", "fraction", "weight"))
    if (length(candidates) < 1L) {
      stop("Could not identify the IGBP class-value column returned by terra::extract().")
    }
    value_col <- candidates[1]
  }

  fraction_candidates <- intersect(c("fraction", "weight"), names(raw_values))
  if (length(fraction_candidates) < 1L) {
    stop("Exact IGBP extraction did not return a polygon-overlap fraction column.")
  }
  fraction_col <- fraction_candidates[1]

  mode_vec <- rep(NA_real_, nrow(cells_r))
  for (i in seq_len(nrow(cells_r))) {
    d <- raw_values[
      raw_values$ID == i &
        !is.na(raw_values[[value_col]]) &
        !is.na(raw_values[[fraction_col]]) &
        raw_values[[fraction_col]] > 0,
      c(value_col, fraction_col),
      drop = FALSE
    ]
    if (nrow(d) < 1L) next

    class_values <- suppressWarnings(as.numeric(as.character(d[[value_col]])))
    cover_values <- suppressWarnings(as.numeric(d[[fraction_col]]))
    keep <- is.finite(class_values) & is.finite(cover_values) & cover_values > 0
    if (!any(keep)) next

    area_by_class <- tapply(cover_values[keep], class_values[keep], sum, na.rm = TRUE)
    max_area <- max(area_by_class, na.rm = TRUE)
    tied_classes <- suppressWarnings(as.numeric(names(area_by_class)[area_by_class == max_area]))
    mode_vec[i] <- min(tied_classes, na.rm = TRUE)
  }

  values <- data.frame(lc_igbp_cell_mode = mode_vec, stringsAsFactors = FALSE)

  valid_r <- terra::ifel(is.na(r), 0, 1)
  valid <- terra::extract(valid_r, cells_r, fun = mean, na.rm = TRUE, exact = TRUE, ID = TRUE)
  valid <- valid[match(seq_len(nrow(cells_r)), valid$ID), , drop = FALSE]
  valid$ID <- NULL
  names(valid) <- "lc_igbp_mode__valid_fraction"

  list(values = values, valid = valid, raster = r, cells = cells_r)
}

environment <- data.frame(grid_cell_id = cells, stringsAsFactors = FALSE)
extraction_manifest <- list()

for (key in c("era5", "evi", "phenology", "topography", "surface_water")) {
  result <- extract_stack_mean(raster_paths[[key]], cell_polygons, key)
  environment <- cbind(environment, result$values, result$valid)
  extraction_manifest[[length(extraction_manifest) + 1L]] <- data.frame(
    product = key,
    raster_path = normalizePath(raster_paths[[key]], winslash = "/", mustWork = TRUE),
    n_layers = terra::nlyr(result$raster),
    layer_names = paste(names(result$raster), collapse = ";"),
    extraction = "exact polygon-weighted mean",
    md5 = md5_file(raster_paths[[key]]),
    stringsAsFactors = FALSE
  )
}

# Broad binary land-cover means become within-cell proportions.
lc_result <- extract_stack_mean(raster_paths[["landcover_broad"]], cell_polygons, "landcover_broad")
lc_values <- lc_result$values
names(lc_values) <- sub("_binary$", "_prop", names(lc_values))
lc_valid <- lc_result$valid
names(lc_valid) <- sub("_binary__valid_fraction$", "_prop__valid_fraction", names(lc_valid))
environment <- cbind(environment, lc_values, lc_valid)
extraction_manifest[[length(extraction_manifest) + 1L]] <- data.frame(
  product = "landcover_broad",
  raster_path = normalizePath(raster_paths[["landcover_broad"]], winslash = "/", mustWork = TRUE),
  n_layers = terra::nlyr(lc_result$raster),
  layer_names = paste(names(lc_result$raster), collapse = ";"),
  extraction = "exact polygon-weighted mean; binary pixels interpreted as cell proportions",
  md5 = md5_file(raster_paths[["landcover_broad"]]),
  stringsAsFactors = FALSE
)

igbp_result <- extract_modal(raster_paths[["landcover_igbp"]], cell_polygons, "landcover_igbp")
environment <- cbind(environment, igbp_result$values, igbp_result$valid)
extraction_manifest[[length(extraction_manifest) + 1L]] <- data.frame(
  product = "landcover_igbp",
  raster_path = normalizePath(raster_paths[["landcover_igbp"]], winslash = "/", mustWork = TRUE),
  n_layers = 1L,
  layer_names = names(igbp_result$raster),
  extraction = "exact area-weighted modal IGBP class within cell",
  md5 = md5_file(raster_paths[["landcover_igbp"]]),
  stringsAsFactors = FALSE
)

# Approximate centroid distance to the nearest MODIS water-class pixel.
water_r <- terra::rast(raster_paths[["landcover_broad"]])[["lc_water_binary"]]
water_mask <- terra::ifel(water_r >= 0.5, 1, NA)
log_msg("Calculating centroid distance to nearest MODIS water-class pixel")
water_distance <- terra::distance(water_mask)
centroid_points <- terra::centroids(
  if (terra::same.crs(cell_polygons, water_distance)) cell_polygons else terra::project(cell_polygons, terra::crs(water_distance, proj = TRUE))
)
distance_values <- terra::extract(water_distance, centroid_points, ID = TRUE)
distance_values <- distance_values[match(seq_len(nrow(centroid_points)), distance_values$ID), , drop = FALSE]
environment$distance_to_modis_water_km <- safe_numeric(distance_values[[2]]) / 1000

# Land-cover checks and approximate land area.
lc_prop_fields <- grep("^lc_.*_prop$", names(environment), value = TRUE)
environment$lc_classified_prop_sum <- rowSums(environment[, lc_prop_fields, drop = FALSE], na.rm = TRUE)
environment$lc_unclassified_prop <- pmax(0, 1 - environment$lc_classified_prop_sum)
if ("lc_water_prop" %in% names(environment)) {
  environment$modis_land_fraction <- pmin(1, pmax(0, 1 - environment$lc_water_prop))
} else {
  environment$modis_land_fraction <- NA_real_
}

# ------------------------------ joins --------------------------------------

model_table <- Reduce(function(x, y) merge(x, y, by = "grid_cell_id", all = TRUE, sort = FALSE),
                      list(cell_info, response, environment))
model_table <- model_table[match(cells, model_table$grid_cell_id), , drop = FALSE]
model_table$estimated_land_area_km2 <- model_table$cell_area_km2 * model_table$modis_land_fraction
model_table$response_denominator_positive <- model_table$classified_genera_primary > 0
model_table$response_denominator_ge5 <- model_table$classified_genera_primary >= 5
model_table$response_denominator_ge10 <- model_table$classified_genera_primary >= 10

core_predictors <- intersect(c(
  "tmean_c", "precip_annual_mean_mm", "wind_speed_mean_ms", "vpd_mean_kpa",
  "soil_water_mean_frac", "evi_mean", "elevation_m", "slope_deg",
  "lc_shrub_savanna_prop", "lc_barren_sparse_prop"
), names(model_table))
model_table$core_environment_complete <- rowSums(is.na(model_table[, core_predictors, drop = FALSE])) == 0L
model_table$land_fraction_ge_0_10 <- model_table$modis_land_fraction >= 0.10
model_table$recommended_primary_model_cell <- with(
  model_table,
  core_environment_complete & response_denominator_positive & land_fraction_ge_0_10
)
model_table$sensitivity_model_cell_ge5 <- with(
  model_table, recommended_primary_model_cell & response_denominator_ge5
)
model_table$sensitivity_model_cell_ge10 <- with(
  model_table, recommended_primary_model_cell & response_denominator_ge10
)

# ------------------------------ summaries ----------------------------------

metadata_fields <- c(
  "grid_cell_order", "grid_row", "grid_column", "centroid_latitude",
  "centroid_longitude", "latitude_band", "centroid_x_m", "centroid_y_m",
  "biodiversity_record_count", "lookup_biodiversity_genus_richness",
  "latitude_centered", "latitude_centered_sq", "latitude_raw_sq",
  "easting_km", "northing_km", "cell_area_km2"
)
response_fields <- grep("_(primary|taxonomy_strict|low_conf_exclusion)$", names(model_table), value = TRUE)
flag_fields <- c(
  "response_denominator_positive", "response_denominator_ge5", "response_denominator_ge10",
  "core_environment_complete", "land_fraction_ge_0_10",
  "recommended_primary_model_cell"
)
valid_fraction_fields <- grep("__valid_fraction$", names(model_table), value = TRUE)
environment_fields <- setdiff(names(environment), "grid_cell_id")

predictor_summary_fields <- setdiff(environment_fields, valid_fraction_fields)
predictor_summary <- do.call(rbind, lapply(predictor_summary_fields, function(nm) {
  summarize_numeric(nm, model_table[[nm]])
}))

missingness <- data.frame(
  field = c(predictor_summary_fields, valid_fraction_fields),
  role = c(rep("environmental_value", length(predictor_summary_fields)), rep("raster_valid_fraction", length(valid_fraction_fields))),
  n_cells = nrow(model_table),
  n_missing = vapply(c(predictor_summary_fields, valid_fraction_fields), function(nm) sum(!is.finite(safe_numeric(model_table[[nm]]))), integer(1)),
  stringsAsFactors = FALSE
)
missingness$missing_percent <- 100 * missingness$n_missing / missingness$n_cells
missingness <- missingness[order(-missingness$missing_percent, missingness$field), ]

response_by_band <- aggregate(
  model_table[, c(
    "ballooning_genera_primary", "non_ballooning_genera_primary",
    "classified_genera_primary", "ballooning_proportion_primary",
    "biodiversity_record_count"
  )],
  by = list(latitude_band = model_table$latitude_band),
  FUN = function(x) c(n = sum(is.finite(x)), mean = mean(x, na.rm = TRUE), median = median(x, na.rm = TRUE), min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE))
)
# Flatten aggregate matrix-columns.
response_band_flat <- data.frame(latitude_band = response_by_band$latitude_band, stringsAsFactors = FALSE)
for (nm in names(response_by_band)[-1]) {
  m <- response_by_band[[nm]]
  for (stat in colnames(m)) response_band_flat[[paste(nm, stat, sep = "_")]] <- m[, stat]
}

trait_summary <- as.data.frame(table(
  evidence_class = evidence_class,
  analysis_class = analysis_class,
  confidence = trait_confidence
), stringsAsFactors = FALSE)

field_dictionary <- data.frame(
  field = names(model_table),
  role = ifelse(names(model_table) == "grid_cell_id", "join_key",
    ifelse(names(model_table) %in% metadata_fields, "cell_metadata",
      ifelse(names(model_table) %in% response_fields, "response_or_sensitivity_response",
        ifelse(names(model_table) %in% valid_fraction_fields, "raster_coverage_diagnostic",
          ifelse(names(model_table) %in% flag_fields, "model_readiness_flag", "environmental_predictor_or_derived_metric")
        )
      )
    )
  ),
  stringsAsFactors = FALSE
)

# ------------------------------ validation ---------------------------------

validation <- list()
add_check <- function(check, passed, detail) {
  validation[[length(validation) + 1L]] <<- data.frame(
    check = check,
    passed = isTRUE(passed),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

add_check("seven_step12B_rasters_present", length(raster_paths) == 7L && all(file.exists(raster_paths)), paste(basename(raster_paths), collapse = "; "))
add_check("unique_primary_cell_ids", !anyDuplicated(cells), length(cells))
add_check("model_table_one_row_per_cell", nrow(model_table) == length(cells) && !anyDuplicated(model_table$grid_cell_id), nrow(model_table))
add_check("model_table_cell_order_matches_incidence", identical(model_table$grid_cell_id, cells), paste(head(model_table$grid_cell_id), collapse = "; "))
add_check("all_traits_resolved", !anyNA(evidence_class), paste(table(analysis_class), collapse = "; "))
add_check("response_counts_sum_primary", all(model_table$ballooning_genera_primary + model_table$non_ballooning_genera_primary == model_table$classified_genera_primary), "ballooning + non-ballooning = classified")
add_check("response_counts_sum_strict", all(model_table$ballooning_genera_taxonomy_strict + model_table$non_ballooning_genera_taxonomy_strict == model_table$classified_genera_taxonomy_strict), "ballooning + non-ballooning = classified")
add_check("lookup_richness_matches_all_primary_genera", all(is.na(model_table$lookup_biodiversity_genus_richness) | model_table$lookup_biodiversity_genus_richness == model_table$all_genera_primary), paste0("mismatches=", sum(!is.na(model_table$lookup_biodiversity_genus_richness) & model_table$lookup_biodiversity_genus_richness != model_table$all_genera_primary)))
add_check("D4_excluded_from_primary_denominator", all(model_table$classified_genera_primary + model_table$excluded_D4_genera_primary == model_table$all_genera_primary), "C3 + N0 + D4 = all genera")
candidate <- model_table$recommended_primary_model_cell %in% TRUE
add_check(
  "primary_candidate_environment_complete",
  sum(candidate) > 0L && all(model_table$core_environment_complete[candidate]),
  paste0(
    sum(candidate), " candidate cells; ",
    sum(!model_table$core_environment_complete, na.rm = TRUE),
    " incomplete cells excluded"
  )
)
add_check("landcover_proportion_sum_reasonable", all(is.na(model_table$lc_classified_prop_sum) | (model_table$lc_classified_prop_sum >= 0.95 & model_table$lc_classified_prop_sum <= 1.05)), paste(range(model_table$lc_classified_prop_sum, na.rm = TRUE), collapse = " to "))
add_check(
  "primary_candidate_denominators_positive",
  sum(candidate) > 0L && all(model_table$classified_genera_primary[candidate] > 0),
  paste0(
    sum(candidate), " candidate cells; ",
    sum(model_table$classified_genera_primary <= 0, na.rm = TRUE),
    " zero-denominator cells excluded"
  )
)
add_check("cell_area_approximately_625_km2", all(model_table$cell_area_km2 > 600 & model_table$cell_area_km2 < 650), paste(round(range(model_table$cell_area_km2), 3), collapse = " to "))
validation_df <- do.call(rbind, validation)

# ------------------------------- outputs -----------------------------------

response_path <- file.path(out_dir, "12C_cell_response_counts.csv")
environment_path <- file.path(out_dir, "12C_cell_environmental_predictors.csv")
model_path <- file.path(out_dir, "12C_cell_environment_model_table.csv")
model_primary_path <- file.path(out_dir, "12C_primary_glm_candidate_table.csv")
model_ge5_path <- file.path(out_dir, "12C_sensitivity_glm_candidate_ge5.csv")
model_ge10_path <- file.path(out_dir, "12C_sensitivity_glm_candidate_ge10.csv")

write.csv(response, response_path, row.names = FALSE, na = "")
write.csv(environment, environment_path, row.names = FALSE, na = "")
write.csv(model_table, model_path, row.names = FALSE, na = "")
write.csv(model_table[model_table$recommended_primary_model_cell, , drop = FALSE], model_primary_path, row.names = FALSE, na = "")
write.csv(model_table[model_table$sensitivity_model_cell_ge5, , drop = FALSE], model_ge5_path, row.names = FALSE, na = "")
write.csv(model_table[model_table$sensitivity_model_cell_ge10, , drop = FALSE], model_ge10_path, row.names = FALSE, na = "")
write.csv(predictor_summary, file.path(out_dir, "12C_predictor_summary.csv"), row.names = FALSE, na = "")
write.csv(missingness, file.path(out_dir, "12C_predictor_missingness.csv"), row.names = FALSE, na = "")
write.csv(response_band_flat, file.path(out_dir, "12C_response_by_latitude_band_summary.csv"), row.names = FALSE, na = "")
write.csv(trait_summary, file.path(out_dir, "12C_trait_classification_summary.csv"), row.names = FALSE, na = "")
write.csv(field_dictionary, file.path(out_dir, "12C_field_dictionary.csv"), row.names = FALSE, na = "")
write.csv(do.call(rbind, extraction_manifest), file.path(out_dir, "12C_extraction_manifest.csv"), row.names = FALSE, na = "")
write.csv(validation_df, file.path(out_dir, "12C_validation.csv"), row.names = FALSE, na = "")

# Geometry-bearing output for mapping/diagnostics.
cell_polygons_out <- cell_polygons
cell_polygons_out$ballooning_genera_primary <- model_table$ballooning_genera_primary
cell_polygons_out$non_ballooning_genera_primary <- model_table$non_ballooning_genera_primary
cell_polygons_out$ballooning_proportion_primary <- model_table$ballooning_proportion_primary
cell_polygons_out$recommended_primary_model_cell <- as.integer(model_table$recommended_primary_model_cell)
terra::writeVector(cell_polygons_out, file.path(out_dir, "12C_occupied_cells_response.geojson"), overwrite = TRUE)

# Freeze the main table and dictionary for subsequent model steps.
file.copy(model_path, file.path(frozen_dir, basename(model_path)), overwrite = TRUE)
file.copy(model_primary_path, file.path(frozen_dir, basename(model_primary_path)), overwrite = TRUE)
file.copy(model_ge5_path, file.path(frozen_dir, basename(model_ge5_path)), overwrite = TRUE)
file.copy(model_ge10_path, file.path(frozen_dir, basename(model_ge10_path)), overwrite = TRUE)
file.copy(file.path(out_dir, "12C_field_dictionary.csv"), file.path(frozen_dir, "12C_field_dictionary.csv"), overwrite = TRUE)
file.copy(file.path(out_dir, "12C_validation.csv"), file.path(frozen_dir, "12C_validation.csv"), overwrite = TRUE)

input_manifest <- data.frame(
  role = c("primary_incidence", "taxonomy_strict_incidence", "cell_lookup", "cell_polygons", "trait_lookup", names(raster_paths)),
  path = c(primary_path, strict_path, cell_lookup_path, cell_geojson_path, trait_path, unname(raster_paths)),
  bytes = file.info(c(primary_path, strict_path, cell_lookup_path, cell_geojson_path, trait_path, unname(raster_paths)))$size,
  md5 = vapply(c(primary_path, strict_path, cell_lookup_path, cell_geojson_path, trait_path, unname(raster_paths)), md5_file, character(1)),
  stringsAsFactors = FALSE
)
write.csv(input_manifest, file.path(out_dir, "12C_input_manifest.csv"), row.names = FALSE)

provenance <- list(
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  script_version = SCRIPT_VERSION,
  project_root = project_root,
  analytical_unit = "one occupied 25-km equal-area grid cell",
  join_key = "grid_cell_id",
  response = "counts of C3 (D1-D3) and fixed N0 genera per cell; D4 excluded from denominator",
  trait_evidence_field = evidence_field,
  trait_evidence_counts = as.list(table(evidence_class)),
  trait_analysis_counts = as.list(table(analysis_class)),
  environmental_extraction = "polygon-weighted cell means; land-cover binary means interpreted as proportions; IGBP class summarized by exact area-weighted modal class",
  latitude_center = lat_center,
  n_cells = nrow(model_table),
  n_genera = length(genera),
  low_confidence_genera_excluded_in_sensitivity = genera[low_rows],
  core_predictors_for_completeness_flag = core_predictors,
  recommended_primary_model_cells = sum(model_table$recommended_primary_model_cell),
  validation_passed = all(validation_df$passed),
  input_manifest = input_manifest
)
if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(provenance, file.path(out_dir, "12C_provenance.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
} else {
  capture.output(dput(provenance), file = file.path(out_dir, "12C_provenance.txt"))
}

readme_lines <- c(
  "STEP 12C — VERIFIED CELL-LEVEL ENVIRONMENTAL MODEL TABLE",
  "========================================================",
  paste0("Rows: ", nrow(model_table), " occupied 25-km cells"),
  paste0("Primary genera: ", length(genera)),
  paste0("Recommended primary candidate cells: ", sum(model_table$recommended_primary_model_cell)),
  paste0("Sensitivity candidate cells (>=5 classified genera): ", sum(model_table$sensitivity_model_cell_ge5)),
  paste0("Sensitivity candidate cells (>=10 classified genera): ", sum(model_table$sensitivity_model_cell_ge10)),
  "",
  "Primary trait definition: C3 = D1 + D2 + D3 versus fixed N0; D4 excluded.",
  paste0("Trait evidence field: ", evidence_field),
  "",
  "Primary response columns:",
  "  ballooning_genera_primary",
  "  non_ballooning_genera_primary",
  "  classified_genera_primary",
  "  ballooning_proportion_primary",
  "",
  "Important interpretation:",
  "  Land-cover *_prop columns are polygon means of binary source pixels and",
  "  therefore estimate within-cell class proportions. distance_to_modis_water_km",
  "  is distance from the cell centroid to the nearest MODIS water-class pixel;",
  "  it is not a pure coastline-only distance.",
  "",
  "Next step:",
  "  Step 12D screens distributions, correlations, collinearity, and candidate",
  "  geography-only, environment-only, and combined model formulas."
)
writeLines(readme_lines, file.path(out_dir, "README_12C_OUTPUTS.txt"))

log_msg("Cells in model table: ", nrow(model_table))
log_msg("Primary genera represented: ", length(genera))
log_msg("Recommended primary candidate cells: ", sum(model_table$recommended_primary_model_cell))
log_msg("Validation checks passed: ", sum(validation_df$passed), "/", nrow(validation_df))

if (!all(validation_df$passed)) {
  log_msg("STEP 12C INCOMPLETE — inspect 12C_validation.csv")
  quit(save = "no", status = 1)
}

writeLines(c(
  "STEP 12C COMPLETED SUCCESSFULLY",
  paste0("Version: ", SCRIPT_VERSION),
  paste0("Cell-level model table: ", model_path),
  paste0("Frozen analysis-ready table: ", file.path(frozen_dir, basename(model_path))),
  "Next: Step 12D predictor screening and candidate model comparison."
), file.path(out_dir, "12C_completion.txt"))

log_msg("STEP 12C COMPLETED SUCCESSFULLY")
log_msg("Next: Step 12D predictor screening and candidate model comparison.")
