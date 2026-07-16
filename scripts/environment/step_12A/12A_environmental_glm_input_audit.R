#!/usr/bin/env Rscript

# Step 12A — Environmental GLM input discovery and audit
# Purpose: identify the current frozen community/trait/spatial/environmental inputs,
#          verify join keys, summarize missingness, and inventory environmental rasters.
# This script does NOT fit a model.

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(path.expand(args[1]), mustWork = TRUE) else stop("Usage: Rscript 12A_environmental_glm_input_audit.R <project_root>")
legacy_root <- if (length(args) >= 2 && nzchar(args[2])) {
  normalizePath(path.expand(args[2]), mustWork = FALSE)
} else {
  NA_character_
}

out_dir <- file.path(project_root, "04_analysis", "12A_environmental_glm_input_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "12A_analysis_log.txt")
log_con <- file(log_file, open = "wt")
on.exit(close(log_con), add = TRUE)
log_msg <- function(...) {
  x <- paste0(...)
  cat(x, "\n")
  cat(x, "\n", file = log_con)
  flush(log_con)
}

version <- "12A_v2_2026-07-13"
log_msg("STEP 12A STARTED")
log_msg("Version: ", version)
log_msg("Project root: ", project_root)
log_msg("Optional legacy search root: ", ifelse(is.na(legacy_root), "DISABLED", legacy_root))

# ---------- helpers ----------
clean_names_simple <- function(x) {
  y <- tolower(x)
  y <- gsub("[^a-z0-9]+", "_", y)
  y <- gsub("^_+|_+$", "", y)
  y
}

read_header <- function(path) {
  ext <- tolower(tools::file_ext(path))
  out <- list(ok = FALSE, nrow = NA_integer_, ncol = NA_integer_, names = character(), error = NA_character_)
  tryCatch({
    if (ext == "csv") {
      d <- utils::read.csv(path, nrows = 2000, check.names = FALSE, stringsAsFactors = FALSE)
      out$ok <- TRUE; out$nrow <- nrow(d); out$ncol <- ncol(d); out$names <- names(d)
    } else if (ext %in% c("tsv", "txt")) {
      d <- utils::read.delim(path, nrows = 2000, check.names = FALSE, stringsAsFactors = FALSE)
      out$ok <- TRUE; out$nrow <- nrow(d); out$ncol <- ncol(d); out$names <- names(d)
    } else if (ext %in% c("xlsx", "xls") && requireNamespace("readxl", quietly = TRUE)) {
      d <- readxl::read_excel(path, n_max = 2000)
      out$ok <- TRUE; out$nrow <- nrow(d); out$ncol <- ncol(d); out$names <- names(d)
    } else {
      out$error <- if (ext %in% c("xlsx", "xls")) "readxl not installed" else "unsupported extension"
    }
  }, error = function(e) out$error <<- conditionMessage(e))
  out
}

score_table <- function(cols, filename) {
  c0 <- clean_names_simple(cols)
  fn <- clean_names_simple(basename(filename))
  has_any <- function(patterns) any(vapply(patterns, function(p) any(grepl(p, c0)), logical(1)))
  has_key <- has_any(c("^site_id$", "grid25", "cell_id", "sampling_unit", "locality_id"))
  has_lat <- has_any(c("latitude", "decimal_lat"))
  has_lon <- has_any(c("longitude", "decimal_lon"))
  has_genus <- has_any(c("^genus$", "genus_name"))
  has_trait <- has_any(c("balloon", "dispersal", "trait"))
  env_hits <- sum(vapply(c("precip", "temperature", "temp", "wind", "evi", "ndvi", "aridity", "pet", "topo", "elevation", "greenup", "landcover", "land_cover"), function(p) any(grepl(p, c0)), logical(1)))
  score_env <- 3 * has_key + 2 * has_lat + 2 * has_lon + env_hits + 2 * grepl("env|environment", fn)
  score_trait <- 3 * has_genus + 3 * has_trait + 2 * grepl("trait|balloon", fn)
  score_comm <- 3 * has_key + 2 * grepl("community|incidence|matrix", fn) + as.integer(length(cols) > 20)
  type <- c(environment = score_env, trait = score_trait, community = score_comm)
  list(best_type = names(which.max(type)), best_score = max(type), score_environment = score_env, score_trait = score_trait, score_community = score_comm,
       has_key = has_key, has_lat = has_lat, has_lon = has_lon, has_genus = has_genus, has_trait = has_trait, env_hits = env_hits)
}

find_files <- function(root, patterns) {
  if (!dir.exists(root)) return(character())
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files[file.exists(files) & grepl(patterns, files, ignore.case = TRUE)]
}

roots <- project_root
if (!is.na(legacy_root) && dir.exists(legacy_root)) roots <- unique(c(roots, legacy_root))

# ---------- expected current frozen inputs ----------
expected <- c(
  primary_incidence = file.path(project_root, "ANALYSIS_READY_INPUTS", "02_incidence_matrices_25km", "10_biodiversity_final_genus_by_grid25km_incidence.csv"),
  taxonomy_strict_incidence = file.path(project_root, "ANALYSIS_READY_INPUTS", "02_incidence_matrices_25km", "10_biodiversity_taxonomy_strict_genus_by_grid25km_incidence.csv"),
  cell_lookup = file.path(project_root, "ANALYSIS_READY_INPUTS", "04_spatial_reference", "10_common_grid25km_cell_lookup.csv"),
  trait_lookup = file.path(project_root, "ANALYSIS_READY_INPUTS", "03_trait_tables", "07_reviewed_genus_trait_lookup_final.csv")
)
expected_status <- data.frame(role = names(expected), path = unname(expected), exists = file.exists(expected), stringsAsFactors = FALSE)
utils::write.csv(expected_status, file.path(out_dir, "12A_expected_frozen_input_status.csv"), row.names = FALSE)

log_msg("Frozen Step 10/11 inputs found: ", sum(expected_status$exists), "/", nrow(expected_status))
for (i in seq_len(nrow(expected_status))) log_msg("  ", expected_status$role[i], ": ", ifelse(expected_status$exists[i], "FOUND", "MISSING"), " — ", expected_status$path[i])

# ---------- tabular candidate inventory ----------
table_files <- unique(unlist(lapply(roots, find_files, patterns = "\\.(csv|tsv|txt|xlsx|xls)$")))
# Avoid scanning obvious huge iteration/output tables when not useful.
keep <- !grepl("(/|^)(\\.git|baja-map-env|__pycache__)(/|$)", table_files, ignore.case = TRUE)
table_files <- table_files[keep]

records <- vector("list", length(table_files))
for (i in seq_along(table_files)) {
  p <- table_files[i]
  fi <- file.info(p)
  h <- read_header(p)
  s <- score_table(h$names, p)
  records[[i]] <- data.frame(
    path = p,
    filename = basename(p),
    size_mb = round(fi$size / 1024^2, 3),
    modified = as.character(fi$mtime),
    readable = h$ok,
    rows_read_max_2000 = h$nrow,
    n_columns = h$ncol,
    best_type = s$best_type,
    best_score = s$best_score,
    score_environment = s$score_environment,
    score_trait = s$score_trait,
    score_community = s$score_community,
    has_join_key = s$has_key,
    has_latitude = s$has_lat,
    has_longitude = s$has_lon,
    has_genus = s$has_genus,
    has_trait = s$has_trait,
    environmental_keyword_hits = s$env_hits,
    columns = paste(h$names, collapse = " | "),
    read_error = h$error,
    stringsAsFactors = FALSE
  )
}
tab_inv <- if (length(records)) do.call(rbind, records) else data.frame()
if (nrow(tab_inv)) {
  tab_inv <- tab_inv[order(-tab_inv$best_score, -tab_inv$environmental_keyword_hits, tab_inv$filename), ]
}
utils::write.csv(tab_inv, file.path(out_dir, "12A_candidate_table_inventory.csv"), row.names = FALSE)

# ---------- raster inventory ----------
raster_files <- unique(unlist(lapply(roots, find_files, patterns = "\\.(tif|tiff|img|grd|nc)$")))
raster_records <- list()
for (i in seq_along(raster_files)) {
  p <- raster_files[i]
  fi <- file.info(p)
  fn <- clean_names_simple(basename(p))
  inferred <- if (grepl("land.?cover|lc_type|mcd12", fn)) "categorical_landcover" else if (grepl("precip", fn)) "precipitation" else if (grepl("wind|u_component|v_component", fn)) "wind" else if (grepl("evi|ndvi|greenup|phenol", fn)) "vegetation_phenology" else if (grepl("arid|pet", fn)) "aridity_moisture" else if (grepl("topo|elev|dem|rugged", fn)) "topography" else if (grepl("temp", fn)) "temperature" else "unclassified"
  nlyr <- nrowr <- ncolr <- resx <- resy <- xmin <- xmax <- ymin <- ymax <- NA
  crs_txt <- NA_character_; read_error <- NA_character_
  if (requireNamespace("terra", quietly = TRUE)) {
    tryCatch({
      r <- terra::rast(p)
      nlyr <- terra::nlyr(r); nrowr <- terra::nrow(r); ncolr <- terra::ncol(r)
      rr <- terra::res(r); resx <- rr[1]; resy <- rr[2]
      ee <- terra::ext(r); xmin <- ee[1]; xmax <- ee[2]; ymin <- ee[3]; ymax <- ee[4]
      crs_txt <- terra::crs(r, proj = TRUE)
    }, error = function(e) read_error <<- conditionMessage(e))
  } else read_error <- "terra not installed"
  raster_records[[i]] <- data.frame(path = p, filename = basename(p), size_mb = round(fi$size / 1024^2, 3), modified = as.character(fi$mtime), inferred_variable = inferred,
                                    n_layers = nlyr, nrow = nrowr, ncol = ncolr, res_x = resx, res_y = resy,
                                    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, crs = crs_txt, read_error = read_error,
                                    stringsAsFactors = FALSE)
}
rast_inv <- if (length(raster_records)) do.call(rbind, raster_records) else data.frame()
utils::write.csv(rast_inv, file.path(out_dir, "12A_environmental_raster_inventory.csv"), row.names = FALSE)

# ---------- focused candidate summaries ----------
if (nrow(tab_inv)) {
  env_candidates <- subset(tab_inv, score_environment >= 5)
  trait_candidates <- subset(tab_inv, score_trait >= 5)
  comm_candidates <- subset(tab_inv, score_community >= 5)
} else {
  env_candidates <- trait_candidates <- comm_candidates <- data.frame()
}
utils::write.csv(env_candidates, file.path(out_dir, "12A_environment_table_candidates.csv"), row.names = FALSE)
utils::write.csv(trait_candidates, file.path(out_dir, "12A_trait_table_candidates.csv"), row.names = FALSE)
utils::write.csv(comm_candidates, file.path(out_dir, "12A_community_table_candidates.csv"), row.names = FALSE)

# ---------- inspect likely old community/env pair if present ----------
legacy_comm <- table_files[basename(table_files) == "community_matrix.csv"]
legacy_env <- table_files[basename(table_files) == "env_variables.csv"]
if (!is.na(legacy_root)) {
  legacy_comm <- unique(c(
    file.path(legacy_root, "Scripts Used for Ballooning", "community_matrix.csv"),
    legacy_comm
  ))
  legacy_env <- unique(c(
    file.path(legacy_root, "Scripts Used for Ballooning", "env_variables.csv"),
    legacy_env
  ))
}
legacy_comm <- legacy_comm[file.exists(legacy_comm)]
legacy_env <- legacy_env[file.exists(legacy_env)]

pair_audit <- data.frame()
if (length(legacy_comm) && length(legacy_env)) {
  cp <- legacy_comm[1]; ep <- legacy_env[1]
  tryCatch({
    comm <- utils::read.csv(cp, check.names = FALSE, stringsAsFactors = FALSE)
    env <- utils::read.csv(ep, check.names = FALSE, stringsAsFactors = FALSE)
    comm_key <- intersect(clean_names_simple(names(comm)), c("site_id", "grid25km_id", "cell_id"))
    env_key <- intersect(clean_names_simple(names(env)), c("site_id", "grid25km_id", "cell_id"))
    pair_audit <- data.frame(
      community_path = cp,
      environment_path = ep,
      community_rows = nrow(comm),
      environment_rows = nrow(env),
      equal_row_counts = nrow(comm) == nrow(env),
      community_has_site_id = "site_id" %in% names(comm),
      environment_has_site_id = "site_id" %in% names(env),
      duplicate_community_site_id = if ("site_id" %in% names(comm)) sum(duplicated(comm$site_id)) else NA,
      duplicate_environment_site_id = if ("site_id" %in% names(env)) sum(duplicated(env$site_id)) else NA,
      identical_site_id_order = if (all(c("site_id") %in% names(comm)) && all(c("site_id") %in% names(env)) && nrow(comm) == nrow(env)) identical(as.character(comm$site_id), as.character(env$site_id)) else NA,
      stringsAsFactors = FALSE
    )
  }, error = function(e) log_msg("Legacy pair audit failed: ", conditionMessage(e)))
}
utils::write.csv(pair_audit, file.path(out_dir, "12A_legacy_community_environment_pair_audit.csv"), row.names = FALSE)

# ---------- recommendations ----------
rec_file <- file.path(out_dir, "12A_recommended_next_inputs.txt")
rec <- c(
  "STEP 12A ENVIRONMENTAL GLM INPUT AUDIT — NEXT ACTIONS",
  "===================================================",
  "",
  "Required model unit:",
  "  One row per occupied 25-km cell, with counts of ballooning and non-ballooning genera",
  "  derived from the frozen genus × grid-cell incidence matrix and reviewed genus trait lookup.",
  "",
  "Required environmental join:",
  "  Environmental values must join explicitly by the same 25-km cell identifier or by cell geometry.",
  "  Do not bind rows by position.",
  "",
  "Preferred predictors to evaluate (subject to availability and collinearity):",
  "  latitude + latitude^2, precipitation, mean/max temperature, temperature variability,",
  "  wind speed, vegetation phenology/productivity, aridity or PET, topographic diversity,",
  "  and a defensible land-cover representation.",
  "",
  "Model diagnostics required before mapping:",
  "  missingness, predictor distributions, correlation/VIF, overdispersion, influential cells,",
  "  residual spatial autocorrelation, and sensitivity to taxonomy/trait confidence.",
  "",
  "Prediction rule:",
  "  A manuscript prediction map must be generated from aligned environmental raster predictors",
  "  using the validated model. Rasterizing fitted values at observed cells is not a prediction surface.",
  "",
  paste0("Candidate environmental tables found: ", nrow(env_candidates)),
  paste0("Candidate community/incidence tables found: ", nrow(comm_candidates)),
  paste0("Candidate trait tables found: ", nrow(trait_candidates)),
  paste0("Environmental rasters found: ", nrow(rast_inv)),
  "",
  "Open these files next:",
  "  12A_environment_table_candidates.csv",
  "  12A_environmental_raster_inventory.csv",
  "  12A_legacy_community_environment_pair_audit.csv",
  "  12A_expected_frozen_input_status.csv"
)
writeLines(rec, rec_file)

# ---------- validation ----------
validation <- data.frame(
  check = c(
    "primary_incidence_exists",
    "cell_lookup_exists",
    "trait_lookup_exists",
    "environment_table_candidate_found",
    "environmental_raster_found"
  ),
  passed = c(
    file.exists(expected[["primary_incidence"]]),
    file.exists(expected[["cell_lookup"]]),
    file.exists(expected[["trait_lookup"]]),
    nrow(env_candidates) > 0,
    nrow(rast_inv) > 0
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(validation, file.path(out_dir, "12A_validation.csv"), row.names = FALSE)

log_msg("Candidate tables inventoried: ", nrow(tab_inv))
log_msg("  Environmental candidates: ", nrow(env_candidates))
log_msg("  Community/incidence candidates: ", nrow(comm_candidates))
log_msg("  Trait candidates: ", nrow(trait_candidates))
log_msg("Environmental rasters inventoried: ", nrow(rast_inv))
log_msg("Validation checks passed: ", sum(validation$passed), "/", nrow(validation))
log_msg("")
log_msg("STEP 12A COMPLETED SUCCESSFULLY")
log_msg("Outputs: ", out_dir)
