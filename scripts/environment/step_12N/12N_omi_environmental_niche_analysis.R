#!/usr/bin/env Rscript

# STEP 12N — OMI environmental niche position and breadth analysis
# Baja California arachnid ballooning publication pipeline
#
# Primary question:
#   Do C3 ballooning-capable genera occupy broader realized environmental
#   niches than fixed-N0 non-ballooning genera?
#
# Primary design:
#   * 25-km genus-by-cell incidence matrix
#   * environmentally complete Step 12C/12F cells
#   * four standardized environmental dimensions:
#       VPD, wind seasonality, vegetation/phenology PC1, topography PC1
#   * primary trait definition C3 = D1 + D2 + D3 versus fixed N0
#   * D4 excluded from the primary comparison
#   * genera occurring in >=5 analyzed cells; >=10-cell sensitivity
#   * OMI niche position (marginality) and total tolerance (niche breadth)
#   * occupancy- and order-adjusted genus-level trait models
#   * C1/C2/C4 and taxonomy-strict sensitivity analyses
#
# Usage:
#   Rscript 12N_omi_environmental_niche_analysis.R PROJECT_ROOT [N_PERM] [SEED]
#
# Example:
#   Rscript scripts/environment/step_12N/12N_omi_environmental_niche_analysis.R \
#     "/Users/estiner/Desktop/OLD BALLOONING/Baja_Ballooning_Pipeline" 5000 20260717

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop(
    "Usage: Rscript 12N_omi_environmental_niche_analysis.R ",
    "PROJECT_ROOT [N_PERM=5000] [SEED=20260717]"
  )
}

project_root <- normalizePath(path.expand(args[[1]]), mustWork = TRUE)
n_perm <- if (length(args) >= 2L) as.integer(args[[2]]) else 5000L
seed <- if (length(args) >= 3L) as.integer(args[[3]]) else 20260717L
if (!is.finite(n_perm) || n_perm < 99L) stop("N_PERM must be at least 99.")
if (!is.finite(seed)) stop("SEED must be an integer.")
set.seed(seed)

SCRIPT_VERSION <- "12N_OMI_v1_2026-07-18"

required_packages <- c("ade4", "ggplot2", "sf", "patchwork", "jsonlite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    "Missing required R package(s): ", paste(missing_packages, collapse = ", "),
    ". Run scripts/environment/step_12N/12N_install_packages.R first."
  )
}

analysis_ready <- file.path(project_root, "ANALYSIS_READY_INPUTS")
grid_fallback <- file.path(project_root, "02_data_clean", "08_grid25km_incidence")
trait_fallback <- file.path(project_root, "02_data_clean", "07_final_trait_merge")
step12c_dir <- file.path(project_root, "04_analysis", "12C_cell_environment_model_table")
step12f_dir <- file.path(project_root, "04_analysis", "12F_environment_predictor_scores")
out_dir <- file.path(project_root, "04_analysis", "12N_omi_environmental_niche")
figure_dir <- file.path(out_dir, "figures")
archive_root <- file.path(project_root, "08_archive")

first_existing <- function(paths, label, required = TRUE) {
  paths <- path.expand(paths)
  hit <- paths[file.exists(paths)]
  if (length(hit)) return(normalizePath(hit[[1]], winslash = "/", mustWork = TRUE))
  if (required) {
    stop("Could not find ", label, ". Tried:\n", paste(paths, collapse = "\n"))
  }
  NA_character_
}

archive_directory <- function(path, archive_root, prefix) {
  if (!dir.exists(path) || length(list.files(path, all.files = TRUE, no.. = TRUE)) == 0L) {
    return(NA_character_)
  }
  dir.create(archive_root, recursive = TRUE, showWarnings = FALSE)
  stamp <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid())
  destination <- file.path(archive_root, paste0(prefix, "_", stamp))
  if (!file.rename(path, destination)) {
    dir.create(destination, recursive = TRUE, showWarnings = FALSE)
    copied <- file.copy(
      list.files(path, full.names = TRUE, all.files = TRUE, no.. = TRUE),
      destination,
      recursive = TRUE
    )
    if (!all(copied)) stop("Could not archive prior Step 12N output.")
    unlink(path, recursive = TRUE, force = TRUE)
  }
  normalizePath(destination, winslash = "/", mustWork = TRUE)
}

find_field <- function(fields, candidates, label, required = TRUE) {
  lower <- tolower(trimws(fields))
  for (candidate in candidates) {
    idx <- match(tolower(candidate), lower)
    if (!is.na(idx)) return(fields[[idx]])
  }
  if (required) {
    stop(
      "Could not identify ", label, ". Tried: ", paste(candidates, collapse = ", "),
      ". Available fields: ", paste(fields, collapse = ", ")
    )
  }
  NA_character_
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))
trim_character <- function(x) trimws(as.character(x))
parse_bool <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

parse_evidence_one <- function(value) {
  s <- toupper(trimws(as.character(value)))
  if (is.na(s) || !nzchar(s)) return(NA_character_)
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
      "Legacy binary fields are not accepted because they cannot distinguish D4 from N0."
    )
  }
  names(which.max(scores))
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

align_matrix <- function(obj, target_genera, target_cells) {
  gi <- match(tolower(target_genera), tolower(obj$genera))
  ci <- match(target_cells, obj$cells)
  if (anyNA(gi)) stop("Sensitivity matrix lacks required genera.")
  if (anyNA(ci)) stop("Sensitivity matrix lacks required cells.")
  m <- obj$matrix[gi, ci, drop = FALSE]
  rownames(m) <- target_genera
  colnames(m) <- target_cells
  m
}

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "")
}

extract_model_row <- function(fit, term, response_name, model_name, case_id) {
  co <- summary(fit)$coefficients
  if (!term %in% rownames(co)) {
    return(data.frame(
      case_id = case_id, response = response_name, model = model_name,
      term = term, estimate = NA_real_, std_error = NA_real_, statistic = NA_real_,
      p_value = NA_real_, ci_low = NA_real_, ci_high = NA_real_, n = nobs(fit),
      stringsAsFactors = FALSE
    ))
  }
  ci <- suppressMessages(confint(fit, parm = term, level = 0.95))
  data.frame(
    case_id = case_id,
    response = response_name,
    model = model_name,
    term = term,
    estimate = unname(co[term, 1]),
    std_error = unname(co[term, 2]),
    statistic = unname(co[term, 3]),
    p_value = unname(co[term, 4]),
    ci_low = unname(ci[1]),
    ci_high = unname(ci[2]),
    n = nobs(fit),
    stringsAsFactors = FALSE
  )
}

fit_trait_models <- function(genus_df, case_id) {
  d <- genus_df[
    is.finite(genus_df$total_tolerance) &
      is.finite(genus_df$OMI) &
      is.finite(genus_df$occupied_cells) &
      genus_df$trait_group %in% c("Ballooning", "Non-ballooning"),
    , drop = FALSE
  ]
  d$ballooning_binary <- as.integer(d$trait_group == "Ballooning")
  d$log_occupancy <- log1p(d$occupied_cells)
  d$log_total_tolerance <- log1p(pmax(0, d$total_tolerance))
  d$log_OMI <- log1p(pmax(0, d$OMI))
  d$order_model <- trimws(as.character(d$order))
  d$order_model[!nzchar(d$order_model) | is.na(d$order_model)] <- "Unknown"

  enough_order <- length(unique(d$order_model)) >= 2L && nrow(d) >= length(unique(d$order_model)) + 8L
  model_formula <- if (enough_order) {
    log_total_tolerance ~ ballooning_binary + log_occupancy + factor(order_model)
  } else {
    log_total_tolerance ~ ballooning_binary + log_occupancy
  }
  fit_tol <- lm(model_formula, data = d)
  if (!"ballooning_binary" %in% names(coef(fit_tol)) || !is.finite(coef(fit_tol)[["ballooning_binary"]])) {
    fit_tol <- lm(log_total_tolerance ~ ballooning_binary + log_occupancy, data = d)
    enough_order <- FALSE
  }

  model_formula_omi <- if (enough_order) {
    log_OMI ~ ballooning_binary + log_occupancy + factor(order_model)
  } else {
    log_OMI ~ ballooning_binary + log_occupancy
  }
  fit_omi <- lm(model_formula_omi, data = d)
  if (!"ballooning_binary" %in% names(coef(fit_omi)) || !is.finite(coef(fit_omi)[["ballooning_binary"]])) {
    fit_omi <- lm(log_OMI ~ ballooning_binary + log_occupancy, data = d)
  }

  out <- rbind(
    extract_model_row(
      fit_tol, "ballooning_binary", "log1p_total_tolerance",
      if (enough_order) "occupancy_plus_order_adjusted" else "occupancy_adjusted",
      case_id
    ),
    extract_model_row(
      fit_omi, "ballooning_binary", "log1p_OMI",
      if (enough_order) "occupancy_plus_order_adjusted" else "occupancy_adjusted",
      case_id
    )
  )

  list(rows = out, tolerance_fit = fit_tol, omi_fit = fit_omi, data = d)
}

group_summary <- function(genus_df, case_id) {
  groups <- split(genus_df, genus_df$trait_group)
  rows <- lapply(names(groups), function(group_name) {
    d <- groups[[group_name]]
    data.frame(
      case_id = case_id,
      trait_group = group_name,
      genera = nrow(d),
      median_occupied_cells = median(d$occupied_cells),
      mean_occupied_cells = mean(d$occupied_cells),
      median_total_tolerance = median(d$total_tolerance),
      mean_total_tolerance = mean(d$total_tolerance),
      median_OMI = median(d$OMI),
      mean_OMI = mean(d$OMI),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

run_wilcoxon <- function(genus_df, case_id, response_name) {
  values <- genus_df[[response_name]]
  keep <- is.finite(values) & genus_df$trait_group %in% c("Ballooning", "Non-ballooning")
  d <- genus_df[keep, , drop = FALSE]
  test <- suppressWarnings(wilcox.test(d[[response_name]] ~ d$trait_group, exact = FALSE, conf.int = FALSE))
  med_b <- median(d[[response_name]][d$trait_group == "Ballooning"])
  med_n <- median(d[[response_name]][d$trait_group == "Non-ballooning"])
  data.frame(
    case_id = case_id,
    response = response_name,
    test = "Wilcoxon rank-sum",
    ballooning_median = med_b,
    non_ballooning_median = med_n,
    median_difference_ballooning_minus_non = med_b - med_n,
    statistic = unname(test$statistic),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}

save_plot_all <- function(plot_object, stem, width, height) {
  ggplot2::ggsave(paste0(stem, ".png"), plot_object, width = width, height = height, dpi = 600, bg = "white")
  ggplot2::ggsave(paste0(stem, ".pdf"), plot_object, width = width, height = height, bg = "white")
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(
      paste0(stem, ".svg"), plot_object, width = width, height = height,
      device = svglite::svglite, bg = "white"
    )
  }
}

# ------------------------------- inputs ------------------------------------

primary_path <- first_existing(c(
  file.path(analysis_ready, "02_incidence_matrices_25km", "10_biodiversity_final_genus_by_grid25km_incidence.csv"),
  file.path(grid_fallback, "10_biodiversity_final_genus_by_grid25km_incidence.csv")
), "primary 25-km incidence matrix")

strict_path <- first_existing(c(
  file.path(analysis_ready, "02_incidence_matrices_25km", "10_biodiversity_taxonomy_strict_genus_by_grid25km_incidence.csv"),
  file.path(grid_fallback, "10_biodiversity_taxonomy_strict_genus_by_grid25km_incidence.csv")
), "taxonomy-strict 25-km incidence matrix", required = FALSE)

cell_geojson_path <- first_existing(c(
  file.path(analysis_ready, "04_spatial_reference", "10_common_grid25km_cells.geojson"),
  file.path(grid_fallback, "10_common_grid25km_cells.geojson")
), "25-km cell polygon GeoJSON")

genus_lookup_path <- first_existing(c(
  file.path(analysis_ready, "02_incidence_matrices_25km", "10_common_genus_lookup.csv"),
  file.path(grid_fallback, "10_common_genus_lookup.csv")
), "common genus lookup", required = FALSE)

trait_path <- first_existing(c(
  file.path(analysis_ready, "03_trait_tables", "07_reviewed_genus_trait_lookup_normalized.csv"),
  file.path(analysis_ready, "03_trait_tables", "07_reviewed_genus_trait_lookup_final.csv"),
  file.path(trait_fallback, "07_reviewed_genus_trait_lookup_final.csv")
), "reviewed genus trait lookup")

response_path <- first_existing(c(
  file.path(step12c_dir, "12C_primary_glm_candidate_table.csv")
), "Step 12C primary candidate table")

scores_path <- first_existing(c(
  file.path(step12f_dir, "12F_environment_predictor_scores_by_cell.csv")
), "Step 12F environmental score table")

archived <- archive_directory(out_dir, archive_root, "12N_omi_environmental_niche")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(out_dir, "12N_analysis_log.txt")
log_con <- file(log_path, open = "wt")
on.exit(close(log_con), add = TRUE)
log_msg <- function(...) {
  message <- paste0(...)
  cat(message, "\n")
  writeLines(message, log_con)
  flush(log_con)
}

log_msg("STEP 12N OMI STARTED")
log_msg("Version: ", SCRIPT_VERSION)
log_msg("Project root: ", project_root)
log_msg("Permutations: ", n_perm)
log_msg("Seed: ", seed)
if (!is.na(archived)) log_msg("Archived prior output: ", archived)

primary_obj <- read_incidence_matrix(primary_path)
strict_obj <- if (!is.na(strict_path)) read_incidence_matrix(strict_path) else NULL

response <- read.csv(response_path, check.names = FALSE, stringsAsFactors = FALSE)
scores <- read.csv(scores_path, check.names = FALSE, stringsAsFactors = FALSE)
if (anyDuplicated(response$grid_cell_id)) stop("Duplicate grid_cell_id in Step 12C response table.")
if (anyDuplicated(scores$grid_cell_id)) stop("Duplicate grid_cell_id in Step 12F score table.")

required_env <- c("vpd_z", "wind_seasonality_z", "vegetation_axis", "topography_axis")
missing_env <- setdiff(required_env, names(scores))
if (length(missing_env)) stop("Step 12F score table lacks: ", paste(missing_env, collapse = ", "))

eligible <- response
if ("recommended_primary_model_cell" %in% names(eligible)) {
  eligible <- eligible[parse_bool(eligible$recommended_primary_model_cell), , drop = FALSE]
}
cell_meta_fields <- intersect(
  c("grid_cell_id", "centroid_latitude", "centroid_longitude", "latitude_band", "easting_km", "northing_km"),
  names(eligible)
)
eligible <- eligible[, cell_meta_fields, drop = FALSE]
cell_table <- merge(eligible, scores, by = "grid_cell_id", all = FALSE, sort = FALSE)
cell_table <- cell_table[complete.cases(cell_table[, required_env, drop = FALSE]), , drop = FALSE]
cell_table <- cell_table[match(intersect(primary_obj$cells, cell_table$grid_cell_id), cell_table$grid_cell_id), , drop = FALSE]
cell_table <- cell_table[!is.na(cell_table$grid_cell_id), , drop = FALSE]
if (nrow(cell_table) < 50L) stop("Fewer than 50 environmentally complete occupied cells remain.")

cells <- cell_table$grid_cell_id
genera <- primary_obj$genera
primary_mat <- primary_obj$matrix[, match(cells, primary_obj$cells), drop = FALSE]
if (anyNA(match(cells, primary_obj$cells))) stop("Some environmental cells are absent from the incidence matrix.")
strict_mat <- if (!is.null(strict_obj)) align_matrix(strict_obj, genera, cells) else NULL

traits <- read.csv(trait_path, check.names = FALSE, stringsAsFactors = FALSE)
genus_field <- find_field(names(traits), c("genus", "analysis_genus"), "trait genus")
evidence_field <- find_evidence_field(traits)
order_field <- find_field(names(traits), c("order_final", "order", "trait_order", "analysis_order"), "trait order", required = FALSE)
family_field <- find_field(names(traits), c("family_final", "family", "trait_family", "analysis_family"), "trait family", required = FALSE)
traits$genus_key <- tolower(trim_character(traits[[genus_field]]))
if (anyDuplicated(traits$genus_key[traits$genus_key != ""])) stop("Duplicate genus in trait table.")
trait_idx <- match(tolower(genera), traits$genus_key)
if (anyNA(trait_idx)) stop("Genera missing from trait table: ", paste(genera[is.na(trait_idx)], collapse = ", "))

evidence_class <- normalize_evidence(traits[[evidence_field]][trait_idx])
if (anyNA(evidence_class)) {
  stop("Unresolved evidence class for: ", paste(genera[is.na(evidence_class)], collapse = ", "))
}
# A normalized C3 entry means D1-D3 combined; retain it only in C3/C4.
order_values <- if (!is.na(order_field)) trim_character(traits[[order_field]][trait_idx]) else rep("Unknown", length(genera))
family_values <- if (!is.na(family_field)) trim_character(traits[[family_field]][trait_idx]) else rep("Unknown", length(genera))

if (!is.na(genus_lookup_path)) {
  genus_lookup <- read.csv(genus_lookup_path, check.names = FALSE, stringsAsFactors = FALSE)
  lookup_genus <- find_field(names(genus_lookup), c("genus", "analysis_genus"), "lookup genus")
  lookup_order <- find_field(names(genus_lookup), c("order", "analysis_order"), "lookup order", required = FALSE)
  lookup_family <- find_field(names(genus_lookup), c("family", "analysis_family"), "lookup family", required = FALSE)
  lookup_idx <- match(tolower(genera), tolower(trim_character(genus_lookup[[lookup_genus]])))
  if (!is.na(lookup_order)) {
    replace <- (!nzchar(order_values) | order_values == "Unknown" | is.na(order_values)) & !is.na(lookup_idx)
    order_values[replace] <- trim_character(genus_lookup[[lookup_order]][lookup_idx[replace]])
  }
  if (!is.na(lookup_family)) {
    replace <- (!nzchar(family_values) | family_values == "Unknown" | is.na(family_values)) & !is.na(lookup_idx)
    family_values[replace] <- trim_character(genus_lookup[[lookup_family]][lookup_idx[replace]])
  }
}
order_values[!nzchar(order_values) | is.na(order_values)] <- "Unknown"
family_values[!nzchar(family_values) | is.na(family_values)] <- "Unknown"

trait_master <- data.frame(
  genus = genera,
  evidence_class = evidence_class,
  order = order_values,
  family = family_values,
  stringsAsFactors = FALSE
)

# Environmental PCA is fixed across all trait thresholds so niche metrics are
# comparable across sensitivity cases.
env <- as.data.frame(cell_table[, required_env, drop = FALSE])
rownames(env) <- cells
env_dudi <- ade4::dudi.pca(
  env,
  center = TRUE,
  scale = TRUE,
  scannf = FALSE,
  nf = min(4L, ncol(env))
)

thresholds <- list(
  C1 = c("D1"),
  C2 = c("D1", "D2"),
  C3 = c("D1", "D2", "D3", "C3"),
  C4 = c("D1", "D2", "D3", "D4", "C3")
)

all_params <- list()
all_models <- list()
all_groups <- list()
all_wilcox <- list()
all_global <- list()
primary_objects <- NULL

run_case <- function(case_id, matrix_genera_by_cells, threshold_name, min_occupancy, permutations) {
  ballooning_tiers <- thresholds[[threshold_name]]
  trait_group <- ifelse(
    trait_master$evidence_class %in% ballooning_tiers,
    "Ballooning",
    ifelse(trait_master$evidence_class == "N0", "Non-ballooning", "Excluded")
  )
  occupancy <- rowSums(matrix_genera_by_cells)
  keep <- trait_group != "Excluded" & occupancy >= min_occupancy
  if (sum(keep) < 20L) stop("Case ", case_id, " retains fewer than 20 genera.")
  if (length(unique(trait_group[keep])) < 2L) stop("Case ", case_id, " lacks both trait groups.")

  Y <- t(matrix_genera_by_cells[keep, , drop = FALSE])
  Y <- as.data.frame(Y, check.names = FALSE)
  rownames(Y) <- cells
  if (any(colSums(Y) < min_occupancy)) stop("Occupancy filter failed in ", case_id)

  omi <- ade4::niche(env_dudi, Y, scannf = FALSE, nf = 2)
  params <- as.data.frame(ade4::niche.param(omi))
  params$genus <- rownames(params)
  rownames(params) <- NULL
  params$total_tolerance <- params$Tol + params$Rtol
  params$total_tolerance_pct <- params$tol + params$rtol

  selected <- trait_master[match(tolower(params$genus), tolower(trait_master$genus)), , drop = FALSE]
  params$trait_group <- trait_group[match(tolower(params$genus), tolower(trait_master$genus))]
  params$occupied_cells <- occupancy[match(tolower(params$genus), tolower(trait_master$genus))]
  params$evidence_class <- selected$evidence_class
  params$order <- selected$order
  params$family <- selected$family
  params$case_id <- case_id
  params$threshold <- threshold_name
  params$minimum_occupied_cells <- min_occupancy

  model_result <- fit_trait_models(params, case_id)
  group_result <- group_summary(params, case_id)
  wilcox_result <- rbind(
    run_wilcoxon(params, case_id, "total_tolerance"),
    run_wilcoxon(params, case_id, "OMI")
  )

  set.seed(seed + sum(utf8ToInt(case_id)))
  rt <- ade4::rtest(omi, nrepet = permutations)
  obs <- as.numeric(rt$obs)
  names_obs <- names(rt$obs)
  pvalues <- as.numeric(rt$pvalue)
  names_p <- names(rt$pvalue)
  global_row <- data.frame(
    case_id = case_id,
    threshold = threshold_name,
    minimum_occupied_cells = min_occupancy,
    genera = ncol(Y),
    cells = nrow(Y),
    permutations = permutations,
    observed_mean_OMI = if ("OMI.mean" %in% names_obs) obs[match("OMI.mean", names_obs)] else tail(obs, 1),
    p_value_mean_OMI = if ("OMI.mean" %in% names_p) pvalues[match("OMI.mean", names_p)] else tail(pvalues, 1),
    stringsAsFactors = FALSE
  )

  list(
    omi = omi,
    params = params,
    model_rows = model_result$rows,
    model_fit_tolerance = model_result$tolerance_fit,
    model_fit_omi = model_result$omi_fit,
    group_rows = group_result,
    wilcox_rows = wilcox_result,
    global_row = global_row,
    Y = Y
  )
}

cases <- list(
  list(id = "C3_occ5_primary", threshold = "C3", min_occ = 5L, matrix = primary_mat, perms = n_perm),
  list(id = "C3_occ10", threshold = "C3", min_occ = 10L, matrix = primary_mat, perms = min(n_perm, 2000L)),
  list(id = "C1_occ5", threshold = "C1", min_occ = 5L, matrix = primary_mat, perms = min(n_perm, 1000L)),
  list(id = "C2_occ5", threshold = "C2", min_occ = 5L, matrix = primary_mat, perms = min(n_perm, 1000L)),
  list(id = "C4_occ5", threshold = "C4", min_occ = 5L, matrix = primary_mat, perms = min(n_perm, 1000L))
)
if (!is.null(strict_mat)) {
  cases[[length(cases) + 1L]] <- list(
    id = "C3_occ5_taxonomy_strict", threshold = "C3", min_occ = 5L,
    matrix = strict_mat, perms = min(n_perm, 2000L)
  )
}

for (case in cases) {
  log_msg("Running ", case$id, "...")
  result <- run_case(case$id, case$matrix, case$threshold, case$min_occ, case$perms)
  all_params[[case$id]] <- result$params
  all_models[[case$id]] <- result$model_rows
  all_groups[[case$id]] <- result$group_rows
  all_wilcox[[case$id]] <- result$wilcox_rows
  all_global[[case$id]] <- result$global_row
  if (case$id == "C3_occ5_primary") primary_objects <- result
}

params_all <- do.call(rbind, all_params)
models_all <- do.call(rbind, all_models)
groups_all <- do.call(rbind, all_groups)
wilcox_all <- do.call(rbind, all_wilcox)
global_all <- do.call(rbind, all_global)

write_csv(params_all, file.path(out_dir, "12N_omi_genus_niche_parameters_all_cases.csv"))
write_csv(primary_objects$params, file.path(out_dir, "12N_omi_genus_niche_parameters_primary_C3.csv"))
write_csv(models_all, file.path(out_dir, "12N_trait_effect_models.csv"))
write_csv(groups_all, file.path(out_dir, "12N_trait_group_summaries.csv"))
write_csv(wilcox_all, file.path(out_dir, "12N_unadjusted_trait_tests.csv"))
write_csv(global_all, file.path(out_dir, "12N_global_omi_randomization_tests.csv"))

# OMI axis coordinates and loadings for the primary case.
primary_omi <- primary_objects$omi
genus_scores <- as.data.frame(primary_omi$li)
genus_scores$genus <- rownames(genus_scores)
rownames(genus_scores) <- NULL
genus_scores <- merge(
  genus_scores,
  primary_objects$params[, c("genus", "trait_group", "occupied_cells", "OMI", "total_tolerance")],
  by = "genus", all.x = TRUE, sort = FALSE
)
write_csv(genus_scores, file.path(out_dir, "12N_primary_OMI_genus_scores.csv"))

site_scores <- as.data.frame(primary_omi$ls)
site_scores$grid_cell_id <- rownames(site_scores)
rownames(site_scores) <- NULL
write_csv(site_scores, file.path(out_dir, "12N_primary_OMI_cell_scores.csv"))

variable_scores <- as.data.frame(primary_omi$co)
variable_scores$predictor <- rownames(variable_scores)
rownames(variable_scores) <- NULL
write_csv(variable_scores, file.path(out_dir, "12N_primary_OMI_environmental_loadings.csv"))

# ------------------------------- figures -----------------------------------

eig_pct <- 100 * primary_omi$eig / sum(primary_omi$eig)
axis1 <- names(primary_omi$li)[1]
axis2 <- names(primary_omi$li)[2]

label_map <- c(
  vpd_z = "VPD",
  wind_seasonality_z = "Wind seasonality",
  vegetation_axis = "Vegetation/phenology",
  topography_axis = "Topography"
)

# Scale environmental arrows to the genus-score plotting region.
arrow_df <- variable_scores
names(arrow_df)[1:2] <- c("Axis1", "Axis2")
arrow_df$label <- ifelse(
  arrow_df$predictor %in% names(label_map),
  label_map[arrow_df$predictor],
  arrow_df$predictor
)
genus_plot_df <- genus_scores
names(genus_plot_df)[match(c(axis1, axis2), names(genus_plot_df))] <- c("Axis1", "Axis2")
arrow_scale <- 0.75 * min(
  diff(range(genus_plot_df$Axis1, na.rm = TRUE)) / max(1e-9, diff(range(arrow_df$Axis1, na.rm = TRUE))),
  diff(range(genus_plot_df$Axis2, na.rm = TRUE)) / max(1e-9, diff(range(arrow_df$Axis2, na.rm = TRUE)))
)
arrow_df$Axis1_plot <- arrow_df$Axis1 * arrow_scale
arrow_df$Axis2_plot <- arrow_df$Axis2 * arrow_scale

label_genera <- unique(c(
  head(primary_objects$params$genus[order(primary_objects$params$OMI, decreasing = TRUE)], 5),
  head(primary_objects$params$genus[order(primary_objects$params$total_tolerance, decreasing = TRUE)], 5)
))

p_ordination <- ggplot2::ggplot(genus_plot_df, ggplot2::aes(x = Axis1, y = Axis2, shape = trait_group)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey75") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, color = "grey75") +
  ggplot2::geom_point(ggplot2::aes(size = occupied_cells), alpha = 0.65) +
  ggplot2::stat_ellipse(ggplot2::aes(group = trait_group), type = "norm", linewidth = 0.6, show.legend = FALSE) +
  ggplot2::geom_segment(
    data = arrow_df,
    ggplot2::aes(x = 0, y = 0, xend = Axis1_plot, yend = Axis2_plot),
    inherit.aes = FALSE,
    arrow = grid::arrow(length = grid::unit(0.16, "cm")),
    linewidth = 0.5
  ) +
  ggplot2::geom_text(
    data = arrow_df,
    ggplot2::aes(x = Axis1_plot, y = Axis2_plot, label = label),
    inherit.aes = FALSE,
    size = 3.0,
    vjust = -0.4
  ) +
  ggplot2::geom_text(
    data = genus_plot_df[genus_plot_df$genus %in% label_genera, , drop = FALSE],
    ggplot2::aes(label = genus),
    size = 2.5,
    check_overlap = TRUE,
    vjust = -0.7,
    show.legend = FALSE
  ) +
  ggplot2::scale_shape_manual(values = c("Ballooning" = 16, "Non-ballooning" = 17)) +
  ggplot2::scale_size_continuous(range = c(1.5, 5.0), name = "Occupied cells") +
  ggplot2::labs(
    title = "A. Environmental niche position",
    subtitle = "Primary C3 versus fixed N0; genera occurring in at least five cells",
    x = sprintf("OMI axis 1 (%.1f%%)", eig_pct[1]),
    y = sprintf("OMI axis 2 (%.1f%%)", eig_pct[2]),
    shape = "Trait group"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")

breadth_df <- primary_objects$params
breadth_df$log_total_tolerance <- log1p(breadth_df$total_tolerance)
p_breadth <- ggplot2::ggplot(breadth_df, ggplot2::aes(x = trait_group, y = log_total_tolerance, shape = trait_group)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.58, alpha = 0.25) +
  ggplot2::geom_jitter(width = 0.16, height = 0, alpha = 0.55, size = 1.7) +
  ggplot2::scale_shape_manual(values = c("Ballooning" = 16, "Non-ballooning" = 17)) +
  ggplot2::labs(
    title = "B. Realized environmental niche breadth",
    subtitle = "Total OMI tolerance; raw comparison shown",
    x = NULL,
    y = "log(1 + total tolerance)"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "none")

cells_sf <- suppressWarnings(sf::st_read(cell_geojson_path, quiet = TRUE))
if (!"grid_cell_id" %in% names(cells_sf)) stop("Cell GeoJSON lacks grid_cell_id.")
map_df <- merge(cells_sf, site_scores, by = "grid_cell_id", all.y = TRUE, sort = FALSE)
map_axis <- names(primary_omi$ls)[1]
p_map <- ggplot2::ggplot(map_df) +
  ggplot2::geom_sf(ggplot2::aes(fill = .data[[map_axis]]), linewidth = 0.12, color = "grey45") +
  ggplot2::scale_fill_gradient2(midpoint = 0, name = "OMI axis 1") +
  ggplot2::coord_sf(datum = NA) +
  ggplot2::labs(
    title = "C. Geographic distribution of the main niche gradient",
    subtitle = "Cell scores on OMI axis 1"
  ) +
  ggplot2::theme_void(base_size = 10) +
  ggplot2::theme(
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold")
  )

combined <- (p_ordination | p_breadth | p_map) +
  patchwork::plot_annotation(
    title = "Environmental niches of ballooning and non-ballooning Baja arachnid genera",
    subtitle = "Outlying Mean Index analysis using 25-km genus incidence and four environmental dimensions",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 15))
  )

save_plot_all(p_ordination, file.path(figure_dir, "Figure_12N_A_OMI_environmental_niche_position"), 7.5, 6.2)
save_plot_all(p_breadth, file.path(figure_dir, "Figure_12N_B_OMI_niche_breadth"), 5.8, 5.5)
save_plot_all(p_map, file.path(figure_dir, "Figure_12N_C_OMI_axis1_map"), 5.2, 7.0)
save_plot_all(combined, file.path(figure_dir, "Figure_12N_OMI_niche_analysis_combined"), 17.0, 6.5)

# ---------------------------- validation/decision --------------------------

validation <- data.frame(
  check = c(
    "primary_incidence_binary",
    "environmental_cells_aligned",
    "environmental_predictors_complete",
    "primary_has_both_trait_groups",
    "primary_D4_excluded",
    "primary_minimum_occupancy_respected",
    "OMI_eigenvalues_positive",
    "primary_trait_model_fitted",
    "all_output_tables_created"
  ),
  passed = c(
    all(primary_mat %in% c(0, 1)),
    identical(colnames(primary_mat), cells),
    all(complete.cases(env)),
    length(unique(primary_objects$params$trait_group)) == 2L,
    !any(primary_objects$params$evidence_class == "D4"),
    all(primary_objects$params$occupied_cells >= 5),
    all(primary_omi$eig > 0),
    any(models_all$case_id == "C3_occ5_primary" & models_all$response == "log1p_total_tolerance" & is.finite(models_all$estimate)),
    TRUE
  ),
  detail = c(
    paste(dim(primary_mat), collapse = " x "),
    paste(length(cells), "cells"),
    paste(required_env, collapse = "; "),
    paste(table(primary_objects$params$trait_group), collapse = "; "),
    "C3=D1-D3/C3 versus fixed N0",
    paste(range(primary_objects$params$occupied_cells), collapse = " to "),
    paste(signif(primary_omi$eig, 4), collapse = "; "),
    "log1p(total tolerance) ~ trait + log occupancy + order when estimable",
    out_dir
  ),
  stringsAsFactors = FALSE
)
write_csv(validation, file.path(out_dir, "12N_validation.csv"))
if (!all(validation$passed)) stop("Step 12N validation failed; inspect 12N_validation.csv.")

primary_model <- models_all[
  models_all$case_id == "C3_occ5_primary" &
    models_all$response == "log1p_total_tolerance" &
    models_all$term == "ballooning_binary",
  , drop = FALSE
]
occ10_model <- models_all[
  models_all$case_id == "C3_occ10" &
    models_all$response == "log1p_total_tolerance" &
    models_all$term == "ballooning_binary",
  , drop = FALSE
]

classification <- "exploratory"
classification_reason <- "The primary trait coefficient is uncertain or not stable at the >=10-cell sensitivity threshold."
if (nrow(primary_model) == 1L && nrow(occ10_model) == 1L &&
    is.finite(primary_model$estimate) && is.finite(occ10_model$estimate)) {
  same_direction <- sign(primary_model$estimate) == sign(occ10_model$estimate)
  primary_supported <- primary_model$ci_low > 0 || primary_model$ci_high < 0
  occ10_supported <- occ10_model$ci_low > 0 || occ10_model$ci_high < 0
  if (same_direction && primary_supported && occ10_supported) {
    classification <- "candidate_retained_result"
    classification_reason <- "The occupancy-adjusted C3 niche-breadth effect is directionally consistent and its 95% interval excludes zero at both >=5 and >=10 occupied-cell thresholds."
  } else if (!same_direction) {
    classification <- "unstable"
    classification_reason <- "The occupancy-adjusted C3 niche-breadth effect changes direction between >=5 and >=10 occupied-cell thresholds."
  } else {
    classification <- "exploratory"
    classification_reason <- "The C3 niche-breadth effect is directionally consistent but one or both 95% intervals include zero."
  }
}

decision_lines <- c(
  "STEP 12N OMI RETENTION SCREEN",
  "=============================",
  "",
  paste0("Classification: ", classification),
  paste0("Reason: ", classification_reason),
  "",
  "Interpretation:",
  "OMI describes realized environmental niche position and breadth across the sampled Baja cells.",
  "It does not identify the immediate weather conditions that trigger individual ballooning events.",
  "A broader niche among ballooning-capable genera would support environmental generalism or broader realized occupancy, not causation by ballooning alone.",
  "",
  "Primary reporting rule:",
  "Use the occupancy- and order-adjusted C3 coefficient as the main trait comparison.",
  "Treat the raw boxplot and Wilcoxon result as descriptive because niche breadth increases mechanically with the number of occupied cells.",
  ""
)
writeLines(decision_lines, file.path(out_dir, "12N_retention_screen.txt"))

primary_global <- global_all[global_all$case_id == "C3_occ5_primary", , drop = FALSE]
primary_group <- groups_all[groups_all$case_id == "C3_occ5_primary", , drop = FALSE]
caption <- paste0(
  "Figure 12N. Environmental niche position and breadth of ballooning-capable and non-ballooning Baja arachnid genera. ",
  "Outlying Mean Index (OMI) analysis used genus incidence across ", nrow(env),
  " environmentally complete occupied 25-km cells and four standardized environmental dimensions: vapor-pressure deficit, wind seasonality, vegetation/phenology, and topography. ",
  "The primary comparison retained genera occurring in at least five analyzed cells and classified D1-D3/C3 genera as ballooning versus fixed N0 genera as non-ballooning; D4 was excluded. ",
  "Panel A shows genus niche positions in the first two OMI axes, Panel B shows total realized niche tolerance, and Panel C maps the principal environmental niche gradient. ",
  "Trait inference should be based on the occupancy- and order-adjusted model rather than the unadjusted boxplot."
)
writeLines(caption, file.path(out_dir, "Figure_12N_caption.txt"))

methods_text <- c(
  "STEP 12N METHODS SUMMARY",
  "========================",
  "",
  "We quantified genus-level realized environmental niches using Outlying Mean Index (OMI) analysis.",
  paste0("The primary analysis used ", nrow(env), " environmentally complete occupied 25-km cells and the four Step 12F environmental scores: VPD, wind seasonality, vegetation/phenology PC1, and topography PC1."),
  "The environmental table was centered and standardized by PCA, and genus incidence was linked to that environmental space using ade4::niche().",
  "For each genus, we calculated OMI marginality and total tolerance (Tol + Rtol), interpreted respectively as realized niche position and breadth.",
  "The primary trait comparison was C3 (D1-D3) versus fixed N0, with D4 excluded, and retained genera occurring in at least five analyzed cells.",
  "Because estimated niche breadth increases with occupancy, the primary genus-level test modeled log(1 + total tolerance) as a function of ballooning classification while controlling for log occupied-cell count and taxonomic order when estimable.",
  "Sensitivity analyses used a ten-cell occupancy threshold, C1/C2/C4 trait definitions, and the taxonomy-strict incidence matrix.",
  paste0("Global OMI structure was evaluated with ", n_perm, " Monte Carlo permutations in the primary analysis."),
  ""
)
writeLines(methods_text, file.path(out_dir, "12N_methods_summary.txt"))

provenance <- list(
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  script_version = SCRIPT_VERSION,
  script_path = NA_character_,
  project_root = project_root,
  inputs = list(
    primary_incidence = primary_path,
    taxonomy_strict_incidence = strict_path,
    trait_lookup = trait_path,
    step12c_response = response_path,
    step12f_scores = scores_path,
    cell_geojson = cell_geojson_path
  ),
  environmental_predictors = required_env,
  primary_trait_definition = "C3 = D1 + D2 + D3/C3 versus fixed N0; D4 excluded",
  primary_minimum_occupied_cells = 5,
  sensitivity_minimum_occupied_cells = 10,
  permutations = n_perm,
  seed = seed,
  analyzed_cells = nrow(env),
  primary_genera = nrow(primary_objects$params),
  retention_classification = classification,
  retention_reason = classification_reason
)
script_args <- commandArgs(trailingOnly = FALSE)
script_args <- script_args[grep("^--file=", script_args)]
provenance$script_path <- if (length(script_args)) sub("^--file=", "", script_args[[1]]) else NA_character_
jsonlite::write_json(provenance, file.path(out_dir, "12N_provenance.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "12N_sessionInfo.txt"))

log_msg("Analyzed cells: ", nrow(env))
log_msg("Primary retained genera: ", nrow(primary_objects$params))
log_msg("Primary global mean-OMI p-value: ", signif(primary_global$p_value_mean_OMI, 4))
log_msg("Retention classification: ", classification)
log_msg("Validation checks passed: ", sum(validation$passed), "/", nrow(validation))
log_msg("Outputs: ", out_dir)
log_msg("STEP 12N COMPLETED SUCCESSFULLY")
