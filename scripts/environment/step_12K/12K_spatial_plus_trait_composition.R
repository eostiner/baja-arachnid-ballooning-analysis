#!/usr/bin/env Rscript

# ============================================================
# STEP 12K — Spatial+ robustness analysis of arachnid
#            ballooning-capable genus representation
#
# Baja Ballooning Publication
#
# Purpose
#   Test whether the retained VPD and wind-seasonality associations
#   persist after broad spatial structure is removed from the
#   environmental predictors using the Spatial+ approach of
#   Dupont, Wood & Augustin (2022).
#
# Primary response
#   cbind(ballooning genera, non-ballooning genera) at the
#   occupied 25-km cell level.
#
# Primary models
#   1. Null
#   2. Spatial-only GAM
#   3. Nonspatial VPD + wind model
#   4. Spatial GAM with raw VPD + wind
#   5. Spatial+ GAM with residualized VPD + wind
#   6. Spatial+ extended GAM adding residualized vegetation
#      and topographic axes
#
# Key inference
#   Spatial+ coefficients describe associations with local
#   departures from each predictor's broad spatial trend.
#   They are not total environmental effects and are not causal.
#
# The script refits the covariate spatial trends inside every
# leave-one-latitude-band-out fold to prevent information leakage.
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) {
  normalizePath(args[1], mustWork = TRUE)
} else {
  normalizePath("~/Desktop/Baja_Ballooning_Pipeline", mustWork = TRUE)
}
run_mode <- if (length(args) >= 2) tolower(args[2]) else "quick"
seed <- if (length(args) >= 3) as.integer(args[3]) else 20260714L
if (!run_mode %in% c("quick", "paper")) stop("Run mode must be 'quick' or 'paper'.")
if (!is.finite(seed)) seed <- 20260714L
set.seed(seed)

version <- "12K_v2_2026-07-14"
primary_k <- 30L
k_sensitivity <- c(20L, 30L, 40L)
moran_permutations <- if (run_mode == "paper") 999L else 199L
band_levels <- c("23-24N", "24-26N", "26-28N", "28-30N", "30-32N")

required_packages <- c("mgcv", "ggplot2", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
  stop(
    "Missing required R package(s): ", paste(missing_packages, collapse = ", "),
    ". Run 12K_install_packages.R first."
  )
}

step12c_dir <- file.path(project_root, "04_analysis", "12C_cell_environment_model_table")
step12f_dir <- file.path(project_root, "04_analysis", "12F_environment_predictor_scores")
out_dir <- file.path(project_root, "04_analysis", "12K_spatial_plus_trait_composition")
archive_root <- file.path(project_root, "08_archive")

archive_existing <- function(path, archive_root, label) {
  if (!dir.exists(path)) return(invisible(NULL))
  dir.create(archive_root, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y%m%dT%H%M%S")
  destination <- file.path(archive_root, paste0(label, "_", stamp))
  if (file.rename(path, destination)) return(invisible(destination))
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  old_files <- list.files(path, full.names = TRUE, recursive = TRUE, all.files = TRUE,
                          no.. = TRUE, include.dirs = FALSE)
  if (length(old_files)) {
    rel <- substring(old_files, nchar(path) + 2L)
    new_files <- file.path(destination, rel)
    dir.create(unique(dirname(new_files)), recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(old_files, new_files, overwrite = TRUE, copy.mode = TRUE,
                    copy.date = TRUE)
    if (!all(ok)) stop("Failed to archive prior Step 12K output.")
  }
  unlink(path, recursive = TRUE, force = TRUE)
  invisible(destination)
}

archived <- archive_existing(
  out_dir, archive_root, "12K_spatial_plus_trait_composition"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "models"), recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(out_dir, "12K_analysis_log.txt")
log_con <- file(log_path, open = "wt")
on.exit(close(log_con), add = TRUE)
log_msg <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n")
  writeLines(txt, log_con)
  flush(log_con)
}

log_msg("STEP 12K SPATIAL+ STARTED")
log_msg("Version: ", version)
log_msg("Run mode: ", run_mode)
log_msg("Seed: ", seed)
log_msg("Project root: ", project_root)
if (!is.null(archived)) log_msg("Archived prior output: ", archived)
log_msg("mgcv version: ", as.character(utils::packageVersion("mgcv")))

response_path <- file.path(step12c_dir, "12C_primary_glm_candidate_table.csv")
score_path <- file.path(step12f_dir, "12F_environment_predictor_scores_by_cell.csv")

required_inputs <- c(response_path, score_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = "; "))
}

inputs <- data.frame(
  input = c("cell_response_table", "Step12F_environment_scores"),
  path = required_inputs,
  stringsAsFactors = FALSE
)
utils::write.csv(inputs, file.path(out_dir, "12K_input_manifest.csv"), row.names = FALSE)

response <- utils::read.csv(response_path, stringsAsFactors = FALSE, check.names = FALSE)
scores <- utils::read.csv(score_path, stringsAsFactors = FALSE, check.names = FALSE)

required_response <- c(
  "grid_cell_id", "ballooning_genera_primary", "non_ballooning_genera_primary",
  "classified_genera_primary", "ballooning_genera_taxonomy_strict",
  "non_ballooning_genera_taxonomy_strict", "classified_genera_taxonomy_strict",
  "ballooning_genera_low_conf_exclusion",
  "non_ballooning_genera_low_conf_exclusion",
  "classified_genera_low_conf_exclusion"
)
required_scores <- c(
  "grid_cell_id", "centroid_latitude", "centroid_longitude", "latitude_band",
  "easting_km", "northing_km", "vpd_z", "wind_seasonality_z",
  "vegetation_axis", "topography_axis"
)
missing_response <- setdiff(required_response, names(response))
missing_scores <- setdiff(required_scores, names(scores))
if (length(missing_response)) {
  stop("Response table missing column(s): ", paste(missing_response, collapse = "; "))
}
if (length(missing_scores)) {
  stop("Environmental score table missing column(s): ", paste(missing_scores, collapse = "; "))
}
if (anyDuplicated(response$grid_cell_id)) stop("Duplicate grid_cell_id in response table.")
if (anyDuplicated(scores$grid_cell_id)) stop("Duplicate grid_cell_id in environmental score table.")

keep_response <- unique(c("grid_cell_id", required_response[-1]))
d <- merge(
  response[, keep_response, drop = FALSE],
  scores[, required_scores, drop = FALSE],
  by = "grid_cell_id", all = FALSE, sort = FALSE
)
d$latitude_band <- factor(d$latitude_band, levels = band_levels)

core_environment <- c("vpd_z", "wind_seasonality_z", "vegetation_axis", "topography_axis")
complete_core <- complete.cases(d[, c(
  "easting_km", "northing_km", "centroid_latitude", core_environment,
  "ballooning_genera_primary", "non_ballooning_genera_primary",
  "classified_genera_primary"
), drop = FALSE])
d <- d[complete_core & d$classified_genera_primary > 0, , drop = FALSE]
d <- d[order(d$centroid_latitude, d$centroid_longitude), , drop = FALSE]
row.names(d) <- NULL

log_msg("Joined complete primary cells: ", nrow(d))
log_msg("Primary classified genus-cell entries: ", sum(d$classified_genera_primary))
log_msg("Latitude bands: ", paste(levels(droplevels(d$latitude_band)), collapse = "; "))

safe_sd <- function(x) {
  ans <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(ans) || ans == 0) 1 else ans
}

adaptive_k <- function(n, requested = primary_k) {
  upper <- max(8L, floor(n / 4))
  max(6L, min(as.integer(requested), as.integer(upper), n - 2L))
}

clip_prob <- function(p, eps = 1e-8) pmin(pmax(as.numeric(p), eps), 1 - eps)

fit_gam_safe <- function(formula, data, family, method = "REML", label = "") {
  warnings_seen <- character()
  fit <- withCallingHandlers(
    tryCatch(
      mgcv::gam(
        formula = formula, data = data, family = family,
        method = method, na.action = na.exclude
      ),
      error = function(e) structure(list(error = conditionMessage(e)), class = "fit_error")
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (length(warnings_seen)) {
    for (msg in unique(warnings_seen)) log_msg("Warning [", label, "]: ", msg)
  }
  fit
}
fit_glm_safe <- function(formula, data, family, label = "") {
  warnings_seen <- character()
  fit <- withCallingHandlers(
    tryCatch(
      stats::glm(
        formula = formula, data = data, family = family,
        na.action = na.exclude
      ),
      error = function(e) structure(list(error = conditionMessage(e)), class = "fit_error")
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (length(warnings_seen)) {
    for (msg in unique(warnings_seen)) log_msg("Warning [", label, "]: ", msg)
  }
  fit
}
is_fit_error <- function(x) inherits(x, "fit_error")
fit_ok <- function(x) !is.null(x) && !is_fit_error(x) &&
  inherits(x, c("glm", "gam")) && isTRUE(x$converged)

response_formula <- function(success, failure, rhs) {
  stats::as.formula(paste0("cbind(", success, ", ", failure, ") ~ ", rhs))
}

raw_term_names <- c(
  vpd_z = "vpd_rawstd",
  wind_seasonality_z = "wind_seasonality_rawstd",
  vegetation_axis = "vegetation_rawstd",
  topography_axis = "topography_rawstd"
)
spplus_term_names <- c(
  vpd_z = "vpd_spplus",
  wind_seasonality_z = "wind_seasonality_spplus",
  vegetation_axis = "vegetation_spplus",
  topography_axis = "topography_spplus"
)

spatial_plus_transform <- function(train, test = train, predictors = core_environment,
                                   requested_k = primary_k, label = "") {
  train_out <- train
  test_out <- test
  k_eff <- adaptive_k(nrow(train), requested_k)
  decomposition <- list()
  covariate_models <- list()

  for (pred in predictors) {
    mu <- mean(train[[pred]], na.rm = TRUE)
    sig <- safe_sd(train[[pred]])
    raw_name <- unname(raw_term_names[pred])
    sp_name <- unname(spplus_term_names[pred])

    train_out[[raw_name]] <- (train[[pred]] - mu) / sig
    test_out[[raw_name]] <- (test[[pred]] - mu) / sig

    fm <- stats::as.formula(
      paste0(raw_name, " ~ s(easting_km, northing_km, bs='tp', k=", k_eff, ")")
    )
    cov_fit <- fit_gam_safe(
      fm, train_out, stats::gaussian(), method = "REML",
      label = paste0(label, " first-stage ", pred)
    )
    if (!fit_ok(cov_fit)) {
      err <- if (is_fit_error(cov_fit)) cov_fit$error else "unusable first-stage model"
      stop("Spatial+ first-stage model failed for ", pred, ": ", err)
    }

    trend_train <- as.numeric(stats::predict(cov_fit, newdata = train_out, type = "response"))
    trend_test <- as.numeric(stats::predict(cov_fit, newdata = test_out, type = "response"))
    resid_train <- train_out[[raw_name]] - trend_train
    resid_test <- test_out[[raw_name]] - trend_test
    resid_mu <- mean(resid_train, na.rm = TRUE)
    resid_sd <- safe_sd(resid_train)
    train_out[[sp_name]] <- (resid_train - resid_mu) / resid_sd
    test_out[[sp_name]] <- (resid_test - resid_mu) / resid_sd

    sm <- summary(cov_fit)
    total_ss <- sum((train_out[[raw_name]] - mean(train_out[[raw_name]]))^2)
    residual_ss <- sum(resid_train^2)
    spatial_r2 <- if (is.finite(total_ss) && total_ss > 0) {
      1 - residual_ss / total_ss
    } else {
      NA_real_
    }
    s_table <- sm$s.table
    smooth_edf <- if (!is.null(s_table) && nrow(s_table)) s_table[1, "edf"] else NA_real_
    smooth_p <- if (!is.null(s_table) && nrow(s_table)) s_table[1, ncol(s_table)] else NA_real_

    decomposition[[pred]] <- data.frame(
      predictor = pred,
      raw_term = raw_name,
      spatial_plus_term = sp_name,
      n_training_cells = nrow(train_out),
      basis_dimension_k = k_eff,
      spatial_r_squared = spatial_r2,
      deviance_explained = sm$dev.expl,
      smooth_edf = smooth_edf,
      smooth_p_value = smooth_p,
      residual_sd_before_standardization = resid_sd,
      raw_trend_correlation = suppressWarnings(stats::cor(
        train_out[[raw_name]], trend_train, use = "complete.obs"
      )),
      raw_residual_correlation = suppressWarnings(stats::cor(
        train_out[[raw_name]], resid_train, use = "complete.obs"
      )),
      stringsAsFactors = FALSE
    )
    covariate_models[[pred]] <- cov_fit
  }

  list(
    train = train_out,
    test = test_out,
    decomposition = do.call(rbind, decomposition),
    covariate_models = covariate_models,
    k = k_eff
  )
}

model_definitions <- function(k_eff) {
  smooth_rhs <- paste0("s(easting_km, northing_km, bs='tp', k=", k_eff, ")")
  list(
    null = list(engine = "glm", rhs = "1"),
    spatial_only = list(engine = "gam", rhs = smooth_rhs),
    nonspatial_core = list(
      engine = "glm",
      rhs = "vpd_rawstd + wind_seasonality_rawstd"
    ),
    spatial_raw_core = list(
      engine = "gam",
      rhs = paste("vpd_rawstd + wind_seasonality_rawstd +", smooth_rhs)
    ),
    spatial_plus_core = list(
      engine = "gam",
      rhs = paste("vpd_spplus + wind_seasonality_spplus +", smooth_rhs)
    ),
    spatial_raw_extended = list(
      engine = "gam",
      rhs = paste(
        "vpd_rawstd + wind_seasonality_rawstd + vegetation_rawstd +",
        "topography_rawstd +", smooth_rhs
      )
    ),
    spatial_plus_extended = list(
      engine = "gam",
      rhs = paste(
        "vpd_spplus + wind_seasonality_spplus + vegetation_spplus +",
        "topography_spplus +", smooth_rhs
      )
    )
  )
}

fit_outcome_model <- function(model_name, data, success, failure,
                              family = stats::binomial(link = "logit"),
                              requested_k = primary_k, label = "") {
  k_eff <- adaptive_k(nrow(data), requested_k)
  defs <- model_definitions(k_eff)
  if (!model_name %in% names(defs)) stop("Unknown model: ", model_name)
  spec <- defs[[model_name]]
  fm <- response_formula(success, failure, spec$rhs)
  if (spec$engine == "gam") {
    fit_gam_safe(fm, data, family, method = "REML",
                 label = paste(label, model_name))
  } else {
    fit_glm_safe(fm, data, family,
                 label = paste(label, model_name))
  }
}

pearson_dispersion <- function(fit) {
  if (!fit_ok(fit)) return(NA_real_)
  r <- stats::residuals(fit, type = "pearson")
  rdf <- stats::df.residual(fit)
  if (!is.finite(rdf) || rdf <= 0) return(NA_real_)
  sum(r^2, na.rm = TRUE) / rdf
}

deviance_explained <- function(fit) {
  if (!fit_ok(fit)) return(NA_real_)
  if (inherits(fit, "gam")) {
    return(summary(fit)$dev.expl)
  }
  if (!is.finite(fit$null.deviance) || fit$null.deviance <= 0) return(NA_real_)
  1 - fit$deviance / fit$null.deviance
}

extract_parametric_coefficients <- function(fit, model_name, dataset = "primary",
                                            response = "primary") {
  if (!fit_ok(fit)) return(data.frame())
  if (inherits(fit, "gam")) {
    tab <- summary(fit)$p.table
  } else {
    tab <- summary(fit)$coefficients
  }
  if (is.null(tab) || !nrow(tab)) return(data.frame())
  p_col <- grep("^Pr\\(", colnames(tab), value = TRUE)
  stat_col <- setdiff(colnames(tab), c("Estimate", "Std. Error", p_col))
  out <- data.frame(
    dataset = dataset,
    response = response,
    model = model_name,
    term = rownames(tab),
    estimate = tab[, "Estimate"],
    std_error = tab[, "Std. Error"],
    statistic = if (length(stat_col)) tab[, stat_col[1]] else NA_real_,
    p_value = if (length(p_col)) tab[, p_col[1]] else NA_real_,
    stringsAsFactors = FALSE
  )
  out$odds_ratio <- exp(out$estimate)
  out$ci_low <- exp(out$estimate - 1.96 * out$std_error)
  out$ci_high <- exp(out$estimate + 1.96 * out$std_error)
  out
}

model_comparison_row <- function(fit, model_name, n_cells) {
  if (!fit_ok(fit)) {
    return(data.frame(
      model = model_name, n_cells = n_cells, fitted = FALSE,
      AIC = NA_real_, deviance_explained_pct = NA_real_,
      pearson_dispersion = NA_real_,
      error = if (is_fit_error(fit)) fit$error else "unusable fit",
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    model = model_name,
    n_cells = n_cells,
    fitted = TRUE,
    AIC = suppressWarnings(tryCatch(stats::AIC(fit), error = function(e) NA_real_)),
    deviance_explained_pct = 100 * deviance_explained(fit),
    pearson_dispersion = pearson_dispersion(fit),
    error = "",
    stringsAsFactors = FALSE
  )
}

extract_smooth_summary <- function(fit, model_name) {
  if (!fit_ok(fit) || !inherits(fit, "gam")) return(data.frame())
  tab <- summary(fit)$s.table
  if (is.null(tab) || !nrow(tab)) return(data.frame())
  data.frame(
    model = model_name,
    smooth = rownames(tab),
    edf = tab[, "edf"],
    reference_df = if ("Ref.df" %in% colnames(tab)) tab[, "Ref.df"] else NA_real_,
    statistic = tab[, ncol(tab) - 1L],
    p_value = tab[, ncol(tab)],
    stringsAsFactors = FALSE
  )
}

extract_concurvity <- function(fit, model_name) {
  if (!fit_ok(fit) || !inherits(fit, "gam")) return(data.frame())
  cc <- tryCatch(mgcv::concurvity(fit, full = TRUE), error = function(e) NULL)
  if (is.null(cc)) return(data.frame())
  if (is.matrix(cc)) {
    out <- as.data.frame(as.table(cc), stringsAsFactors = FALSE)
    names(out) <- c("measure", "term", "concurvity")
    out$model <- model_name
    out[, c("model", "measure", "term", "concurvity")]
  } else {
    data.frame()
  }
}

moran_permutation <- function(residuals, x, y, neighbors = 4L,
                              permutations = moran_permutations) {
  good <- is.finite(residuals) & is.finite(x) & is.finite(y)
  z <- residuals[good]
  coords <- cbind(x[good], y[good])
  n <- length(z)
  if (n < 10) {
    return(data.frame(
      n = n, neighbors = NA_integer_, moran_i = NA_real_,
      permutation_p_value = NA_real_, permutations = permutations
    ))
  }
  k <- min(as.integer(neighbors), n - 1L)
  distances <- as.matrix(stats::dist(coords))
  diag(distances) <- Inf
  W <- matrix(0, n, n)
  for (i in seq_len(n)) {
    nn <- order(distances[i, ])[seq_len(k)]
    W[i, nn] <- 1
  }
  W <- (W + t(W) > 0) * 1
  rs <- rowSums(W)
  W[rs > 0, ] <- W[rs > 0, , drop = FALSE] / rs[rs > 0]
  z <- z - mean(z)
  s0 <- sum(W)
  denom <- sum(z^2)
  moran_stat <- function(v) {
    if (!is.finite(denom) || denom == 0 || s0 == 0) return(NA_real_)
    (n / s0) * sum(W * outer(v, v)) / sum(v^2)
  }
  observed <- moran_stat(z)
  permuted <- replicate(permutations, moran_stat(sample(z, replace = FALSE)))
  p <- (1 + sum(abs(permuted) >= abs(observed), na.rm = TRUE)) /
    (1 + sum(is.finite(permuted)))
  data.frame(
    n = n, neighbors = k, moran_i = observed,
    permutation_p_value = p, permutations = permutations
  )
}

# ------------------------------------------------------------
# Full-data Spatial+ transformation and model fitting
# ------------------------------------------------------------
sp_full <- spatial_plus_transform(
  train = d, test = d, predictors = core_environment,
  requested_k = primary_k, label = "full data"
)
analysis_data <- sp_full$train
decomposition <- sp_full$decomposition
utils::write.csv(
  decomposition,
  file.path(out_dir, "12K_covariate_spatial_decomposition.csv"),
  row.names = FALSE
)

primary_model_names <- names(model_definitions(sp_full$k))
primary_fits <- list()
comparison_rows <- list()
coef_rows <- list()
smooth_rows <- list()
concurvity_rows <- list()
moran_rows <- list()

for (nm in primary_model_names) {
  fit <- fit_outcome_model(
    nm, analysis_data,
    "ballooning_genera_primary", "non_ballooning_genera_primary",
    stats::binomial(link = "logit"), primary_k,
    label = "full primary"
  )
  primary_fits[[nm]] <- fit
  comparison_rows[[nm]] <- model_comparison_row(fit, nm, nrow(analysis_data))
  cc <- extract_parametric_coefficients(fit, nm)
  if (nrow(cc)) coef_rows[[nm]] <- cc
  ss <- extract_smooth_summary(fit, nm)
  if (nrow(ss)) smooth_rows[[nm]] <- ss
  cv <- extract_concurvity(fit, nm)
  if (nrow(cv)) concurvity_rows[[nm]] <- cv

  mr <- if (fit_ok(fit)) {
    moran_permutation(
      stats::residuals(fit, type = "pearson"),
      analysis_data$easting_km, analysis_data$northing_km
    )
  } else {
    data.frame(
      n = nrow(analysis_data), neighbors = NA_integer_,
      moran_i = NA_real_, permutation_p_value = NA_real_,
      permutations = moran_permutations
    )
  }
  mr$model <- nm
  moran_rows[[nm]] <- mr
}

model_comparison <- do.call(rbind, comparison_rows)
if (any(is.finite(model_comparison$AIC))) {
  model_comparison$delta_AIC <- model_comparison$AIC -
    min(model_comparison$AIC, na.rm = TRUE)
} else {
  model_comparison$delta_AIC <- NA_real_
}
model_coefficients <- if (length(coef_rows)) do.call(rbind, coef_rows) else data.frame()
smooth_summary <- if (length(smooth_rows)) do.call(rbind, smooth_rows) else data.frame()
concurvity_summary <- if (length(concurvity_rows)) {
  do.call(rbind, concurvity_rows)
} else {
  data.frame()
}
moran_summary <- do.call(rbind, moran_rows)

utils::write.csv(model_comparison,
                 file.path(out_dir, "12K_primary_model_comparison.csv"),
                 row.names = FALSE)
utils::write.csv(model_coefficients,
                 file.path(out_dir, "12K_primary_model_coefficients.csv"),
                 row.names = FALSE)
utils::write.csv(smooth_summary,
                 file.path(out_dir, "12K_spatial_smooth_summary.csv"),
                 row.names = FALSE)
utils::write.csv(concurvity_summary,
                 file.path(out_dir, "12K_concurvity_comparison.csv"),
                 row.names = FALSE)
utils::write.csv(moran_summary,
                 file.path(out_dir, "12K_residual_moran_permutation.csv"),
                 row.names = FALSE)

# Quasibinomial uncertainty companion.
quasi_fit <- fit_outcome_model(
  "spatial_plus_core", analysis_data,
  "ballooning_genera_primary", "non_ballooning_genera_primary",
  stats::quasibinomial(link = "logit"), primary_k,
  label = "quasibinomial companion"
)
quasi_coef <- extract_parametric_coefficients(
  quasi_fit, "spatial_plus_core_quasibinomial"
)
if (nrow(quasi_coef)) quasi_coef$pearson_dispersion <- pearson_dispersion(quasi_fit)
utils::write.csv(
  quasi_coef,
  file.path(out_dir, "12K_spatial_plus_quasibinomial_coefficients.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Leave-one-latitude-band-out cross-validation
# Spatial+ covariate models are refitted inside each fold.
# ------------------------------------------------------------
cv_model_names <- c(
  "null", "spatial_only", "nonspatial_core", "spatial_raw_core",
  "spatial_plus_core", "spatial_plus_extended"
)
cv_rows <- list()
cv_prediction_rows <- list()

for (held_out in band_levels) {
  train0 <- d[as.character(d$latitude_band) != held_out, , drop = FALSE]
  test0 <- d[as.character(d$latitude_band) == held_out, , drop = FALSE]
  if (!nrow(train0) || !nrow(test0)) next

  fold_sp <- spatial_plus_transform(
    train0, test0, core_environment, primary_k,
    label = paste0("CV held out ", held_out)
  )
  train <- fold_sp$train
  test <- fold_sp$test

  for (nm in cv_model_names) {
    fit <- fit_outcome_model(
      nm, train,
      "ballooning_genera_primary", "non_ballooning_genera_primary",
      stats::binomial(link = "logit"), primary_k,
      label = paste0("CV ", held_out)
    )
    if (!fit_ok(fit)) {
      cv_rows[[length(cv_rows) + 1L]] <- data.frame(
        model = nm, held_out_band = held_out,
        n_cells = nrow(test),
        total_classified_genera = sum(test$classified_genera_primary),
        weighted_brier = NA_real_, weighted_mae = NA_real_,
        binomial_log_score_per_genus = NA_real_,
        error = if (is_fit_error(fit)) fit$error else "unusable fit",
        stringsAsFactors = FALSE
      )
      next
    }

    pred <- tryCatch(
      clip_prob(stats::predict(fit, newdata = test, type = "response")),
      error = function(e) rep(NA_real_, nrow(test))
    )
    y <- test$ballooning_genera_primary
    n <- test$classified_genera_primary
    obs <- y / n
    valid <- is.finite(pred) & is.finite(obs) & is.finite(n) & n > 0
    w <- n[valid]

    cv_rows[[length(cv_rows) + 1L]] <- data.frame(
      model = nm,
      held_out_band = held_out,
      n_cells = sum(valid),
      total_classified_genera = sum(w),
      weighted_brier = if (sum(valid)) {
        sum(w * (pred[valid] - obs[valid])^2) / sum(w)
      } else NA_real_,
      weighted_mae = if (sum(valid)) {
        sum(w * abs(pred[valid] - obs[valid])) / sum(w)
      } else NA_real_,
      binomial_log_score_per_genus = if (sum(valid)) {
        -sum(stats::dbinom(
          y[valid], size = n[valid], prob = pred[valid], log = TRUE
        )) / sum(w)
      } else NA_real_,
      error = "",
      stringsAsFactors = FALSE
    )

    cv_prediction_rows[[length(cv_prediction_rows) + 1L]] <- data.frame(
      grid_cell_id = test$grid_cell_id,
      latitude_band = as.character(test$latitude_band),
      model = nm,
      observed_ballooning = y,
      denominator = n,
      observed_proportion = obs,
      predicted_proportion = pred,
      stringsAsFactors = FALSE
    )
  }
}

cv_by_fold <- do.call(rbind, cv_rows)
cv_predictions <- if (length(cv_prediction_rows)) {
  do.call(rbind, cv_prediction_rows)
} else {
  data.frame()
}
cv_summary <- do.call(rbind, lapply(
  split(cv_by_fold, cv_by_fold$model),
  function(x) {
    valid <- is.finite(x$weighted_brier) &
      is.finite(x$binomial_log_score_per_genus) &
      is.finite(x$total_classified_genera) &
      x$total_classified_genera > 0
    w <- x$total_classified_genera[valid]
    data.frame(
      model = unique(x$model),
      folds_completed = sum(valid),
      total_classified_genera = sum(w),
      weighted_brier = if (length(w)) {
        stats::weighted.mean(x$weighted_brier[valid], w)
      } else NA_real_,
      weighted_mae = if (length(w)) {
        stats::weighted.mean(x$weighted_mae[valid], w)
      } else NA_real_,
      binomial_log_score_per_genus = if (length(w)) {
        stats::weighted.mean(x$binomial_log_score_per_genus[valid], w)
      } else NA_real_,
      stringsAsFactors = FALSE
    )
  }
))
cv_summary <- cv_summary[match(cv_model_names, cv_summary$model), , drop = FALSE]
null_cv <- cv_summary[cv_summary$model == "null", , drop = FALSE]
if (nrow(null_cv)) {
  cv_summary$delta_log_score_vs_null <-
    cv_summary$binomial_log_score_per_genus -
    null_cv$binomial_log_score_per_genus
  cv_summary$delta_brier_vs_null <-
    cv_summary$weighted_brier - null_cv$weighted_brier
} else {
  cv_summary$delta_log_score_vs_null <- NA_real_
  cv_summary$delta_brier_vs_null <- NA_real_
}

utils::write.csv(cv_by_fold,
                 file.path(out_dir, "12K_leave_one_band_out_cv_by_fold.csv"),
                 row.names = FALSE)
utils::write.csv(cv_predictions,
                 file.path(out_dir, "12K_leave_one_band_out_cv_predictions.csv"),
                 row.names = FALSE)
utils::write.csv(cv_summary,
                 file.path(out_dir, "12K_leave_one_band_out_cv_summary.csv"),
                 row.names = FALSE)

# ------------------------------------------------------------
# Response and denominator sensitivity.
# Environmental residuals use the full current primary-cell decomposition
# so the coefficient scales remain comparable.
# ------------------------------------------------------------
response_specs <- data.frame(
  response = c("primary", "taxonomy_strict", "low_conf_exclusion"),
  success = c(
    "ballooning_genera_primary",
    "ballooning_genera_taxonomy_strict",
    "ballooning_genera_low_conf_exclusion"
  ),
  failure = c(
    "non_ballooning_genera_primary",
    "non_ballooning_genera_taxonomy_strict",
    "non_ballooning_genera_low_conf_exclusion"
  ),
  denominator = c(
    "classified_genera_primary",
    "classified_genera_taxonomy_strict",
    "classified_genera_low_conf_exclusion"
  ),
  stringsAsFactors = FALSE
)
thresholds <- c(all_positive = 1L, denominator_ge5 = 5L, denominator_ge10 = 10L)
sensitivity_rows <- list()

for (i in seq_len(nrow(response_specs))) {
  rs <- response_specs[i, ]
  for (threshold_name in names(thresholds)) {
    threshold <- thresholds[[threshold_name]]
    keep <- is.finite(analysis_data[[rs$denominator]]) &
      analysis_data[[rs$denominator]] >= threshold
    ds <- analysis_data[keep, , drop = FALSE]
    fit <- fit_outcome_model(
      "spatial_plus_core", ds, rs$success, rs$failure,
      stats::binomial(link = "logit"), primary_k,
      label = paste("sensitivity", rs$response, threshold_name)
    )
    cc <- extract_parametric_coefficients(
      fit, "spatial_plus_core", threshold_name, rs$response
    )
    cc <- cc[cc$term %in% c("vpd_spplus", "wind_seasonality_spplus"), , drop = FALSE]
    if (nrow(cc)) {
      cc$n_cells <- nrow(ds)
      cc$denominator_threshold <- threshold
      cc$converged <- fit_ok(fit)
      sensitivity_rows[[length(sensitivity_rows) + 1L]] <- cc
    } else {
      sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.frame(
        dataset = threshold_name, response = rs$response,
        model = "spatial_plus_core",
        term = c("vpd_spplus", "wind_seasonality_spplus"),
        estimate = NA_real_, std_error = NA_real_,
        statistic = NA_real_, p_value = NA_real_,
        odds_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        n_cells = nrow(ds), denominator_threshold = threshold,
        converged = FALSE, stringsAsFactors = FALSE
      )
    }
  }
}
sensitivity_coefficients <- do.call(rbind, sensitivity_rows)
utils::write.csv(
  sensitivity_coefficients,
  file.path(out_dir, "12K_spatial_plus_key_effect_sensitivity.csv"),
  row.names = FALSE
)

stability_summary <- do.call(rbind, lapply(
  split(sensitivity_coefficients, sensitivity_coefficients$term),
  function(x) {
    x <- x[is.finite(x$estimate), , drop = FALSE]
    reference <- x$estimate[
      x$dataset == "all_positive" & x$response == "primary"
    ]
    reference <- if (length(reference)) reference[1] else NA_real_
    data.frame(
      term = unique(x$term),
      fits_with_estimate = nrow(x),
      negative_in_all_fits = nrow(x) > 0 && all(x$estimate < 0),
      same_direction_as_primary = nrow(x) > 0 && is.finite(reference) &&
        all(sign(x$estimate) == sign(reference)),
      confidence_interval_excludes_one = sum(
        x$ci_high < 1 | x$ci_low > 1, na.rm = TRUE
      ),
      min_odds_ratio = if (nrow(x)) min(x$odds_ratio, na.rm = TRUE) else NA_real_,
      max_odds_ratio = if (nrow(x)) max(x$odds_ratio, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
))
utils::write.csv(
  stability_summary,
  file.path(out_dir, "12K_spatial_plus_key_effect_stability_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Basis-dimension sensitivity for the primary Spatial+ model.
# ------------------------------------------------------------
k_rows <- list()
for (kk in k_sensitivity) {
  sp_k <- spatial_plus_transform(
    d, d, core_environment, requested_k = kk,
    label = paste0("k sensitivity ", kk)
  )
  fit <- fit_outcome_model(
    "spatial_plus_core", sp_k$train,
    "ballooning_genera_primary", "non_ballooning_genera_primary",
    stats::binomial(link = "logit"), requested_k = kk,
    label = paste0("k sensitivity ", kk)
  )
  cc <- extract_parametric_coefficients(
    fit, paste0("spatial_plus_core_k", kk)
  )
  cc <- cc[cc$term %in% c("vpd_spplus", "wind_seasonality_spplus"), , drop = FALSE]
  if (nrow(cc)) {
    cc$requested_k <- kk
    cc$actual_k <- sp_k$k
    cc$deviance_explained_pct <- 100 * deviance_explained(fit)
    k_rows[[length(k_rows) + 1L]] <- cc
  }
}
k_sensitivity_table <- if (length(k_rows)) do.call(rbind, k_rows) else data.frame()
utils::write.csv(
  k_sensitivity_table,
  file.path(out_dir, "12K_spatial_basis_dimension_sensitivity.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Figures
# ------------------------------------------------------------
library(ggplot2)
library(patchwork)

term_labels <- c(
  vpd_rawstd = "Vapor-pressure deficit",
  vpd_spplus = "Vapor-pressure deficit",
  wind_seasonality_rawstd = "Wind seasonality",
  wind_seasonality_spplus = "Wind seasonality",
  vegetation_rawstd = "Vegetation structure",
  vegetation_spplus = "Vegetation structure",
  topography_rawstd = "Topographic heterogeneity",
  topography_spplus = "Topographic heterogeneity"
)
model_labels <- c(
  nonspatial_core = "Nonspatial GLM",
  spatial_raw_core = "Spatial GAM, raw predictors",
  spatial_plus_core = "Spatial+ GAM"
)

forest <- model_coefficients[
  model_coefficients$model %in% names(model_labels) &
    model_coefficients$term %in% c(
      "vpd_rawstd", "vpd_spplus",
      "wind_seasonality_rawstd", "wind_seasonality_spplus"
    ),
  , drop = FALSE
]
forest$predictor <- unname(term_labels[forest$term])
forest$model_label <- unname(model_labels[forest$model])
forest$model_label <- factor(
  forest$model_label,
  levels = rev(unname(model_labels))
)
forest$predictor <- factor(
  forest$predictor,
  levels = c("Wind seasonality", "Vapor-pressure deficit")
)

p_effect <- ggplot(forest, aes(
  x = odds_ratio, y = predictor, shape = model_label
)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.45) +
  geom_segment(aes(
    x = ci_low, xend = ci_high, y = predictor, yend = predictor
  ), position = position_dodge(width = 0.55), linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.55), size = 2.8) +
  scale_x_log10() +
  labs(
    x = "Odds ratio per 1 SD increase",
    y = NULL,
    shape = NULL,
    title = "Environmental associations before and after Spatial+"
  ) +
  theme_bw(base_size = 10.5) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

cv_plot_data <- cv_summary[is.finite(cv_summary$binomial_log_score_per_genus), , drop = FALSE]
cv_plot_data$model_label <- c(
  null = "Null",
  spatial_only = "Spatial only",
  nonspatial_core = "VPD + wind GLM",
  spatial_raw_core = "Raw spatial GAM",
  spatial_plus_core = "Spatial+ core",
  spatial_plus_extended = "Spatial+ extended"
)[cv_plot_data$model]
cv_plot_data$model_label <- factor(
  cv_plot_data$model_label,
  levels = rev(cv_plot_data$model_label[
    order(cv_plot_data$binomial_log_score_per_genus)
  ])
)

fold_plot <- cv_by_fold[
  is.finite(cv_by_fold$binomial_log_score_per_genus) &
    cv_by_fold$model %in% cv_plot_data$model,
  , drop = FALSE
]
fold_plot$model_label <- c(
  null = "Null",
  spatial_only = "Spatial only",
  nonspatial_core = "VPD + wind GLM",
  spatial_raw_core = "Raw spatial GAM",
  spatial_plus_core = "Spatial+ core",
  spatial_plus_extended = "Spatial+ extended"
)[fold_plot$model]
fold_plot$model_label <- factor(
  fold_plot$model_label,
  levels = levels(cv_plot_data$model_label)
)

p_cv <- ggplot() +
  geom_point(
    data = fold_plot,
    aes(x = binomial_log_score_per_genus, y = model_label),
    alpha = 0.35, size = 1.8, position = position_jitter(height = 0.08)
  ) +
  geom_point(
    data = cv_plot_data,
    aes(x = binomial_log_score_per_genus, y = model_label),
    size = 3
  ) +
  labs(
    x = "Leave-one-band-out log score per genus\n(lower is better)",
    y = NULL,
    title = "Spatial transferability"
  ) +
  theme_bw(base_size = 10.5) +
  theme(panel.grid.minor = element_blank())

main_figure <- p_effect | p_cv +
  plot_annotation(tag_levels = "A")

fig_dir <- file.path(out_dir, "figures")
ggsave(
  file.path(fig_dir, "Figure_3_spatial_plus_environmental_effects.png"),
  main_figure, width = 11.2, height = 5.8, units = "in", dpi = 500,
  bg = "white"
)
ggsave(
  file.path(fig_dir, "Figure_3_spatial_plus_environmental_effects.pdf"),
  main_figure, width = 11.2, height = 5.8, units = "in",
  device = "pdf", bg = "white"
)
ggsave(
  file.path(fig_dir, "Figure_3_spatial_plus_environmental_effects.svg"),
  main_figure, width = 11.2, height = 5.8, units = "in",
  device = "svg", bg = "white"
)
ggsave(
  file.path(fig_dir, "Figure_3_spatial_plus_environmental_effects.tif"),
  main_figure, width = 11.2, height = 5.8, units = "in",
  dpi = 600, compression = "lzw", bg = "white"
)

# Original spatial trends and Spatial+ residuals for the key predictors.
map_data <- rbind(
  data.frame(
    analysis_data[, c("easting_km", "northing_km")],
    predictor = "VPD — raw standardized",
    value = analysis_data$vpd_rawstd
  ),
  data.frame(
    analysis_data[, c("easting_km", "northing_km")],
    predictor = "VPD — Spatial+ residual",
    value = analysis_data$vpd_spplus
  ),
  data.frame(
    analysis_data[, c("easting_km", "northing_km")],
    predictor = "Wind seasonality — raw standardized",
    value = analysis_data$wind_seasonality_rawstd
  ),
  data.frame(
    analysis_data[, c("easting_km", "northing_km")],
    predictor = "Wind seasonality — Spatial+ residual",
    value = analysis_data$wind_seasonality_spplus
  )
)
p_maps <- ggplot(map_data, aes(easting_km, northing_km, color = value)) +
  geom_point(size = 2.2) +
  facet_wrap(~ predictor, ncol = 2) +
  scale_color_gradient2(midpoint = 0) +
  coord_equal() +
  labs(
    x = "Easting (km)", y = "Northing (km)",
    color = "Standardized\nvalue",
    title = "Spatial decomposition of the retained environmental predictors"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())

ggsave(
  file.path(fig_dir, "Figure_S12K_spatial_decomposition.png"),
  p_maps, width = 9.5, height = 8, units = "in", dpi = 400, bg = "white"
)
ggsave(
  file.path(fig_dir, "Figure_S12K_spatial_decomposition.pdf"),
  p_maps, width = 9.5, height = 8, units = "in",
  device = "pdf", bg = "white"
)

if (nrow(k_sensitivity_table)) {
  k_plot <- k_sensitivity_table
  k_plot$predictor <- unname(term_labels[k_plot$term])
  k_plot$k_label <- paste0("k = ", k_plot$requested_k)
  p_k <- ggplot(k_plot, aes(
    x = odds_ratio, y = k_label, shape = predictor
  )) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.45) +
    geom_segment(aes(
      x = ci_low, xend = ci_high, y = k_label, yend = k_label
    ), position = position_dodge(width = 0.45), linewidth = 0.7) +
    geom_point(position = position_dodge(width = 0.45), size = 2.7) +
    scale_x_log10() +
    labs(
      x = "Spatial+ odds ratio per 1 SD residual increase",
      y = NULL, shape = NULL,
      title = "Sensitivity to spatial smooth basis dimension"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  ggsave(
    file.path(fig_dir, "Figure_S12K_basis_dimension_sensitivity.png"),
    p_k, width = 8, height = 4.5, units = "in", dpi = 400, bg = "white"
  )
}

# ------------------------------------------------------------
# Recommendation and manuscript text.
# ------------------------------------------------------------
key_sp <- model_coefficients[
  model_coefficients$model == "spatial_plus_core" &
    model_coefficients$term %in% c("vpd_spplus", "wind_seasonality_spplus"),
  , drop = FALSE
]
both_estimated <- nrow(key_sp) == 2 && all(is.finite(key_sp$estimate))
both_negative <- both_estimated && all(key_sp$estimate < 0)
both_exclude_one <- both_estimated && all(key_sp$ci_high < 1)

robustness_class <- if (both_negative && both_exclude_one) {
  "Strongly robust after Spatial+: both effects remain negative and both 95% confidence intervals exclude an odds ratio of 1."
} else if (both_negative) {
  "Directionally robust but more uncertain after Spatial+: both effects remain negative, but one or both confidence intervals include an odds ratio of 1."
} else if (both_estimated && sum(key_sp$estimate < 0) == 1) {
  "Partially robust after Spatial+: only one retained environmental effect remains negative."
} else {
  "Not robust after Spatial+: the retained effects are absent or reverse after broad spatial structure is removed."
}

sp_cv <- cv_summary[cv_summary$model == "spatial_plus_core", , drop = FALSE]
raw_cv <- cv_summary[cv_summary$model == "spatial_raw_core", , drop = FALSE]
nonsp_cv <- cv_summary[cv_summary$model == "nonspatial_core", , drop = FALSE]
best_cv_model <- if (any(is.finite(cv_summary$binomial_log_score_per_genus))) {
  cv_summary$model[which.min(cv_summary$binomial_log_score_per_genus)]
} else {
  NA_character_
}

decomp_key <- decomposition[
  decomposition$predictor %in% c("vpd_z", "wind_seasonality_z"),
  , drop = FALSE
]

recommendation <- c(
  "STEP 12K SPATIAL+ MODEL RECOMMENDATION",
  "=======================================",
  "",
  paste0("Primary classification: ", robustness_class),
  paste0("Best leave-one-latitude-band-out model by log score: ", best_cv_model),
  "",
  "Spatial decomposition of key predictors:",
  paste0(
    "  VPD broad spatial component R-squared: ",
    signif(decomp_key$spatial_r_squared[decomp_key$predictor == "vpd_z"], 4)
  ),
  paste0(
    "  Wind-seasonality broad spatial component R-squared: ",
    signif(decomp_key$spatial_r_squared[
      decomp_key$predictor == "wind_seasonality_z"
    ], 4)
  ),
  "",
  "Interpretation:",
  "  Spatial+ coefficients quantify associations with local departures from broad spatial environmental trends.",
  "  A retained negative effect indicates that cells locally drier or more wind-seasonal than expected from their position contain a lower relative representation of ballooning-capable genera.",
  "  This is a community-composition association, not evidence about individual ballooning events or a causal environmental mechanism.",
  "  Repeated genera across cells mean that genus-cell entries should not be interpreted as independent evolutionary trials.",
  "",
  "Decision rule for the manuscript:",
  if (both_negative) {
    "  Retain VPD and wind seasonality as spatially deconfounded environmental associations, while reporting confidence intervals and blocked cross-validation."
  } else {
    "  Reframe the earlier VPD/wind result as primarily reflecting broad spatial gradients; do not present it as an independent local environmental association."
  },
  "",
  "Figure use:",
  "  Figure_3_spatial_plus_environmental_effects is the candidate replacement for the original latitude-band GLM and heat map.",
  "  The fine-scale prediction heat map should not be generated from Spatial+ residual effects without a separate prediction objective and validation."
)
writeLines(recommendation, file.path(out_dir, "12K_model_recommendation.txt"))

caption <- paste(
  "Figure 3. Spatially deconfounded environmental associations with the relative",
  "representation of ballooning-capable arachnid genera across occupied 25-km cells",
  "of the Baja California Peninsula. (A) Odds ratios and 95% confidence intervals",
  "for vapor-pressure deficit and wind seasonality from the nonspatial binomial model,",
  "a spatial generalized additive model using the original predictors, and a Spatial+",
  "model in which broad spatial structure was first removed from each predictor.",
  "Spatial+ estimates represent a one-standard-deviation increase in the local",
  "predictor departure from its fitted broad spatial trend. (B) Leave-one-latitude-band-out",
  "binomial log scores for candidate models; smaller values indicate better transfer to",
  "withheld regions. Faint points show individual held-out bands and larger points show",
  "denominator-weighted means. The response was the number of ballooning-capable versus",
  "non-ballooning genera per cell. Associations are interpreted at the community level",
  "and do not measure individual ballooning behavior.",
  sep = " "
)
writeLines(caption, file.path(out_dir, "Figure_3_publication_caption.txt"))

methods_text <- c(
  "SPATIAL+ METHODS REPLACEMENT TEXT",
  "=================================",
  "",
  paste(
    "To evaluate whether environmental associations with ballooning-capable genus",
    "representation were independent of broad spatial gradients, we applied the",
    "Spatial+ approach. Each standardized environmental predictor was first modeled",
    "as a two-dimensional thin-plate regression spline of projected cell coordinates",
    "using restricted maximum likelihood. Residuals from these covariate models",
    "represent local departures from the broad spatial trend and were standardized",
    "before inclusion in a binomial generalized additive model of ballooning-capable",
    "and non-ballooning genus counts. The outcome model retained a two-dimensional",
    "spatial smooth, thereby separating covariate associations from residual spatial",
    "structure. The primary model included vapor-pressure deficit and wind seasonality;",
    "an extended model also included vegetation and topographic axes. Spatial",
    "transferability was assessed by leaving out each latitude band in turn. Within",
    "each fold, covariate spatial trends, residual scaling, and the outcome model were",
    "estimated using training cells only before prediction to the withheld band."
  ),
  "",
  "RESULTS TEMPLATE",
  "================",
  "",
  paste(
    "Broad spatial trends explained",
    paste0(
      round(100 * decomp_key$spatial_r_squared[
        decomp_key$predictor == "vpd_z"
      ], 1), "%"
    ),
    "of standardized VPD variation and",
    paste0(
      round(100 * decomp_key$spatial_r_squared[
        decomp_key$predictor == "wind_seasonality_z"
      ], 1), "%"
    ),
    "of wind-seasonality variation.",
    robustness_class
  ),
  "",
  paste0(
    "The best leave-one-band-out model by log score was ",
    best_cv_model, "."
  ),
  "",
  "The final manuscript wording should be completed after reviewing the coefficient, sensitivity, residual Moran, and cross-validation tables."
)
writeLines(methods_text, file.path(out_dir, "12K_manuscript_replacement_guide.txt"))

saveRDS(
  list(
    version = version,
    run_mode = run_mode,
    primary_k = primary_k,
    analysis_data = analysis_data,
    covariate_models = sp_full$covariate_models,
    outcome_models = primary_fits,
    quasi_model = quasi_fit,
    decomposition = decomposition,
    model_comparison = model_comparison,
    coefficients = model_coefficients,
    cv_summary = cv_summary,
    sensitivity_coefficients = sensitivity_coefficients
  ),
  file.path(out_dir, "models", "12K_spatial_plus_model_objects.rds")
)

# ------------------------------------------------------------
# Validation.
# Noncritical optional failures are recorded but do not halt the run.
# ------------------------------------------------------------
required_primary <- c(
  "null", "nonspatial_core", "spatial_raw_core", "spatial_plus_core"
)
primary_status <- vapply(primary_fits[required_primary], fit_ok, logical(1))
cv_counts <- table(cv_by_fold$model[is.finite(cv_by_fold$weighted_brier)])

validation <- data.frame(
  check = c(
    "required_inputs_found",
    "joined_cell_count_sufficient",
    "grid_cell_ids_unique",
    "five_latitude_bands_present",
    "four_covariate_spatial_models_completed",
    "required_primary_outcome_models_completed",
    "spatial_plus_key_effects_extracted",
    "five_cv_folds_completed_for_required_models",
    "nine_response_denominator_sensitivities_attempted",
    "three_basis_dimension_sensitivities_completed",
    "residual_moran_diagnostics_written",
    "primary_figure_written",
    "recommendation_written"
  ),
  severity = c(
    "critical", "critical", "critical", "critical", "critical",
    "critical", "critical", "critical", "noncritical", "noncritical",
    "noncritical", "critical", "critical"
  ),
  passed = c(
    all(file.exists(required_inputs)),
    nrow(d) >= 150,
    !anyDuplicated(d$grid_cell_id),
    length(unique(stats::na.omit(as.character(d$latitude_band)))) == 5,
    length(sp_full$covariate_models) == 4 &&
      all(vapply(sp_full$covariate_models, fit_ok, logical(1))),
    all(primary_status),
    nrow(key_sp) == 2,
    all(cv_counts[required_primary] == 5),
    length(unique(paste(
      sensitivity_coefficients$dataset,
      sensitivity_coefficients$response
    ))) == 9,
    length(unique(k_sensitivity_table$requested_k)) == 3,
    file.exists(file.path(out_dir, "12K_residual_moran_permutation.csv")),
    file.exists(file.path(
      fig_dir, "Figure_3_spatial_plus_environmental_effects.png"
    )),
    file.exists(file.path(out_dir, "12K_model_recommendation.txt"))
  ),
  detail = c(
    paste(required_inputs, collapse = "; "),
    paste0("joined complete cells=", nrow(d)),
    paste0("duplicate IDs=", anyDuplicated(d$grid_cell_id)),
    paste(unique(as.character(d$latitude_band)), collapse = "; "),
    paste(names(sp_full$covariate_models), collapse = "; "),
    paste(names(primary_status), primary_status, collapse = "; "),
    paste(key_sp$term, collapse = "; "),
    paste(names(cv_counts), as.integer(cv_counts), collapse = "; "),
    paste0("rows=", nrow(sensitivity_coefficients)),
    paste0("k values=", paste(unique(k_sensitivity_table$requested_k), collapse = "; ")),
    paste0("models=", nrow(moran_summary)),
    file.path(fig_dir, "Figure_3_spatial_plus_environmental_effects.png"),
    robustness_class
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(validation, file.path(out_dir, "12K_validation.csv"), row.names = FALSE)

critical_failed <- validation$severity == "critical" & !validation$passed
log_msg("Validation checks passed: ", sum(validation$passed), "/", nrow(validation))
log_msg("Spatial+ classification: ", robustness_class)

if (any(critical_failed)) {
  log_msg(
    "STEP 12K COMPLETED WITH CRITICAL VALIDATION FAILURE: ",
    paste(validation$check[critical_failed], collapse = "; ")
  )
  stop(
    "Step 12K critical validation failed: ",
    paste(validation$check[critical_failed], collapse = "; ")
  )
} else {
  if (all(validation$passed)) {
    log_msg("STEP 12K COMPLETED SUCCESSFULLY")
  } else {
    log_msg(
      "STEP 12K COMPLETED WITH NONCRITICAL WARNINGS: ",
      paste(validation$check[!validation$passed], collapse = "; ")
    )
  }
}
