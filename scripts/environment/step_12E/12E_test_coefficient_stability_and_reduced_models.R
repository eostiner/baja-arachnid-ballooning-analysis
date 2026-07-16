#!/usr/bin/env Rscript

# ============================================================
# STEP 12E — Coefficient stability, parsimonious models, and
#             prediction-readiness diagnostics
# Baja Ballooning Publication
#
# Purpose
#   1. Refit the Step 12D combined model across denominator,
#      taxonomy, and trait-confidence sensitivity datasets using
#      the SAME transformations and scaling parameters from the
#      current Step 12C primary dataset.
#   2. Quantify coefficient direction and magnitude stability.
#   3. Compare transparent parsimonious models focused on VPD and
#      wind seasonality with the null, geography, and full models.
#   4. Re-run leave-one-latitude-band-out validation by held-out band.
#   5. Check quasibinomial uncertainty and influential-cell sensitivity.
#
# No prediction map is produced in Step 12E.
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) {
  normalizePath(args[1], mustWork = TRUE)
} else {
  normalizePath("~/Desktop/Baja_Ballooning_Pipeline", mustWork = TRUE)
}

version <- "12E_v2_2026-07-16"
step12c_dir <- file.path(project_root, "04_analysis", "12C_cell_environment_model_table")
step12d_dir <- file.path(project_root, "04_analysis", "12D_predictor_screening_model_comparison")
out_dir <- file.path(project_root, "04_analysis", "12E_coefficient_stability_diagnostics")
archive_dir <- file.path(project_root, "08_archive")

if (dir.exists(out_dir)) {
  dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y%m%dT%H%M%S")
  archived <- file.path(archive_dir, paste0("12E_coefficient_stability_diagnostics_", stamp))
  ok <- file.rename(out_dir, archived)
  if (!ok) stop("Could not archive prior Step 12E output: ", out_dir)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "models"), recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(out_dir, "12E_analysis_log.txt")
log_con <- file(log_path, open = "wt")
on.exit(close(log_con), add = TRUE)
log_msg <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n")
  writeLines(txt, log_con)
  flush(log_con)
}

log_msg("STEP 12E STARTED")
log_msg("Version: ", version)
log_msg("Project root: ", project_root)

input_primary <- file.path(step12c_dir, "12C_primary_glm_candidate_table.csv")
input_ge5 <- file.path(step12c_dir, "12C_sensitivity_glm_candidate_ge5.csv")
input_ge10 <- file.path(step12c_dir, "12C_sensitivity_glm_candidate_ge10.csv")
model_rds <- file.path(step12d_dir, "models", "12D_model_objects.rds")
influence_path <- file.path(step12d_dir, "12D_combined_model_influence_audit.csv")

required_files <- c(input_primary, input_ge5, input_ge10, model_rds)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing required input(s): ", paste(missing_files, collapse = "; "))

read_candidate <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"grid_cell_id" %in% names(d)) stop("grid_cell_id missing from ", path)
  if (anyDuplicated(d$grid_cell_id)) stop("Duplicate grid_cell_id values in ", path)
  d
}

primary_raw <- read_candidate(input_primary)
ge5_raw <- read_candidate(input_ge5)
ge10_raw <- read_candidate(input_ge10)
obj12d <- readRDS(model_rds)

actual_counts <- c(primary = nrow(primary_raw), denominator_ge5 = nrow(ge5_raw), denominator_ge10 = nrow(ge10_raw))
log_msg("Candidate cells — primary: ", actual_counts[["primary"]],
        "; >=5: ", actual_counts[["denominator_ge5"]],
        "; >=10: ", actual_counts[["denominator_ge10"]])

for (candidate_table in list(primary_raw, ge5_raw, ge10_raw)) {
  missing_candidate_fields <- setdiff(c("grid_cell_id", "classified_genera_primary"), names(candidate_table))
  if (length(missing_candidate_fields)) {
    stop("A Step 12C candidate table is missing required fields: ",
         paste(missing_candidate_fields, collapse = "; "))
  }
}

candidate_sets_consistent <- isTRUE(
  nrow(primary_raw) > 0L &&
  nrow(ge5_raw) > 0L &&
  nrow(ge10_raw) > 0L &&
  all(ge5_raw$grid_cell_id %in% primary_raw$grid_cell_id) &&
  all(ge10_raw$grid_cell_id %in% ge5_raw$grid_cell_id) &&
  all(ge5_raw$classified_genera_primary >= 5L) &&
  all(ge10_raw$classified_genera_primary >= 10L)
)

band_levels <- c("23-24N", "24-26N", "26-28N", "28-30N", "30-32N")

# ----------------------------
# Predictor transformations and fixed primary scaling
# ----------------------------
derive_predictors <- function(d) {
  d$latitude_band <- factor(d$latitude_band, levels = band_levels)
  d$log_precip <- log1p(pmax(d$precip_annual_mean_mm, 0))
  d$log_relief <- log1p(pmax(d$relief_5km_m, 0))
  d$asin_barren <- asin(sqrt(pmin(pmax(d$lc_barren_sparse_prop, 0), 1)))
  d$log_distance_water <- log1p(pmax(d$distance_to_modis_water_km, 0))
  d
}

primary_scaler <- obj12d$scaler
if (is.null(primary_scaler$env) || is.null(primary_scaler$lat_center)) {
  stop("Step 12D model object does not contain the expected primary scaler")
}

apply_fixed_primary_scaler <- function(d) {
  d <- derive_predictors(d)
  d$z_lat <- (d$centroid_latitude - primary_scaler$lat_center) / primary_scaler$lat_scale
  d$z_lat2 <- (d$z_lat^2 - primary_scaler$lat2_center) / primary_scaler$lat2_scale
  for (i in seq_len(nrow(primary_scaler$env))) {
    v <- primary_scaler$env$variable[i]
    if (!v %in% names(d)) stop("Transformed predictor missing: ", v)
    d[[paste0("z_", v)]] <- (d[[v]] - primary_scaler$env$center[i]) / primary_scaler$env$scale[i]
  }
  d
}

fit_training_scaler <- function(d, env_vars) {
  d <- derive_predictors(d)
  lat_center <- mean(d$centroid_latitude, na.rm = TRUE)
  lat_scale <- sd(d$centroid_latitude, na.rm = TRUE)
  if (!is.finite(lat_scale) || lat_scale == 0) lat_scale <- 1
  d$z_lat <- (d$centroid_latitude - lat_center) / lat_scale
  lat2_center <- mean(d$z_lat^2, na.rm = TRUE)
  lat2_scale <- sd(d$z_lat^2, na.rm = TRUE)
  if (!is.finite(lat2_scale) || lat2_scale == 0) lat2_scale <- 1
  env_center <- sapply(env_vars, function(v) mean(d[[v]], na.rm = TRUE))
  env_scale <- sapply(env_vars, function(v) sd(d[[v]], na.rm = TRUE))
  env_scale[!is.finite(env_scale) | env_scale == 0] <- 1
  list(
    lat_center = lat_center, lat_scale = lat_scale,
    lat2_center = lat2_center, lat2_scale = lat2_scale,
    env = data.frame(variable = env_vars, center = as.numeric(env_center),
                     scale = as.numeric(env_scale), stringsAsFactors = FALSE)
  )
}

apply_any_scaler <- function(d, scaler) {
  d <- derive_predictors(d)
  d$z_lat <- (d$centroid_latitude - scaler$lat_center) / scaler$lat_scale
  d$z_lat2 <- (d$z_lat^2 - scaler$lat2_center) / scaler$lat2_scale
  for (i in seq_len(nrow(scaler$env))) {
    v <- scaler$env$variable[i]
    d[[paste0("z_", v)]] <- (d[[v]] - scaler$env$center[i]) / scaler$env$scale[i]
  }
  d
}

full_combined_terms <- obj12d$model_terms$combined
full_environment_terms <- obj12d$model_terms$environment
if (is.null(full_combined_terms) || is.null(full_environment_terms)) {
  stop("Step 12D model terms missing from model object")
}

key_terms <- c("z_vpd_mean_kpa", "z_wind_monthly_sd_ms")
missing_key <- setdiff(key_terms, full_combined_terms)
if (length(missing_key)) stop("Expected key terms not retained in Step 12D: ", paste(missing_key, collapse = "; "))

env_raw_union <- unique(primary_scaler$env$variable)

response_specs <- data.frame(
  response = c("primary", "taxonomy_strict", "low_conf_exclusion"),
  success = c("ballooning_genera_primary", "ballooning_genera_taxonomy_strict", "ballooning_genera_low_conf_exclusion"),
  failure = c("non_ballooning_genera_primary", "non_ballooning_genera_taxonomy_strict", "non_ballooning_genera_low_conf_exclusion"),
  denominator = c("classified_genera_primary", "classified_genera_taxonomy_strict", "classified_genera_low_conf_exclusion"),
  stringsAsFactors = FALSE
)

dataset_specs <- list(
  primary = primary_raw,
  denominator_ge5 = ge5_raw,
  denominator_ge10 = ge10_raw
)

dataset_labels <- c(
  primary = "Primary cells",
  denominator_ge5 = "Cells with >=5 genera",
  denominator_ge10 = "Cells with >=10 genera"
)
response_labels <- c(
  primary = "Primary traits",
  taxonomy_strict = "Taxonomy-strict",
  low_conf_exclusion = "LOW-confidence excluded"
)

response_formula <- function(success, failure, terms) {
  rhs <- if (length(terms)) paste(terms, collapse = " + ") else "1"
  as.formula(paste0("cbind(", success, ", ", failure, ") ~ ", rhs))
}

fit_safe <- function(formula, data, family = binomial(link = "logit")) {
  tryCatch(glm(formula, data = data, family = family),
           error = function(e) structure(list(error = conditionMessage(e)), class = "fit_error"))
}
is_fit_error <- function(x) inherits(x, "fit_error")

model_converged <- function(fit) !is_fit_error(fit) && isTRUE(fit$converged)
pearson_dispersion <- function(fit) {
  if (is_fit_error(fit)) return(NA_real_)
  r <- residuals(fit, type = "pearson")
  df <- df.residual(fit)
  if (!is.finite(df) || df <= 0) return(NA_real_)
  sum(r^2, na.rm = TRUE) / df
}
model_aicc <- function(fit) {
  if (is_fit_error(fit)) return(NA_real_)
  a <- tryCatch(AIC(fit), error = function(e) NA_real_)
  ll <- tryCatch(logLik(fit), error = function(e) NULL)
  n <- tryCatch(nobs(fit), error = function(e) NA_integer_)
  if (!is.finite(a) || is.null(ll) || !is.finite(n)) return(NA_real_)
  k <- as.numeric(attr(ll, "df"))
  if (!is.finite(k) || n <= k + 1) return(NA_real_)
  a + (2 * k * (k + 1)) / (n - k - 1)
}

extract_coef <- function(fit, dataset, response, model) {
  if (is_fit_error(fit)) return(data.frame())
  sm <- summary(fit)$coefficients
  out <- data.frame(
    dataset = dataset,
    response = response,
    model = model,
    term = rownames(sm),
    estimate = sm[, "Estimate"],
    std_error = sm[, "Std. Error"],
    statistic = sm[, if ("z value" %in% colnames(sm)) "z value" else "t value"],
    p_value = sm[, if ("Pr(>|z|)" %in% colnames(sm)) "Pr(>|z|)" else "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
  out$odds_ratio <- exp(out$estimate)
  out$ci_low <- exp(out$estimate - 1.96 * out$std_error)
  out$ci_high <- exp(out$estimate + 1.96 * out$std_error)
  out
}

# ----------------------------
# Full combined-model coefficient stability using FIXED primary scaling
# ----------------------------
stability_fits <- list()
stability_summary_rows <- list()
stability_coef_rows <- list()

for (ds_name in names(dataset_specs)) {
  d_raw <- dataset_specs[[ds_name]]
  for (i in seq_len(nrow(response_specs))) {
    rs <- response_specs[i, ]
    keep <- is.finite(d_raw[[rs$denominator]]) & d_raw[[rs$denominator]] > 0
    d0 <- d_raw[keep, , drop = FALSE]
    d <- apply_fixed_primary_scaler(d0)
    complete <- complete.cases(d[, full_combined_terms, drop = FALSE])
    d <- d[complete, , drop = FALSE]
    fm <- response_formula(rs$success, rs$failure, full_combined_terms)
    fit <- fit_safe(fm, d, binomial(link = "logit"))
    key <- paste(ds_name, rs$response, sep = "__")
    stability_fits[[key]] <- fit
    stability_summary_rows[[length(stability_summary_rows) + 1L]] <- data.frame(
      dataset = ds_name,
      dataset_label = unname(dataset_labels[ds_name]),
      response = rs$response,
      response_label = unname(response_labels[rs$response]),
      n_cells = nrow(d),
      converged = model_converged(fit),
      AICc = model_aicc(fit),
      pearson_dispersion = pearson_dispersion(fit),
      error = if (is_fit_error(fit)) fit$error else "",
      stringsAsFactors = FALSE
    )
    cc <- extract_coef(fit, ds_name, rs$response, "full_combined_fixed_scaling")
    if (nrow(cc)) stability_coef_rows[[length(stability_coef_rows) + 1L]] <- cc
  }
}

stability_model_summary <- do.call(rbind, stability_summary_rows)
stability_coef <- do.call(rbind, stability_coef_rows)
write.csv(stability_model_summary, file.path(out_dir, "12E_fixed_scaling_model_summary.csv"), row.names = FALSE)
write.csv(stability_coef, file.path(out_dir, "12E_fixed_scaling_coefficients.csv"), row.names = FALSE)

reference_coef <- stability_coef[
  stability_coef$dataset == "primary" & stability_coef$response == "primary",
  c("term", "estimate", "odds_ratio"), drop = FALSE
]
names(reference_coef)[2:3] <- c("reference_estimate", "reference_odds_ratio")
stability_coef <- merge(stability_coef, reference_coef, by = "term", all.x = TRUE, sort = FALSE)
stability_coef$estimate_difference_from_reference <- stability_coef$estimate - stability_coef$reference_estimate
stability_coef$relative_estimate_change_pct <- ifelse(
  is.finite(stability_coef$reference_estimate) & abs(stability_coef$reference_estimate) > 1e-8,
  100 * (stability_coef$estimate - stability_coef$reference_estimate) / abs(stability_coef$reference_estimate),
  NA_real_
)
stability_coef$same_direction_as_reference <- sign(stability_coef$estimate) == sign(stability_coef$reference_estimate)
write.csv(stability_coef, file.path(out_dir, "12E_coefficient_stability_all_fits.csv"), row.names = FALSE)

non_intercept <- stability_coef[stability_coef$term != "(Intercept)", , drop = FALSE]
stability_by_term <- do.call(rbind, lapply(split(non_intercept, non_intercept$term), function(d) {
  ref <- unique(d$reference_estimate[is.finite(d$reference_estimate)])
  ref <- if (length(ref)) ref[1] else NA_real_
  data.frame(
    term = unique(d$term),
    n_fits = nrow(d),
    n_converged_coefficients = sum(is.finite(d$estimate)),
    reference_estimate = ref,
    all_same_direction_as_reference = all(d$same_direction_as_reference[is.finite(d$estimate)], na.rm = TRUE),
    estimate_min = min(d$estimate, na.rm = TRUE),
    estimate_max = max(d$estimate, na.rm = TRUE),
    odds_ratio_min = min(d$odds_ratio, na.rm = TRUE),
    odds_ratio_max = max(d$odds_ratio, na.rm = TRUE),
    n_ci_excluding_one = sum(d$ci_high < 1 | d$ci_low > 1, na.rm = TRUE),
    n_p_below_0_05 = sum(d$p_value < 0.05, na.rm = TRUE),
    max_absolute_relative_change_pct = if (all(is.na(d$relative_estimate_change_pct))) NA_real_ else max(abs(d$relative_estimate_change_pct), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write.csv(stability_by_term, file.path(out_dir, "12E_coefficient_stability_by_term.csv"), row.names = FALSE)

key_stability <- stability_by_term[stability_by_term$term %in% key_terms, , drop = FALSE]
write.csv(key_stability, file.path(out_dir, "12E_key_effect_stability_summary.csv"), row.names = FALSE)

# ----------------------------
# Parsimonious primary models
# ----------------------------
primary_scaled <- apply_fixed_primary_scaler(primary_raw)

candidate_model_terms <- list(
  null = character(),
  geography = c("z_lat", "z_lat2"),
  vpd_only = "z_vpd_mean_kpa",
  wind_seasonality_only = "z_wind_monthly_sd_ms",
  reduced_environment = key_terms,
  reduced_combined = c("z_lat", "z_lat2", key_terms),
  full_environment = full_environment_terms,
  full_combined = full_combined_terms
)

primary_models <- list()
comparison_rows <- list()
coef_rows <- list()
for (nm in names(candidate_model_terms)) {
  fm <- response_formula("ballooning_genera_primary", "non_ballooning_genera_primary", candidate_model_terms[[nm]])
  fit <- fit_safe(fm, primary_scaled, binomial(link = "logit"))
  primary_models[[nm]] <- fit
  comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
    model = nm,
    n_cells = if (is_fit_error(fit)) NA_integer_ else nobs(fit),
    n_terms = length(candidate_model_terms[[nm]]),
    converged = model_converged(fit),
    logLik = if (is_fit_error(fit)) NA_real_ else as.numeric(logLik(fit)),
    AIC = if (is_fit_error(fit)) NA_real_ else AIC(fit),
    AICc = model_aicc(fit),
    pearson_dispersion = pearson_dispersion(fit),
    error = if (is_fit_error(fit)) fit$error else "",
    stringsAsFactors = FALSE
  )
  cc <- extract_coef(fit, "primary", "primary", nm)
  if (nrow(cc)) coef_rows[[length(coef_rows) + 1L]] <- cc
}
primary_comparison <- do.call(rbind, comparison_rows)
if (any(is.finite(primary_comparison$AICc))) {
  primary_comparison$delta_AICc <- primary_comparison$AICc - min(primary_comparison$AICc, na.rm = TRUE)
} else {
  primary_comparison$delta_AICc <- NA_real_
}
primary_coef <- do.call(rbind, coef_rows)
write.csv(primary_comparison, file.path(out_dir, "12E_primary_parsimonious_model_comparison.csv"), row.names = FALSE)
write.csv(primary_coef, file.path(out_dir, "12E_primary_parsimonious_model_coefficients.csv"), row.names = FALSE)

# ----------------------------
# Leave-one-latitude-band-out CV for parsimonious and benchmark models
# ----------------------------
clip_prob <- function(p, eps = 1e-8) pmin(pmax(p, eps), 1 - eps)

cv_rows <- list()
cv_prediction_rows <- list()
for (model_name in names(candidate_model_terms)) {
  terms <- candidate_model_terms[[model_name]]
  for (b in band_levels) {
    train_raw <- primary_raw[as.character(primary_raw$latitude_band) != b, , drop = FALSE]
    test_raw <- primary_raw[as.character(primary_raw$latitude_band) == b, , drop = FALSE]
    if (!nrow(train_raw) || !nrow(test_raw)) next
    scaler_fold <- fit_training_scaler(train_raw, env_raw_union)
    train <- apply_any_scaler(train_raw, scaler_fold)
    test <- apply_any_scaler(test_raw, scaler_fold)
    fm <- response_formula("ballooning_genera_primary", "non_ballooning_genera_primary", terms)
    fit <- fit_safe(fm, train, binomial(link = "logit"))
    if (is_fit_error(fit)) {
      cv_rows[[length(cv_rows) + 1L]] <- data.frame(
        model = model_name, held_out_band = b, n_cells = nrow(test),
        total_classified_genera = sum(test$classified_genera_primary),
        weighted_brier = NA_real_, weighted_mae = NA_real_,
        binomial_log_score_per_genus = NA_real_, error = fit$error,
        stringsAsFactors = FALSE
      )
      next
    }
    p <- clip_prob(as.numeric(predict(fit, newdata = test, type = "response")))
    y <- test$ballooning_genera_primary
    n <- test$classified_genera_primary
    obs <- y / n
    valid <- is.finite(p) & is.finite(obs) & is.finite(n) & n > 0
    w <- n[valid]
    cv_rows[[length(cv_rows) + 1L]] <- data.frame(
      model = model_name,
      held_out_band = b,
      n_cells = sum(valid),
      total_classified_genera = sum(w),
      weighted_brier = sum(w * (p[valid] - obs[valid])^2) / sum(w),
      weighted_mae = sum(w * abs(p[valid] - obs[valid])) / sum(w),
      binomial_log_score_per_genus = -sum(dbinom(y[valid], size = n[valid], prob = p[valid], log = TRUE)) / sum(w),
      error = "",
      stringsAsFactors = FALSE
    )
    cv_prediction_rows[[length(cv_prediction_rows) + 1L]] <- data.frame(
      grid_cell_id = test$grid_cell_id,
      latitude_band = as.character(test$latitude_band),
      model = model_name,
      observed_ballooning = y,
      denominator = n,
      observed_proportion = obs,
      predicted_proportion = p,
      stringsAsFactors = FALSE
    )
  }
}

cv_by_fold <- do.call(rbind, cv_rows)
cv_predictions <- do.call(rbind, cv_prediction_rows)
cv_summary <- do.call(rbind, lapply(split(cv_by_fold, cv_by_fold$model), function(d) {
  valid <- is.finite(d$weighted_brier) & is.finite(d$total_classified_genera)
  w <- d$total_classified_genera[valid]
  data.frame(
    model = unique(d$model),
    folds_completed = sum(valid),
    total_classified_genera = sum(w),
    weighted_brier = weighted.mean(d$weighted_brier[valid], w),
    weighted_mae = weighted.mean(d$weighted_mae[valid], w),
    binomial_log_score_per_genus = weighted.mean(d$binomial_log_score_per_genus[valid], w),
    stringsAsFactors = FALSE
  )
}))
cv_summary <- cv_summary[match(names(candidate_model_terms), cv_summary$model), , drop = FALSE]
null_row <- cv_summary[cv_summary$model == "null", , drop = FALSE]
cv_summary$delta_log_score_vs_null <- cv_summary$binomial_log_score_per_genus - null_row$binomial_log_score_per_genus
cv_summary$delta_brier_vs_null <- cv_summary$weighted_brier - null_row$weighted_brier
write.csv(cv_by_fold, file.path(out_dir, "12E_leave_one_band_out_cv_by_fold.csv"), row.names = FALSE)
write.csv(cv_predictions, file.path(out_dir, "12E_leave_one_band_out_cv_predictions.csv"), row.names = FALSE)
write.csv(cv_summary, file.path(out_dir, "12E_leave_one_band_out_cv_summary.csv"), row.names = FALSE)

# ----------------------------
# Quasibinomial uncertainty check on primary models
# ----------------------------
quasi_models <- c("reduced_environment", "reduced_combined", "full_environment", "full_combined")
quasi_rows <- list()
for (nm in quasi_models) {
  fm <- response_formula("ballooning_genera_primary", "non_ballooning_genera_primary", candidate_model_terms[[nm]])
  fit <- fit_safe(fm, primary_scaled, quasibinomial(link = "logit"))
  cc <- extract_coef(fit, "primary", "primary", paste0(nm, "_quasibinomial"))
  if (nrow(cc)) {
    cc$pearson_dispersion <- pearson_dispersion(fit)
    quasi_rows[[length(quasi_rows) + 1L]] <- cc
  }
}
quasi_coef <- do.call(rbind, quasi_rows)
write.csv(quasi_coef, file.path(out_dir, "12E_quasibinomial_robustness_coefficients.csv"), row.names = FALSE)

# ----------------------------
# Influential-cell sensitivity for key effects
# ----------------------------
influence_leave_one_out <- data.frame()
influence_all_flagged <- data.frame()
if (file.exists(influence_path)) {
  infl <- read.csv(influence_path, stringsAsFactors = FALSE, check.names = FALSE)
  cook_col <- "cooks_distance_binomial_companion"
  if (all(c("grid_cell_id", cook_col) %in% names(infl))) {
    infl <- infl[order(infl[[cook_col]], decreasing = TRUE), , drop = FALSE]
    top_ids <- head(infl$grid_cell_id[is.finite(infl[[cook_col]])], 10)
    loo_rows <- list()
    for (id in top_ids) {
      d <- primary_scaled[primary_scaled$grid_cell_id != id, , drop = FALSE]
      fm <- response_formula("ballooning_genera_primary", "non_ballooning_genera_primary", full_combined_terms)
      fit <- fit_safe(fm, d, binomial(link = "logit"))
      cc <- extract_coef(fit, paste0("exclude_", id), "primary", "full_combined_leave_one_influential_out")
      cc <- cc[cc$term %in% key_terms, , drop = FALSE]
      if (nrow(cc)) {
        cc$excluded_grid_cell_id <- id
        cc$excluded_cooks_distance <- infl[[cook_col]][match(id, infl$grid_cell_id)]
        loo_rows[[length(loo_rows) + 1L]] <- cc
      }
    }
    if (length(loo_rows)) influence_leave_one_out <- do.call(rbind, loo_rows)

    if ("flag_cook" %in% names(infl)) {
      flag <- as.logical(infl$flag_cook)
      flag[is.na(flag)] <- FALSE
      flagged_ids <- infl$grid_cell_id[flag]
      if (length(flagged_ids) > 0 && length(flagged_ids) <= 30) {
        d <- primary_scaled[!primary_scaled$grid_cell_id %in% flagged_ids, , drop = FALSE]
        fm <- response_formula("ballooning_genera_primary", "non_ballooning_genera_primary", full_combined_terms)
        fit <- fit_safe(fm, d, binomial(link = "logit"))
        influence_all_flagged <- extract_coef(fit, "exclude_all_cook_flags", "primary", "full_combined_exclude_all_cook_flags")
        influence_all_flagged$n_excluded <- length(flagged_ids)
      }
    }
  }
}
write.csv(influence_leave_one_out, file.path(out_dir, "12E_influential_cell_leave_one_out_key_effects.csv"), row.names = FALSE)
write.csv(influence_all_flagged, file.path(out_dir, "12E_exclude_all_cook_flags_coefficients.csv"), row.names = FALSE)

# ----------------------------
# Figures
# ----------------------------
fig_dir <- file.path(out_dir, "figures")
save_plot <- function(stem, width, height, plot_fun) {
  png(file.path(fig_dir, paste0(stem, ".png")), width = width, height = height, units = "in", res = 400)
  plot_fun()
  dev.off()
  pdf(file.path(fig_dir, paste0(stem, ".pdf")), width = width, height = height, useDingbats = FALSE)
  plot_fun()
  dev.off()
}

term_labels <- c(
  z_lat = "Latitude",
  z_lat2 = "Latitude squared",
  z_vpd_mean_kpa = "Vapor-pressure deficit",
  z_log_relief = "Topographic relief",
  z_asin_barren = "Barren/sparse cover",
  z_upward_sensible_heat_mean_wm2 = "Upward sensible heat",
  z_log_distance_water = "Distance to mapped water",
  z_wind_monthly_sd_ms = "Wind seasonality"
)

# Key-effect stability forest plot across all nine sensitivity fits.
plot_key <- stability_coef[stability_coef$term %in% key_terms, , drop = FALSE]
plot_key$fit_label <- paste(unname(dataset_labels[plot_key$dataset]),
                            unname(response_labels[plot_key$response]), sep = " — ")
plot_key$term_label <- unname(term_labels[plot_key$term])
plot_key <- plot_key[order(plot_key$term_label, plot_key$dataset, plot_key$response), , drop = FALSE]

if (nrow(plot_key)) {
  save_plot("Figure_12E_key_effect_stability", 9.5, 8.5, function() {
    op <- par(mar = c(5, 15, 3, 2))
    on.exit(par(op))
    gap <- 2
    terms_unique <- unique(plot_key$term_label)
    y <- numeric(nrow(plot_key))
    pos <- 1
    for (tt in terms_unique) {
      idx <- which(plot_key$term_label == tt)
      y[idx] <- seq(pos, length.out = length(idx))
      pos <- max(y[idx]) + gap
    }
    x_rng <- range(c(plot_key$ci_low, plot_key$ci_high), finite = TRUE)
    plot(plot_key$odds_ratio, y, log = "x", xlim = x_rng,
         ylim = c(0.5, max(y) + 0.8), yaxt = "n", pch = 19,
         xlab = "Odds ratio per 1 SD increase (log scale)", ylab = "",
         main = "Sensitivity of key environmental effects")
    segments(plot_key$ci_low, y, plot_key$ci_high, y)
    axis(2, at = y, labels = plot_key$fit_label, las = 2, cex.axis = 0.68)
    abline(v = 1, lty = 2)
    for (tt in terms_unique) {
      idx <- which(plot_key$term_label == tt)
      text(x_rng[1], max(y[idx]) + 0.75, tt, adj = 0, font = 2, cex = 0.85)
    }
  })
}

# CV comparison plot.
if (nrow(cv_summary)) {
  save_plot("Figure_12E_spatial_cv_model_comparison", 8.5, 6.5, function() {
    d <- cv_summary[order(cv_summary$binomial_log_score_per_genus), , drop = FALSE]
    op <- par(mar = c(5, 10, 3, 1))
    on.exit(par(op))
    barplot(d$binomial_log_score_per_genus, names.arg = d$model, horiz = TRUE, las = 1,
            xlab = "Leave-one-band-out log score per classified genus (lower is better)",
            main = "Spatial transferability of candidate models")
  })
}

# Coefficient stability heatmap for all full-combined terms.
heat <- non_intercept
if (nrow(heat)) {
  heat$fit_label <- paste(heat$dataset, heat$response, sep = "__")
  terms <- unique(heat$term)
  fits <- unique(heat$fit_label)
  mat <- matrix(NA_real_, nrow = length(terms), ncol = length(fits), dimnames = list(terms, fits))
  for (i in seq_len(nrow(heat))) mat[heat$term[i], heat$fit_label[i]] <- heat$estimate[i]
  max_abs <- max(abs(mat), na.rm = TRUE)
  save_plot("Figure_12E_all_coefficient_stability_heatmap", 11, 7.5, function() {
    op <- par(mar = c(12, 12, 3, 2))
    on.exit(par(op))
    image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat), zlim = c(-max_abs, max_abs),
          axes = FALSE, xlab = "", ylab = "", main = "Standardized coefficient stability")
    axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.62)
    axis(2, at = seq_len(nrow(mat)), labels = rownames(mat), las = 2, cex.axis = 0.7)
    box()
  })
}

# ----------------------------
# Prediction-readiness recommendation
# ----------------------------
null_cv <- cv_summary[cv_summary$model == "null", , drop = FALSE]
reduced_env_cv <- cv_summary[cv_summary$model == "reduced_environment", , drop = FALSE]
reduced_combined_cv <- cv_summary[cv_summary$model == "reduced_combined", , drop = FALSE]

best_model <- cv_summary$model[which.min(cv_summary$binomial_log_score_per_genus)]
best_score <- min(cv_summary$binomial_log_score_per_genus, na.rm = TRUE)
key_direction_stable <- nrow(key_stability) == length(key_terms) && all(key_stability$all_same_direction_as_reference)
reduced_beats_null <- FALSE
if (nrow(null_cv) && nrow(reduced_env_cv) && nrow(reduced_combined_cv)) {
  reduced_beats_null <- any(c(
    reduced_env_cv$binomial_log_score_per_genus < null_cv$binomial_log_score_per_genus &&
      reduced_env_cv$weighted_brier < null_cv$weighted_brier,
    reduced_combined_cv$binomial_log_score_per_genus < null_cv$binomial_log_score_per_genus &&
      reduced_combined_cv$weighted_brier < null_cv$weighted_brier
  ))
}

prediction_ready <- key_direction_stable && reduced_beats_null
recommendation <- c(
  "STEP 12E MODEL DIAGNOSTIC INTERPRETATION",
  "========================================",
  "",
  paste0("Key effects directionally stable across all fitted sensitivity combinations: ", key_direction_stable),
  paste0("A reduced environmental model beats the null on both leave-one-band-out log score and Brier score: ", reduced_beats_null),
  paste0("Best leave-one-band-out model by log score: ", best_model, " (", signif(best_score, 5), ")"),
  paste0("Prediction-readiness criterion met: ", prediction_ready),
  "",
  if (prediction_ready) {
    "The reduced model may proceed to Spatial+ deconfounding. This criterion does not authorize an extrapolative habitat-suitability or prediction heat map."
  } else {
    "Do not present an extrapolative heat map. Continue only with clearly labeled within-sample diagnostics and spatial deconfounding."
  },
  "",
  "Interpretation rules:",
  "  * Directional stability is more important than repeated p < 0.05 in the smaller sensitivity subsets.",
  paste0("  * The >=10-genus subset has ", actual_counts[["denominator_ge10"]], " cells and should be treated as a conservative sensitivity check."),
  paste0("  * All coefficient magnitudes are comparable because every sensitivity fit uses transformations and scaling parameters from the current ", actual_counts[["primary"]], "-cell primary dataset."),
  "  * Leave-one-band-out validation remains the primary test of spatial transferability."
)
writeLines(recommendation, file.path(out_dir, "12E_model_recommendation.txt"))

# Save model objects for possible Step 12F use.
saveRDS(list(
  version = version,
  primary_scaler = primary_scaler,
  full_combined_terms = full_combined_terms,
  full_environment_terms = full_environment_terms,
  key_terms = key_terms,
  stability_fits = stability_fits,
  primary_candidate_models = primary_models,
  candidate_model_terms = candidate_model_terms,
  cv_summary = cv_summary,
  prediction_ready = prediction_ready
), file.path(out_dir, "models", "12E_model_objects.rds"))

# ----------------------------
# Validation and provenance
# ----------------------------
validation <- data.frame(
  check = c(
    "candidate_sets_consistent_with_step12C",
    "step12D_primary_scaler_loaded",
    "nine_fixed_scaling_sensitivity_models_attempted",
    "all_fixed_scaling_models_converged",
    "key_terms_present",
    "all_candidate_primary_models_converged",
    "five_cv_folds_per_candidate_model",
    "quasibinomial_robustness_written",
    "coefficient_stability_outputs_written",
    "prediction_readiness_decision_written"
  ),
  passed = c(
    candidate_sets_consistent,
    !is.null(primary_scaler$env),
    nrow(stability_model_summary) == 9,
    all(stability_model_summary$converged),
    all(key_terms %in% stability_coef$term),
    all(vapply(primary_models, model_converged, logical(1))),
    all(table(cv_by_fold$model) == 5),
    file.exists(file.path(out_dir, "12E_quasibinomial_robustness_coefficients.csv")),
    file.exists(file.path(out_dir, "12E_key_effect_stability_summary.csv")),
    file.exists(file.path(out_dir, "12E_model_recommendation.txt"))
  ),
  detail = c(
    paste(names(actual_counts), actual_counts, collapse = "; "),
    paste0("environment scaling rows=", nrow(primary_scaler$env)),
    paste0("models=", nrow(stability_model_summary)),
    paste(stability_model_summary$dataset, stability_model_summary$response, stability_model_summary$converged, collapse = "; "),
    paste(key_terms, collapse = "; "),
    paste(names(primary_models), vapply(primary_models, model_converged, logical(1)), collapse = "; "),
    paste(names(table(cv_by_fold$model)), as.integer(table(cv_by_fold$model)), collapse = "; "),
    paste0("rows=", nrow(quasi_coef)),
    paste0("key rows=", nrow(key_stability)),
    paste0("prediction_ready=", prediction_ready)
  ),
  stringsAsFactors = FALSE
)
write.csv(validation, file.path(out_dir, "12E_validation.csv"), row.names = FALSE)

provenance <- c(
  paste0("Step version: ", version),
  paste0("Primary candidate table: ", input_primary),
  paste0("GE5 candidate table: ", input_ge5),
  paste0("GE10 candidate table: ", input_ge10),
  paste0("Step 12D model object: ", model_rds),
  paste0("Full combined terms: ", paste(full_combined_terms, collapse = "; ")),
  paste0("Key diagnostic terms: ", paste(key_terms, collapse = "; ")),
  paste0("Sensitivity coefficient comparisons use the fixed transformations and scaling parameters from the current ", actual_counts[["primary"]], "-cell Step 12C primary dataset."),
  "Cross-validation transformations and scaling are estimated within each training fold to avoid leakage."
)
writeLines(provenance, file.path(out_dir, "12E_provenance.txt"))

sink(file.path(out_dir, "12E_session_info.txt"))
print(sessionInfo())
sink()

readme <- c(
  "STEP 12E OUTPUTS",
  "================",
  "",
  "This step evaluates whether Step 12D coefficient directions and magnitudes are stable across",
  "denominator, taxonomy, and trait-confidence sensitivity analyses using one fixed primary scaler.",
  "It also tests transparent reduced models and repeats leave-one-latitude-band-out validation.",
  "",
  "Review first:",
  "  12E_key_effect_stability_summary.csv",
  "  12E_coefficient_stability_all_fits.csv",
  "  12E_primary_parsimonious_model_comparison.csv",
  "  12E_leave_one_band_out_cv_summary.csv",
  "  12E_leave_one_band_out_cv_by_fold.csv",
  "  12E_quasibinomial_robustness_coefficients.csv",
  "  12E_model_recommendation.txt",
  "  figures/Figure_12E_key_effect_stability.png",
  "",
  "No prediction map is produced by Step 12E."
)
writeLines(readme, file.path(out_dir, "README_12E_OUTPUTS.txt"))

log_msg("Fixed-scaling sensitivity models: ", nrow(stability_model_summary))
log_msg("Key effects directionally stable: ", key_direction_stable)
log_msg("Best leave-one-band-out model: ", best_model)
log_msg("Reduced model beats null on both primary CV metrics: ", reduced_beats_null)
log_msg("Prediction-readiness criterion met: ", prediction_ready)
log_msg("Validation checks passed: ", sum(validation$passed), "/", nrow(validation))
if (all(validation$passed)) {
  log_msg("STEP 12E COMPLETED SUCCESSFULLY")
} else {
  log_msg("STEP 12E COMPLETED WITH FLAGS — inspect 12E_validation.csv")
}
