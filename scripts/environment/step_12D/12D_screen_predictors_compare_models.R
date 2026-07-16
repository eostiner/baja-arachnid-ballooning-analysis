#!/usr/bin/env Rscript

# ============================================================
# STEP 12D — Predictor screening and environmental model comparison
# Baja Ballooning Publication
#
# Primary response:
#   cbind(ballooning_genera_primary, non_ballooning_genera_primary)
#
# Candidate model families:
#   null
#   geography only
#   environment only
#   geography + environment
#
# Design principles:
#   * one row per occupied 25-km cell
#   * explicit use of the revised Step 12C candidate sets
#   * a priori ecological predictor pool
#   * response-independent correlation/VIF screening
#   * binomial overdispersion audit
#   * beta-binomial fit when glmmTMB is available and warranted
#   * quasibinomial fallback when beta-binomial is unavailable
#   * leave-one-latitude-band-out spatial cross-validation
#   * primary, denominator, taxonomy, and trait-confidence sensitivities
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) {
  normalizePath(args[1], mustWork = TRUE)
} else {
  normalizePath("~/Desktop/Baja_Ballooning_Pipeline", mustWork = TRUE)
}

version <- "12D_v2_2026-07-16"
step12c_dir <- file.path(project_root, "04_analysis", "12C_cell_environment_model_table")
input_primary <- file.path(step12c_dir, "12C_primary_glm_candidate_table.csv")
input_ge5 <- file.path(step12c_dir, "12C_sensitivity_glm_candidate_ge5.csv")
input_ge10 <- file.path(step12c_dir, "12C_sensitivity_glm_candidate_ge10.csv")

out_dir <- file.path(project_root, "04_analysis", "12D_predictor_screening_model_comparison")
archive_dir <- file.path(project_root, "08_archive")
if (dir.exists(out_dir)) {
  dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y%m%dT%H%M%S")
  archived <- file.path(archive_dir, paste0("12D_predictor_screening_model_comparison_", stamp))
  ok <- file.rename(out_dir, archived)
  if (!ok) stop("Could not archive prior Step 12D output: ", out_dir)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "models"), recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(out_dir, "12D_analysis_log.txt")
log_con <- file(log_path, open = "wt")
on.exit(close(log_con), add = TRUE)
log_msg <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n")
  writeLines(txt, log_con)
  flush(log_con)
}

log_msg("STEP 12D STARTED")
log_msg("Version: ", version)
log_msg("Project root: ", project_root)
log_msg("Primary input: ", input_primary)
log_msg("glmmTMB available: ", requireNamespace("glmmTMB", quietly = TRUE))

for (f in c(input_primary, input_ge5, input_ge10)) {
  if (!file.exists(f)) stop("Missing Step 12C candidate table: ", f)
}

read_candidate <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (anyDuplicated(x$grid_cell_id)) stop("Duplicate grid_cell_id values in: ", path)
  x
}

primary_raw <- read_candidate(input_primary)
ge5_raw <- read_candidate(input_ge5)
ge10_raw <- read_candidate(input_ge10)

actual_counts <- c(primary = nrow(primary_raw), ge5 = nrow(ge5_raw), ge10 = nrow(ge10_raw))
log_msg("Candidate cells — primary: ", actual_counts[["primary"]],
        "; >=5: ", actual_counts[["ge5"]],
        "; >=10: ", actual_counts[["ge10"]])

required_fields <- c(
  "grid_cell_id", "centroid_latitude", "centroid_longitude", "latitude_band",
  "centroid_x_m", "centroid_y_m",
  "ballooning_genera_primary", "non_ballooning_genera_primary",
  "classified_genera_primary",
  "ballooning_genera_taxonomy_strict", "non_ballooning_genera_taxonomy_strict",
  "classified_genera_taxonomy_strict",
  "ballooning_genera_low_conf_exclusion", "non_ballooning_genera_low_conf_exclusion",
  "classified_genera_low_conf_exclusion",
  "tmean_c", "precip_annual_mean_mm", "wind_speed_mean_ms",
  "wind_monthly_sd_ms", "vpd_mean_kpa", "soil_water_mean_frac",
  "upward_sensible_heat_mean_wm2", "evi_mean", "relief_5km_m",
  "lc_barren_sparse_prop", "distance_to_modis_water_km"
)
missing_required <- setdiff(required_fields, names(primary_raw))
if (length(missing_required)) {
  stop("Primary candidate table is missing required fields: ",
       paste(missing_required, collapse = "; "))
}

for (candidate_table in list(ge5_raw, ge10_raw)) {
  missing_candidate_fields <- setdiff(c("grid_cell_id", "classified_genera_primary"), names(candidate_table))
  if (length(missing_candidate_fields)) {
    stop("A Step 12C sensitivity table is missing required fields: ",
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
# Predictor transformation
# ----------------------------
derive_predictors <- function(d) {
  d$latitude_band <- factor(d$latitude_band, levels = band_levels)
  d$log_precip <- log1p(pmax(d$precip_annual_mean_mm, 0))
  d$log_relief <- log1p(pmax(d$relief_5km_m, 0))
  d$asin_barren <- asin(sqrt(pmin(pmax(d$lc_barren_sparse_prop, 0), 1)))
  d$log_distance_water <- log1p(pmax(d$distance_to_modis_water_km, 0))
  d
}

primary_t <- derive_predictors(primary_raw)
ge5_t <- derive_predictors(ge5_raw)
ge10_t <- derive_predictors(ge10_raw)

predictor_manifest <- data.frame(
  variable = c(
    "wind_speed_mean_ms", "vpd_mean_kpa", "log_precip", "log_relief",
    "asin_barren", "tmean_c", "upward_sensible_heat_mean_wm2",
    "log_distance_water", "soil_water_mean_frac", "evi_mean",
    "wind_monthly_sd_ms"
  ),
  source_field = c(
    "wind_speed_mean_ms", "vpd_mean_kpa", "precip_annual_mean_mm", "relief_5km_m",
    "lc_barren_sparse_prop", "tmean_c", "upward_sensible_heat_mean_wm2",
    "distance_to_modis_water_km", "soil_water_mean_frac", "evi_mean",
    "wind_monthly_sd_ms"
  ),
  mechanism = c(
    "mean wind regime", "atmospheric dryness", "water input", "topographic heterogeneity",
    "open or sparsely vegetated habitat", "thermal regime", "surface heating and convective uplift proxy",
    "proximity to mapped surface water", "near-surface soil moisture", "vegetation productivity",
    "wind seasonality"
  ),
  transformation = c(
    "none", "none", "log1p", "log1p", "arcsine square root", "none", "none",
    "log1p", "none", "none", "none"
  ),
  priority = seq_len(11),
  stringsAsFactors = FALSE
)
write.csv(predictor_manifest,
          file.path(out_dir, "12D_predictor_manifest.csv"), row.names = FALSE)

candidate_vars <- predictor_manifest$variable

# ----------------------------
# Basic predictor audit
# ----------------------------
predictor_audit <- do.call(rbind, lapply(candidate_vars, function(v) {
  z <- primary_t[[v]]
  finite <- is.finite(z)
  data.frame(
    variable = v,
    n = length(z),
    n_finite = sum(finite),
    missing_pct = 100 * mean(!finite),
    mean = if (any(finite)) mean(z[finite]) else NA_real_,
    sd = if (sum(finite) > 1) sd(z[finite]) else NA_real_,
    min = if (any(finite)) min(z[finite]) else NA_real_,
    median = if (any(finite)) median(z[finite]) else NA_real_,
    max = if (any(finite)) max(z[finite]) else NA_real_,
    stringsAsFactors = FALSE
  )
}))
write.csv(predictor_audit,
          file.path(out_dir, "12D_predictor_distribution_audit.csv"), row.names = FALSE)

eligible_vars <- predictor_audit$variable[
  predictor_audit$missing_pct <= 5 &
    is.finite(predictor_audit$sd) & predictor_audit$sd > 0
]
if (length(eligible_vars) < 3) stop("Fewer than three environmental predictors passed the basic audit")

cor_mat <- cor(primary_t[, eligible_vars, drop = FALSE], use = "pairwise.complete.obs")
write.csv(cbind(variable = rownames(cor_mat), as.data.frame(cor_mat, check.names = FALSE)),
          file.path(out_dir, "12D_predictor_correlation_matrix.csv"), row.names = FALSE)

high_cor_pairs <- data.frame()
if (length(eligible_vars) >= 2) {
  tmp <- list()
  idx <- 1L
  for (i in seq_len(length(eligible_vars) - 1L)) {
    for (j in seq.int(i + 1L, length(eligible_vars))) {
      r <- cor_mat[i, j]
      if (is.finite(r) && abs(r) >= 0.75) {
        tmp[[idx]] <- data.frame(
          variable_1 = eligible_vars[i], variable_2 = eligible_vars[j],
          correlation = r, abs_correlation = abs(r), stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  if (length(tmp)) high_cor_pairs <- do.call(rbind, tmp)
}
if (!nrow(high_cor_pairs)) {
  high_cor_pairs <- data.frame(
    variable_1 = character(), variable_2 = character(),
    correlation = numeric(), abs_correlation = numeric(),
    stringsAsFactors = FALSE
  )
}
write.csv(high_cor_pairs,
          file.path(out_dir, "12D_high_correlation_pairs.csv"), row.names = FALSE)

# Priority-based, response-independent correlation pruning.
# Earlier variables in predictor_manifest have higher retention priority.
priority_order <- predictor_manifest$variable[predictor_manifest$variable %in% eligible_vars]
kept_corr <- character()
cor_decisions <- list()
for (v in priority_order) {
  conflict <- character()
  if (length(kept_corr)) {
    rvals <- cor_mat[v, kept_corr]
    conflict <- kept_corr[is.finite(rvals) & abs(rvals) >= 0.75]
  }
  if (!length(conflict)) {
    kept_corr <- c(kept_corr, v)
    cor_decisions[[length(cor_decisions) + 1L]] <- data.frame(
      variable = v, decision = "retain", reason = "no retained predictor at |r| >= 0.75",
      stringsAsFactors = FALSE
    )
  } else {
    cor_decisions[[length(cor_decisions) + 1L]] <- data.frame(
      variable = v, decision = "remove",
      reason = paste0("correlated with higher-priority retained predictor(s): ",
                      paste(conflict, collapse = "; ")),
      stringsAsFactors = FALSE
    )
  }
}
cor_decision_tab <- do.call(rbind, cor_decisions)
write.csv(cor_decision_tab,
          file.path(out_dir, "12D_correlation_screen_decisions.csv"), row.names = FALSE)

calc_vif <- function(d, vars) {
  if (length(vars) <= 1L) return(setNames(rep(1, length(vars)), vars))
  ans <- setNames(rep(NA_real_, length(vars)), vars)
  for (v in vars) {
    others <- setdiff(vars, v)
    fit <- lm(reformulate(others, response = v), data = d[, vars, drop = FALSE])
    r2 <- summary(fit)$r.squared
    ans[v] <- if (is.finite(r2) && r2 < 1) 1 / (1 - r2) else Inf
  }
  ans
}

prune_vif <- function(d, vars, forced = character(), threshold = 5,
                      priority_lookup = predictor_manifest$priority) {
  names(priority_lookup) <- predictor_manifest$variable
  current <- vars
  history <- list()
  repeat {
    vifs <- calc_vif(d, current)
    history[[length(history) + 1L]] <- data.frame(
      iteration = length(history) + 1L,
      variable = names(vifs), vif = as.numeric(vifs),
      stringsAsFactors = FALSE
    )
    bad <- names(vifs)[!is.na(vifs) & vifs > threshold]
    bad <- setdiff(bad, forced)
    if (!length(bad)) break
    pr <- priority_lookup[bad]
    pr[is.na(pr)] <- max(predictor_manifest$priority) + 100
    remove_v <- bad[order(pr, vifs[bad], decreasing = TRUE)][1]
    current <- setdiff(current, remove_v)
    if (length(current) <= length(forced) + 1L) break
  }
  list(vars = current, history = do.call(rbind, history), vif = calc_vif(d, current))
}

# Environment-only VIF screen.
env_vif_screen <- prune_vif(primary_t, kept_corr, forced = character(), threshold = 5)
env_vars <- env_vif_screen$vars
write.csv(env_vif_screen$history,
          file.path(out_dir, "12D_environment_vif_iterations.csv"), row.names = FALSE)

# Geography variables are constructed before combined-model VIF screening.
make_geo_raw <- function(d) {
  lat_center <- mean(d$centroid_latitude, na.rm = TRUE)
  d$lat_centered_for_vif <- d$centroid_latitude - lat_center
  d$lat_centered_sq_for_vif <- d$lat_centered_for_vif^2
  d
}
primary_vif <- make_geo_raw(primary_t)
combined_initial <- c("lat_centered_for_vif", "lat_centered_sq_for_vif", env_vars)
priority_combined <- c(
  lat_centered_for_vif = -100,
  lat_centered_sq_for_vif = -99,
  setNames(predictor_manifest$priority, predictor_manifest$variable)
)
combined_vif_screen <- prune_vif(
  primary_vif, combined_initial,
  forced = c("lat_centered_for_vif", "lat_centered_sq_for_vif"),
  threshold = 5,
  priority_lookup = priority_combined
)
combined_env_vars <- setdiff(
  combined_vif_screen$vars,
  c("lat_centered_for_vif", "lat_centered_sq_for_vif")
)
write.csv(combined_vif_screen$history,
          file.path(out_dir, "12D_combined_vif_iterations.csv"), row.names = FALSE)

selection_tab <- merge(
  predictor_manifest,
  cor_decision_tab,
  by = "variable", all.x = TRUE, sort = FALSE
)
selection_tab$retained_after_correlation <- selection_tab$variable %in% kept_corr
selection_tab$retained_environment_model <- selection_tab$variable %in% env_vars
selection_tab$retained_combined_model <- selection_tab$variable %in% combined_env_vars
write.csv(selection_tab,
          file.path(out_dir, "12D_selected_predictors.csv"), row.names = FALSE)

log_msg("Predictors eligible after basic audit: ", paste(eligible_vars, collapse = "; "))
log_msg("Predictors retained after correlation screen: ", paste(kept_corr, collapse = "; "))
log_msg("Environment-model predictors after VIF: ", paste(env_vars, collapse = "; "))
log_msg("Combined-model environmental predictors after VIF: ", paste(combined_env_vars, collapse = "; "))

if (length(env_vars) < 2 || length(combined_env_vars) < 2) {
  stop("Predictor screening retained too few variables; inspect Step 12D screening tables")
}

# ----------------------------
# Standardization helpers
# ----------------------------
fit_scaler <- function(d, vars) {
  centers <- sapply(vars, function(v) mean(d[[v]], na.rm = TRUE))
  scales <- sapply(vars, function(v) sd(d[[v]], na.rm = TRUE))
  scales[!is.finite(scales) | scales == 0] <- 1
  data.frame(variable = vars, center = as.numeric(centers), scale = as.numeric(scales),
             stringsAsFactors = FALSE)
}

apply_scaler <- function(d, scaler) {
  out <- d
  for (i in seq_len(nrow(scaler))) {
    v <- scaler$variable[i]
    out[[paste0("z_", v)]] <- (out[[v]] - scaler$center[i]) / scaler$scale[i]
  }
  out
}

prepare_scaled <- function(d, env_union, scaler = NULL) {
  d <- derive_predictors(d)
  # Latitude is scaled first; its square is then centered and scaled separately.
  if (is.null(scaler)) {
    lat_center <- mean(d$centroid_latitude, na.rm = TRUE)
    lat_scale <- sd(d$centroid_latitude, na.rm = TRUE)
    if (!is.finite(lat_scale) || lat_scale == 0) lat_scale <- 1
    d$z_lat <- (d$centroid_latitude - lat_center) / lat_scale
    lat2_center <- mean(d$z_lat^2, na.rm = TRUE)
    lat2_scale <- sd(d$z_lat^2, na.rm = TRUE)
    if (!is.finite(lat2_scale) || lat2_scale == 0) lat2_scale <- 1
    d$z_lat2 <- (d$z_lat^2 - lat2_center) / lat2_scale
    env_scaler <- fit_scaler(d, env_union)
    d <- apply_scaler(d, env_scaler)
    scaler <- list(
      lat_center = lat_center, lat_scale = lat_scale,
      lat2_center = lat2_center, lat2_scale = lat2_scale,
      env = env_scaler
    )
  } else {
    d$z_lat <- (d$centroid_latitude - scaler$lat_center) / scaler$lat_scale
    d$z_lat2 <- (d$z_lat^2 - scaler$lat2_center) / scaler$lat2_scale
    d <- apply_scaler(d, scaler$env)
  }
  list(data = d, scaler = scaler)
}

env_union <- unique(c(env_vars, combined_env_vars))
prepared_primary <- prepare_scaled(primary_raw, env_union)
model_data <- prepared_primary$data
full_scaler <- prepared_primary$scaler

scaler_tab <- rbind(
  data.frame(variable = "centroid_latitude", center = full_scaler$lat_center,
             scale = full_scaler$lat_scale, stringsAsFactors = FALSE),
  data.frame(variable = "z_lat_squared", center = full_scaler$lat2_center,
             scale = full_scaler$lat2_scale, stringsAsFactors = FALSE),
  full_scaler$env
)
write.csv(scaler_tab, file.path(out_dir, "12D_primary_scaling_parameters.csv"), row.names = FALSE)

# ----------------------------
# Model fitting helpers
# ----------------------------
response_formula <- function(success, failure, terms) {
  rhs <- if (length(terms)) paste(terms, collapse = " + ") else "1"
  as.formula(paste0("cbind(", success, ", ", failure, ") ~ ", rhs))
}

model_terms <- list(
  null = character(),
  geography = c("z_lat", "z_lat2"),
  environment = paste0("z_", env_vars),
  combined = c("z_lat", "z_lat2", paste0("z_", combined_env_vars))
)

formulas_primary <- lapply(model_terms, function(tt) {
  response_formula("ballooning_genera_primary", "non_ballooning_genera_primary", tt)
})

fit_glm_safe <- function(formula, data, family_name) {
  tryCatch({
    if (family_name == "beta_binomial") {
      if (!requireNamespace("glmmTMB", quietly = TRUE)) stop("glmmTMB not installed")
      glmmTMB::glmmTMB(
        formula = formula, data = data,
        family = glmmTMB::betabinomial(link = "logit")
      )
    } else if (family_name == "quasibinomial") {
      glm(formula = formula, data = data, family = quasibinomial(link = "logit"))
    } else {
      glm(formula = formula, data = data, family = binomial(link = "logit"))
    }
  }, error = function(e) structure(list(error = conditionMessage(e)), class = "fit_error"))
}

is_fit_error <- function(x) inherits(x, "fit_error")

pearson_dispersion <- function(fit) {
  if (is_fit_error(fit)) return(NA_real_)
  r <- tryCatch(residuals(fit, type = "pearson"), error = function(e) numeric())
  df <- tryCatch(df.residual(fit), error = function(e) NA_real_)
  if (!length(r) || !is.finite(df) || df <= 0) return(NA_real_)
  sum(r^2, na.rm = TRUE) / df
}

binomial_companion <- lapply(formulas_primary, function(fm) {
  fit_glm_safe(fm, model_data, "binomial")
})
phi_combined <- pearson_dispersion(binomial_companion$combined)

if (is.finite(phi_combined) && phi_combined > 1.5) {
  if (requireNamespace("glmmTMB", quietly = TRUE)) {
    family_selected <- "beta_binomial"
    family_reason <- paste0("Combined binomial Pearson dispersion = ", round(phi_combined, 3),
                            "; beta-binomial selected")
  } else {
    family_selected <- "quasibinomial"
    family_reason <- paste0("Combined binomial Pearson dispersion = ", round(phi_combined, 3),
                            "; glmmTMB unavailable, quasibinomial fallback selected")
  }
} else {
  family_selected <- "binomial"
  family_reason <- paste0("Combined binomial Pearson dispersion = ", round(phi_combined, 3),
                          "; standard binomial retained")
}
log_msg("Selected model family: ", family_selected)
log_msg("Family decision: ", family_reason)

primary_models <- lapply(formulas_primary, function(fm) {
  fit_glm_safe(fm, model_data, family_selected)
})

# If beta-binomial failed globally, fall back transparently to quasibinomial.
if (family_selected == "beta_binomial" && any(vapply(primary_models, is_fit_error, logical(1)))) {
  log_msg("WARNING: At least one beta-binomial model failed. Falling back to quasibinomial for all primary models.")
  family_selected <- "quasibinomial"
  family_reason <- paste0(family_reason, "; beta-binomial fitting failure triggered quasibinomial fallback")
  primary_models <- lapply(formulas_primary, function(fm) fit_glm_safe(fm, model_data, family_selected))
}

model_nobs <- function(fit) {
  if (is_fit_error(fit)) return(NA_integer_)
  tryCatch(nobs(fit), error = function(e) NA_integer_)
}

model_loglik <- function(fit) {
  if (is_fit_error(fit)) return(NA_real_)
  tryCatch(as.numeric(logLik(fit)), error = function(e) NA_real_)
}

model_k <- function(fit) {
  if (is_fit_error(fit)) return(NA_real_)
  ll <- tryCatch(logLik(fit), error = function(e) NULL)
  if (is.null(ll)) return(NA_real_)
  as.numeric(attr(ll, "df"))
}

model_aic <- function(fit) {
  if (is_fit_error(fit)) return(NA_real_)
  tryCatch(AIC(fit), error = function(e) NA_real_)
}

model_aicc <- function(fit) {
  a <- model_aic(fit)
  k <- model_k(fit)
  n <- model_nobs(fit)
  if (!is.finite(a) || !is.finite(k) || !is.finite(n) || n <= k + 1) return(NA_real_)
  a + (2 * k * (k + 1)) / (n - k - 1)
}

model_converged <- function(fit) {
  if (is_fit_error(fit)) return(FALSE)
  if (inherits(fit, "glmmTMB")) {
    ok <- tryCatch(isTRUE(fit$sdr$pdHess), error = function(e) FALSE)
    return(ok)
  }
  if (inherits(fit, "glm")) return(isTRUE(fit$converged))
  TRUE
}

model_comparison <- do.call(rbind, lapply(names(primary_models), function(nm) {
  fit <- primary_models[[nm]]
  data.frame(
    model = nm,
    family = family_selected,
    n_cells = model_nobs(fit),
    n_terms = length(model_terms[[nm]]),
    converged = model_converged(fit),
    logLik = model_loglik(fit),
    AIC = model_aic(fit),
    AICc = model_aicc(fit),
    pearson_dispersion = pearson_dispersion(fit),
    error = if (is_fit_error(fit)) fit$error else "",
    stringsAsFactors = FALSE
  )
}))
if (any(is.finite(model_comparison$AICc))) {
  model_comparison$delta_AICc <- model_comparison$AICc - min(model_comparison$AICc, na.rm = TRUE)
} else {
  model_comparison$delta_AICc <- NA_real_
}
write.csv(model_comparison,
          file.path(out_dir, "12D_primary_model_comparison.csv"), row.names = FALSE)

extract_coef <- function(fit, model_name, dataset_name = "primary", response_name = "primary") {
  if (is_fit_error(fit)) return(data.frame())
  if (inherits(fit, "glmmTMB")) {
    sm <- summary(fit)$coefficients$cond
  } else {
    sm <- summary(fit)$coefficients
  }
  out <- data.frame(
    dataset = dataset_name,
    response = response_name,
    model = model_name,
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

coef_primary <- do.call(rbind, lapply(names(primary_models), function(nm) {
  extract_coef(primary_models[[nm]], nm)
}))
write.csv(coef_primary,
          file.path(out_dir, "12D_primary_model_coefficients.csv"), row.names = FALSE)

# ----------------------------
# Leave-one-latitude-band-out cross-validation
# ----------------------------
clip_prob <- function(p, eps = 1e-8) pmin(pmax(p, eps), 1 - eps)

cv_one_model <- function(raw_data, model_name, raw_env_vars, family_name) {
  rows <- list()
  for (b in band_levels) {
    train_raw <- raw_data[as.character(raw_data$latitude_band) != b, , drop = FALSE]
    test_raw <- raw_data[as.character(raw_data$latitude_band) == b, , drop = FALSE]
    if (!nrow(train_raw) || !nrow(test_raw)) next

    prep_train <- prepare_scaled(train_raw, env_union)
    train <- prep_train$data
    test <- prepare_scaled(test_raw, env_union, scaler = prep_train$scaler)$data

    terms <- switch(
      model_name,
      null = character(),
      geography = c("z_lat", "z_lat2"),
      environment = paste0("z_", env_vars),
      combined = c("z_lat", "z_lat2", paste0("z_", combined_env_vars))
    )
    fm <- response_formula("ballooning_genera_primary", "non_ballooning_genera_primary", terms)
    fit <- fit_glm_safe(fm, train, family_name)
    if (is_fit_error(fit)) {
      rows[[length(rows) + 1L]] <- data.frame(
        model = model_name, held_out_band = b, n_cells = nrow(test),
        total_classified_genera = sum(test$classified_genera_primary),
        weighted_brier = NA_real_, weighted_mae = NA_real_,
        binomial_log_score_per_genus = NA_real_, error = fit$error,
        stringsAsFactors = FALSE
      )
      next
    }
    p <- tryCatch(predict(fit, newdata = test, type = "response"), error = function(e) rep(NA_real_, nrow(test)))
    p <- clip_prob(as.numeric(p))
    y <- test$ballooning_genera_primary
    n <- test$classified_genera_primary
    obs <- y / n
    valid <- is.finite(p) & is.finite(obs) & is.finite(n) & n > 0
    if (!any(valid)) next
    w <- n[valid]
    rows[[length(rows) + 1L]] <- data.frame(
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

    pred_rows <- data.frame(
      grid_cell_id = test$grid_cell_id,
      latitude_band = as.character(test$latitude_band),
      model = model_name,
      observed_ballooning = y,
      denominator = n,
      observed_proportion = obs,
      predicted_proportion = p,
      stringsAsFactors = FALSE
    )
    attr(rows[[length(rows)]], "predictions") <- pred_rows
  }
  rows
}

cv_nested <- lapply(names(model_terms), function(nm) {
  cv_one_model(primary_raw, nm, env_union, family_selected)
})
names(cv_nested) <- names(model_terms)

cv_folds <- do.call(rbind, unlist(lapply(cv_nested, function(x) {
  lapply(x, function(z) {
    attr(z, "predictions") <- NULL
    z
  })
}), recursive = FALSE))
write.csv(cv_folds, file.path(out_dir, "12D_leave_one_band_out_cv_by_fold.csv"), row.names = FALSE)

cv_predictions <- do.call(rbind, unlist(lapply(cv_nested, function(x) {
  lapply(x, function(z) attr(z, "predictions"))
}), recursive = FALSE))
if (is.null(cv_predictions)) cv_predictions <- data.frame()
write.csv(cv_predictions, file.path(out_dir, "12D_leave_one_band_out_cv_predictions.csv"), row.names = FALSE)

cv_summary <- do.call(rbind, lapply(split(cv_folds, cv_folds$model), function(d) {
  valid <- is.finite(d$weighted_brier) & is.finite(d$total_classified_genera)
  if (!any(valid)) {
    return(data.frame(model = unique(d$model), folds_completed = 0,
                      total_classified_genera = NA_real_, weighted_brier = NA_real_,
                      weighted_mae = NA_real_, binomial_log_score_per_genus = NA_real_,
                      stringsAsFactors = FALSE))
  }
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
cv_summary <- cv_summary[match(names(model_terms), cv_summary$model), , drop = FALSE]
write.csv(cv_summary, file.path(out_dir, "12D_leave_one_band_out_cv_summary.csv"), row.names = FALSE)

# ----------------------------
# Residual spatial autocorrelation and influence
# ----------------------------
moran_knn <- function(resid, x, y, k = 4L, permutations = 999L, seed = 1204L) {
  ok <- is.finite(resid) & is.finite(x) & is.finite(y)
  resid <- resid[ok]
  x <- x[ok]
  y <- y[ok]
  n <- length(resid)
  if (n <= k + 2L) return(data.frame(n = n, k = k, moran_I = NA_real_, p_permutation = NA_real_))
  coords <- cbind(x, y)
  dm <- as.matrix(dist(coords))
  diag(dm) <- Inf
  W <- matrix(0, n, n)
  for (i in seq_len(n)) {
    nn <- order(dm[i, ])[seq_len(k)]
    W[i, nn] <- 1 / k
  }
  z <- resid - mean(resid)
  S0 <- sum(W)
  calc_I <- function(v) (n / S0) * (sum(W * outer(v, v)) / sum(v^2))
  obs <- calc_I(z)
  set.seed(seed)
  perm <- replicate(permutations, calc_I(sample(z, replace = FALSE)))
  p <- (1 + sum(abs(perm) >= abs(obs))) / (permutations + 1)
  data.frame(n = n, k = k, moran_I = obs, p_permutation = p,
             permutations = permutations, stringsAsFactors = FALSE)
}

combined_fit <- primary_models$combined
if (!is_fit_error(combined_fit)) {
  combined_pred <- as.numeric(predict(combined_fit, newdata = model_data, type = "response"))
  combined_resid <- as.numeric(residuals(combined_fit, type = "pearson"))
  moran_tab <- moran_knn(
    combined_resid, model_data$centroid_x_m, model_data$centroid_y_m,
    k = 4L, permutations = 999L
  )
  write.csv(moran_tab, file.path(out_dir, "12D_combined_model_residual_moran.csv"), row.names = FALSE)

  fitted_tab <- data.frame(
    grid_cell_id = model_data$grid_cell_id,
    latitude_band = as.character(model_data$latitude_band),
    centroid_latitude = model_data$centroid_latitude,
    centroid_longitude = model_data$centroid_longitude,
    ballooning_genera = model_data$ballooning_genera_primary,
    non_ballooning_genera = model_data$non_ballooning_genera_primary,
    classified_genera = model_data$classified_genera_primary,
    observed_proportion = model_data$ballooning_genera_primary / model_data$classified_genera_primary,
    fitted_proportion = combined_pred,
    pearson_residual = combined_resid,
    stringsAsFactors = FALSE
  )
  write.csv(fitted_tab, file.path(out_dir, "12D_combined_model_fitted_cells.csv"), row.names = FALSE)
} else {
  combined_pred <- rep(NA_real_, nrow(model_data))
  combined_resid <- rep(NA_real_, nrow(model_data))
  moran_tab <- data.frame(n = nrow(model_data), k = 4, moran_I = NA_real_, p_permutation = NA_real_, permutations = 999)
  write.csv(moran_tab, file.path(out_dir, "12D_combined_model_residual_moran.csv"), row.names = FALSE)
}

# Influence is evaluated from the same-formula binomial companion because standard
# Cook's distance is not generally defined for a beta-binomial glmmTMB fit.
companion_combined <- binomial_companion$combined
if (!is_fit_error(companion_combined)) {
  influence_tab <- data.frame(
    grid_cell_id = model_data$grid_cell_id,
    latitude_band = as.character(model_data$latitude_band),
    classified_genera = model_data$classified_genera_primary,
    cooks_distance_binomial_companion = cooks.distance(companion_combined),
    leverage_binomial_companion = hatvalues(companion_combined),
    pearson_residual_binomial_companion = residuals(companion_combined, type = "pearson"),
    stringsAsFactors = FALSE
  )
  influence_tab$cook_threshold_4_over_n <- 4 / nrow(influence_tab)
  influence_tab$flag_cook <- influence_tab$cooks_distance_binomial_companion > influence_tab$cook_threshold_4_over_n
  influence_tab <- influence_tab[order(influence_tab$cooks_distance_binomial_companion, decreasing = TRUE), ]
  write.csv(influence_tab, file.path(out_dir, "12D_combined_model_influence_audit.csv"), row.names = FALSE)
}

# ----------------------------
# Sensitivity fits using the primary-selected combined formula
# ----------------------------
response_specs <- data.frame(
  response = c("primary", "taxonomy_strict", "low_conf_exclusion"),
  success = c("ballooning_genera_primary", "ballooning_genera_taxonomy_strict", "ballooning_genera_low_conf_exclusion"),
  failure = c("non_ballooning_genera_primary", "non_ballooning_genera_taxonomy_strict", "non_ballooning_genera_low_conf_exclusion"),
  denominator = c("classified_genera_primary", "classified_genera_taxonomy_strict", "classified_genera_low_conf_exclusion"),
  stringsAsFactors = FALSE
)

dataset_specs <- list(primary = primary_raw, denominator_ge5 = ge5_raw, denominator_ge10 = ge10_raw)

sensitivity_models <- list()
sensitivity_summary_rows <- list()
sensitivity_coef_rows <- list()
for (ds_name in names(dataset_specs)) {
  raw_ds <- dataset_specs[[ds_name]]
  for (i in seq_len(nrow(response_specs))) {
    rs <- response_specs[i, ]
    keep <- is.finite(raw_ds[[rs$denominator]]) & raw_ds[[rs$denominator]] > 0
    d0 <- raw_ds[keep, , drop = FALSE]
    if (nrow(d0) < length(model_terms$combined) + 8L) {
      sensitivity_summary_rows[[length(sensitivity_summary_rows) + 1L]] <- data.frame(
        dataset = ds_name, response = rs$response, n_cells = nrow(d0),
        family = family_selected, converged = FALSE, AICc = NA_real_,
        pearson_dispersion = NA_real_, error = "insufficient cells for combined model",
        stringsAsFactors = FALSE
      )
      next
    }
    prep <- prepare_scaled(d0, env_union)
    d <- prep$data
    fm <- response_formula(rs$success, rs$failure, model_terms$combined)
    fit <- fit_glm_safe(fm, d, family_selected)
    key <- paste(ds_name, rs$response, sep = "__")
    sensitivity_models[[key]] <- fit
    sensitivity_summary_rows[[length(sensitivity_summary_rows) + 1L]] <- data.frame(
      dataset = ds_name, response = rs$response, n_cells = nrow(d),
      family = family_selected, converged = model_converged(fit),
      AICc = model_aicc(fit), pearson_dispersion = pearson_dispersion(fit),
      error = if (is_fit_error(fit)) fit$error else "",
      stringsAsFactors = FALSE
    )
    cc <- extract_coef(fit, "combined", dataset_name = ds_name, response_name = rs$response)
    if (nrow(cc)) sensitivity_coef_rows[[length(sensitivity_coef_rows) + 1L]] <- cc
  }
}

sensitivity_summary <- do.call(rbind, sensitivity_summary_rows)
sensitivity_coef <- if (length(sensitivity_coef_rows)) do.call(rbind, sensitivity_coef_rows) else data.frame()
write.csv(sensitivity_summary, file.path(out_dir, "12D_sensitivity_model_summary.csv"), row.names = FALSE)
write.csv(sensitivity_coef, file.path(out_dir, "12D_sensitivity_model_coefficients.csv"), row.names = FALSE)

# ----------------------------
# Figures
# ----------------------------
fig_dir <- file.path(out_dir, "figures")

save_both <- function(stem, width = 8, height = 6, plot_fun) {
  png(file.path(fig_dir, paste0(stem, ".png")), width = width, height = height,
      units = "in", res = 400)
  plot_fun()
  dev.off()
  pdf(file.path(fig_dir, paste0(stem, ".pdf")), width = width, height = height,
      useDingbats = FALSE)
  plot_fun()
  dev.off()
}

save_both("12D_predictor_correlation_heatmap", 10, 9, function() {
  par(mar = c(10, 10, 3, 2))
  cm <- cor_mat[rev(rownames(cor_mat)), colnames(cor_mat), drop = FALSE]
  image(seq_len(ncol(cm)), seq_len(nrow(cm)), t(cm), zlim = c(-1, 1),
        axes = FALSE, xlab = "", ylab = "", main = "Environmental predictor correlations")
  axis(1, at = seq_len(ncol(cm)), labels = colnames(cm), las = 2, cex.axis = 0.65)
  axis(2, at = seq_len(nrow(cm)), labels = rownames(cm), las = 2, cex.axis = 0.65)
  box()
})

if (nrow(cv_summary) && any(is.finite(cv_summary$binomial_log_score_per_genus))) {
  save_both("12D_spatial_cv_model_comparison", 8, 6, function() {
    d <- cv_summary[is.finite(cv_summary$binomial_log_score_per_genus), , drop = FALSE]
    ord <- order(d$binomial_log_score_per_genus)
    d <- d[ord, ]
    par(mar = c(5, 9, 3, 1))
    barplot(d$binomial_log_score_per_genus, names.arg = d$model, horiz = TRUE, las = 1,
            xlab = "Leave-one-band-out binomial log score per classified genus\n(lower is better)",
            main = "Spatial cross-validation model comparison")
  })
}

if (!is_fit_error(combined_fit)) {
  observed_prop <- model_data$ballooning_genera_primary / model_data$classified_genera_primary
  save_both("12D_combined_observed_vs_fitted", 7.5, 6.5, function() {
    cex <- 0.7 + 1.8 * sqrt(model_data$classified_genera_primary / max(model_data$classified_genera_primary))
    plot(observed_prop, combined_pred, pch = 21, cex = cex,
         xlim = c(0, 1), ylim = c(0, 1),
         xlab = "Observed ballooning-genus proportion",
         ylab = "Fitted ballooning-genus proportion",
         main = "Combined environmental model")
    abline(0, 1, lty = 2)
  })

  save_both("12D_combined_residual_spatial_pattern", 7.5, 8, function() {
    pal <- hcl.colors(100, "Blue-Red 3")
    cuts <- cut(combined_resid, breaks = 100, include.lowest = TRUE, labels = FALSE)
    plot(model_data$centroid_longitude, model_data$centroid_latitude,
         pch = 21, bg = pal[cuts], cex = 1.1,
         xlab = "Longitude", ylab = "Latitude",
         main = "Pearson residuals from combined model")
  })

  cc <- coef_primary[coef_primary$model == "combined" & coef_primary$term != "(Intercept)", , drop = FALSE]
  if (nrow(cc)) {
    cc <- cc[order(cc$odds_ratio), ]
    save_both("12D_combined_model_odds_ratios", 8, 6.5, function() {
      par(mar = c(5, 11, 3, 2))
      yy <- seq_len(nrow(cc))
      xlim <- range(c(cc$ci_low, cc$ci_high), finite = TRUE)
      plot(cc$odds_ratio, yy, log = "x", xlim = xlim, yaxt = "n", pch = 19,
           xlab = "Odds ratio per 1 SD increase (log scale)", ylab = "",
           main = "Combined model standardized effects")
      segments(cc$ci_low, yy, cc$ci_high, yy)
      axis(2, at = yy, labels = cc$term, las = 2, cex.axis = 0.75)
      abline(v = 1, lty = 2)
    })
  }
}

# ----------------------------
# Save models, formulas, and interpretation aids
# ----------------------------
saveRDS(
  list(
    version = version,
    family_selected = family_selected,
    family_reason = family_reason,
    selected_environment_predictors = env_vars,
    selected_combined_predictors = combined_env_vars,
    model_terms = model_terms,
    scaler = full_scaler,
    models = primary_models,
    binomial_companion_models = binomial_companion,
    sensitivity_models = sensitivity_models
  ),
  file.path(out_dir, "models", "12D_model_objects.rds")
)

formula_lines <- c(
  paste0("Selected family: ", family_selected),
  paste0("Family reason: ", family_reason),
  "",
  paste0("Null: ", deparse(formulas_primary$null)),
  paste0("Geography: ", deparse(formulas_primary$geography)),
  paste0("Environment: ", deparse(formulas_primary$environment)),
  paste0("Combined: ", deparse(formulas_primary$combined)),
  "",
  paste0("Environment variables: ", paste(env_vars, collapse = "; ")),
  paste0("Combined environmental variables: ", paste(combined_env_vars, collapse = "; "))
)
writeLines(formula_lines, file.path(out_dir, "12D_model_formulas.txt"))

# Rank by spatial CV log score; do not treat a small numerical difference as biological certainty.
cv_valid <- cv_summary[is.finite(cv_summary$binomial_log_score_per_genus), , drop = FALSE]
if (nrow(cv_valid)) {
  cv_valid <- cv_valid[order(cv_valid$binomial_log_score_per_genus), ]
  best_cv <- cv_valid$model[1]
  best_score <- cv_valid$binomial_log_score_per_genus[1]
} else {
  best_cv <- NA_character_
  best_score <- NA_real_
}

recommendation <- c(
  "STEP 12D MODEL INTERPRETATION NOTE",
  "==================================",
  "",
  paste0("Primary model family: ", family_selected),
  paste0("Reason: ", family_reason),
  paste0("Lowest leave-one-latitude-band-out log score: ", best_cv,
         if (is.finite(best_score)) paste0(" (", signif(best_score, 5), ")") else ""),
  "",
  "The spatial cross-validation ranking should be considered together with:",
  "  * coefficient uncertainty and effect direction;",
  "  * residual Moran's I;",
  "  * stability across denominator, taxonomy-strict, and LOW-confidence sensitivities;",
  "  * biological interpretability; and",
  "  * the fact that cells share genera, so this is a community-composition association model rather than independent genus trials.",
  "",
  "Do not produce an extrapolative manuscript prediction map. Step 12E and Spatial+ evaluate whether any environmental association remains after model reduction and spatial deconfounding."
)
writeLines(recommendation, file.path(out_dir, "12D_model_interpretation_note.txt"))

validation <- data.frame(
  check = c(
    "candidate_sets_consistent_with_step12C",
    "primary_unique_cell_ids",
    "primary_response_counts_valid",
    "environment_predictors_retained",
    "combined_predictors_retained",
    "all_primary_models_fit",
    "all_primary_models_converged",
    "five_spatial_cv_folds_per_model",
    "combined_residual_moran_written",
    "sensitivity_models_attempted"
  ),
  passed = c(
    candidate_sets_consistent,
    !anyDuplicated(primary_raw$grid_cell_id),
    all(primary_raw$ballooning_genera_primary + primary_raw$non_ballooning_genera_primary == primary_raw$classified_genera_primary),
    length(env_vars) >= 2,
    length(combined_env_vars) >= 2,
    all(!vapply(primary_models, is_fit_error, logical(1))),
    all(vapply(primary_models, model_converged, logical(1))),
    all(table(cv_folds$model) == 5),
    file.exists(file.path(out_dir, "12D_combined_model_residual_moran.csv")),
    nrow(sensitivity_summary) >= 6
  ),
  detail = c(
    paste(names(actual_counts), actual_counts, collapse = "; "),
    paste0("rows=", nrow(primary_raw), "; unique_ids=", length(unique(primary_raw$grid_cell_id))),
    paste0("denominator range=", min(primary_raw$classified_genera_primary), " to ", max(primary_raw$classified_genera_primary)),
    paste(env_vars, collapse = "; "),
    paste(combined_env_vars, collapse = "; "),
    paste(names(primary_models), vapply(primary_models, function(x) if (is_fit_error(x)) x$error else "OK", character(1)), collapse = "; "),
    paste(names(primary_models), vapply(primary_models, model_converged, logical(1)), collapse = "; "),
    paste(names(table(cv_folds$model)), as.integer(table(cv_folds$model)), collapse = "; "),
    paste0("Moran file exists=", file.exists(file.path(out_dir, "12D_combined_model_residual_moran.csv"))),
    paste0("sensitivity rows=", nrow(sensitivity_summary))
  ),
  stringsAsFactors = FALSE
)
write.csv(validation, file.path(out_dir, "12D_validation.csv"), row.names = FALSE)

sink(file.path(out_dir, "12D_session_info.txt"))
print(sessionInfo())
sink()

readme <- c(
  "STEP 12D OUTPUTS",
  "================",
  "",
  "This step screens environmental predictors without using the response, compares null, geography-only,",
  "environment-only, and combined grouped-binomial models, audits overdispersion, and evaluates spatial",
  "transferability by leaving out one manuscript latitude band at a time.",
  "",
  "Key files:",
  "  12D_selected_predictors.csv",
  "  12D_primary_model_comparison.csv",
  "  12D_primary_model_coefficients.csv",
  "  12D_leave_one_band_out_cv_summary.csv",
  "  12D_combined_model_residual_moran.csv",
  "  12D_combined_model_influence_audit.csv",
  "  12D_sensitivity_model_summary.csv",
  "  12D_sensitivity_model_coefficients.csv",
  "  12D_model_interpretation_note.txt",
  "  12D_validation.csv",
  "",
  "Model-family behavior:",
  "  A beta-binomial model is used when the companion binomial is overdispersed and glmmTMB is installed.",
  "  If glmmTMB is unavailable, the script completes with quasibinomial inference and relies on spatial CV",
  "  rather than AIC for model comparison.",
  "",
  "No prediction surface is created in Step 12D. Mapping follows only after diagnostic review."
)
writeLines(readme, file.path(out_dir, "README_12D_OUTPUTS.txt"))

log_msg("Primary binomial dispersion before family selection: ", signif(phi_combined, 5))
log_msg("Selected family: ", family_selected)
log_msg("Primary models fit: ", paste(names(primary_models), collapse = "; "))
log_msg("Leave-one-band-out best numerical score: ", best_cv)
log_msg("Validation checks passed: ", sum(validation$passed), "/", nrow(validation))
if (!all(validation$passed)) {
  log_msg("STEP 12D COMPLETED WITH FLAGS — inspect 12D_validation.csv")
} else {
  log_msg("STEP 12D COMPLETED SUCCESSFULLY")
}
