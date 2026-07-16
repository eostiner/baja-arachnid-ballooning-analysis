#!/usr/bin/env Rscript

# ==============================================================================
# STEP 12K1 — SINGLE-PANEL OBSERVED BALLOONING + SPATIAL+ CONTRIBUTION MAP
#
# Purpose:
#   Create one publication map in which:
#     1) 25-km cell fill = observed proportion of classified genera capable
#        of ballooning.
#     2) centroid ring colour = signed combined Spatial+ contribution from
#        local VPD and wind-seasonality departures.
#     3) centroid ring size = absolute magnitude of that contribution.
#
# The environmental contribution is calculated on the probability scale as:
#
#   delta_p = P(full Spatial+ model) -
#             P(same model with VPD and wind Spatial+ terms set to zero)
#
# This is a descriptive model contribution within sampled cells. It is not a
# habitat-suitability surface, a refugia map, or a causal estimate.
#
# Usage:
#   Rscript 12K1_make_combined_spatialplus_map.R \
#     "$HOME/Desktop/Baja_Ballooning_Pipeline"
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) {
  normalizePath(path.expand(args[[1]]), mustWork = TRUE)
} else {
  normalizePath(path.expand("~/Desktop/Baja_Ballooning_Pipeline"),
                mustWork = TRUE)
}

required_packages <- c(
  "sf", "dplyr", "readr", "ggplot2", "mgcv", "maps", "scales", "viridis"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing R packages: ", paste(missing_packages, collapse = ", "), "\n",
      "Install them with:\n",
      "install.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(mgcv)
  library(maps)
  library(scales)
  library(viridis)
})

options(stringsAsFactors = FALSE)
set.seed(20260714)

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

step12c_dir <- file.path(
  project_root, "04_analysis", "12C_cell_environment_model_table"
)
step12f_dir <- file.path(
  project_root, "04_analysis", "12F_environment_predictor_scores"
)
step12k_dir <- file.path(
  project_root, "04_analysis", "12K_spatial_plus_trait_composition"
)

response_path <- file.path(
  step12c_dir, "12C_primary_glm_candidate_table.csv"
)
environment_path <- file.path(
  step12f_dir, "12F_environment_predictor_scores_by_cell.csv"
)
cells_path <- file.path(
  step12c_dir, "12C_occupied_cells_response.geojson"
)
archived_coefficients_path <- file.path(
  step12k_dir, "12K_primary_model_coefficients.csv"
)
archived_decomposition_path <- file.path(
  step12k_dir, "12K_covariate_spatial_decomposition.csv"
)

required_inputs <- c(
  response_path,
  environment_path,
  cells_path,
  archived_coefficients_path,
  archived_decomposition_path
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop(
    paste(
      "Required input files were not found:",
      paste(missing_inputs, collapse = "\n"),
      sep = "\n"
    ),
    call. = FALSE
  )
}

output_dir <- file.path(
  project_root, "04_analysis", "12K1_combined_spatialplus_map"
)
figure_dir <- file.path(output_dir, "figures")
archive_root <- file.path(project_root, "08_archive")

# Preserve any previous 12K1 output.
if (dir.exists(output_dir)) {
  dir.create(archive_root, recursive = TRUE, showWarnings = FALSE)
  timestamp <- format(Sys.time(), "%Y%m%dT%H%M%S")
  archive_path <- file.path(
    archive_root, paste0("12K1_combined_spatialplus_map_", timestamp)
  )
  if (!file.rename(output_dir, archive_path)) {
    stop("Could not archive the previous 12K1 output.", call. = FALSE)
  }
  message("Archived prior 12K1 output: ", archive_path)
}

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(output_dir, "12K1_analysis_log.txt")
log_connection <- file(log_path, open = "wt")
sink(log_connection, type = "output", split = TRUE)
sink(log_connection, type = "message", append = TRUE)

on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(type = "output"), silent = TRUE)
  try(close(log_connection), silent = TRUE)
}, add = TRUE)

cat("STEP 12K1 COMBINED SPATIAL+ MAP STARTED\n")
cat("Version: 12K1_v2_2026-07-16\n")
cat("Project root:", project_root, "\n")
cat("Output directory:", output_dir, "\n\n")

# ------------------------------------------------------------------------------
# Read and align model data
# ------------------------------------------------------------------------------

response <- read_csv(response_path, show_col_types = FALSE)
environment <- read_csv(environment_path, show_col_types = FALSE)

needed_response <- c(
  "grid_cell_id",
  "latitude_band",
  "easting_km",
  "northing_km",
  "ballooning_genera_primary",
  "non_ballooning_genera_primary",
  "classified_genera_primary",
  "ballooning_proportion_primary"
)

needed_environment <- c(
  "grid_cell_id",
  "vpd_z",
  "wind_seasonality_z"
)

missing_response_fields <- setdiff(needed_response, names(response))
missing_environment_fields <- setdiff(needed_environment, names(environment))

if (length(missing_response_fields) > 0) {
  stop(
    "Missing response fields: ",
    paste(missing_response_fields, collapse = ", "),
    call. = FALSE
  )
}
if (length(missing_environment_fields) > 0) {
  stop(
    "Missing environmental fields: ",
    paste(missing_environment_fields, collapse = ", "),
    call. = FALSE
  )
}

dat <- response |>
  select(all_of(needed_response)) |>
  inner_join(
    environment |> select(all_of(needed_environment)),
    by = "grid_cell_id"
  ) |>
  filter(
    is.finite(easting_km),
    is.finite(northing_km),
    is.finite(vpd_z),
    is.finite(wind_seasonality_z),
    classified_genera_primary > 0
  ) |>
  arrange(grid_cell_id)

if (anyDuplicated(dat$grid_cell_id)) {
  stop("Duplicate grid_cell_id values remain after the join.", call. = FALSE)
}

cat("Complete modeled cells:", nrow(dat), "\n")
cat(
  "Classified genus-cell denominator:",
  sum(dat$classified_genera_primary),
  "\n\n"
)

if (nrow(dat) == 0L) {
  stop("No complete modeled cells remain after the Step 12K1 join.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# Reconstruct the Spatial+ covariates
# ------------------------------------------------------------------------------

basis_k <- 30

vpd_spatial_model <- gam(
  vpd_z ~ s(easting_km, northing_km, bs = "tp", k = basis_k),
  data = dat,
  method = "REML"
)

wind_spatial_model <- gam(
  wind_seasonality_z ~
    s(easting_km, northing_km, bs = "tp", k = basis_k),
  data = dat,
  method = "REML"
)

standardize_residual <- function(model) {
  x <- residuals(model, type = "response")
  as.numeric((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

dat$vpd_spplus <- standardize_residual(vpd_spatial_model)
dat$wind_seasonality_spplus <- standardize_residual(wind_spatial_model)

# ------------------------------------------------------------------------------
# Refit the primary Spatial+ outcome model
# ------------------------------------------------------------------------------

spatial_plus_model <- gam(
  cbind(
    ballooning_genera_primary,
    non_ballooning_genera_primary
  ) ~
    s(easting_km, northing_km, bs = "tp", k = basis_k) +
    vpd_spplus +
    wind_seasonality_spplus,
  family = binomial(link = "logit"),
  data = dat,
  method = "REML"
)

model_coefficients <- coef(spatial_plus_model)

required_terms <- c("vpd_spplus", "wind_seasonality_spplus")
if (!all(required_terms %in% names(model_coefficients))) {
  stop(
    "The required Spatial+ coefficients were not recovered.",
    call. = FALSE
  )
}

beta_vpd <- unname(model_coefficients[["vpd_spplus"]])
beta_wind <- unname(model_coefficients[["wind_seasonality_spplus"]])

# Contribution on the logit scale.
dat$spatialplus_contribution_logit <-
  beta_vpd * dat$vpd_spplus +
  beta_wind * dat$wind_seasonality_spplus

# Translate the contribution into the change in modeled probability while
# retaining each cell's fitted spatial smooth and intercept.
full_link <- as.numeric(
  predict(spatial_plus_model, newdata = dat, type = "link")
)
base_link <- full_link - dat$spatialplus_contribution_logit

dat$predicted_probability_full <- plogis(full_link)
dat$predicted_probability_without_local_environment <- plogis(base_link)

dat$spatialplus_contribution_probability <-
  dat$predicted_probability_full -
  dat$predicted_probability_without_local_environment

dat$spatialplus_contribution_percentage_points <-
  100 * dat$spatialplus_contribution_probability

dat$absolute_spatialplus_contribution_percentage_points <-
  abs(dat$spatialplus_contribution_percentage_points)

# ------------------------------------------------------------------------------
# Validate against archived Step 12K results
# ------------------------------------------------------------------------------

archived_coefficients <- read_csv(
  archived_coefficients_path,
  show_col_types = FALSE
) |>
  filter(
    model == "spatial_plus_core",
    term %in% required_terms
  ) |>
  select(term, archived_estimate = estimate)

current_coefficients <- tibble(
  term = required_terms,
  current_estimate = c(beta_vpd, beta_wind)
)

coefficient_validation <- current_coefficients |>
  left_join(archived_coefficients, by = "term") |>
  mutate(
    absolute_difference = abs(current_estimate - archived_estimate),
    within_tolerance_0_01 = absolute_difference <= 0.01
  )

archived_decomposition <- read_csv(
  archived_decomposition_path,
  show_col_types = FALSE
) |>
  filter(predictor %in% c("vpd_z", "wind_seasonality_z")) |>
  select(predictor, archived_spatial_r_squared = spatial_r_squared)

current_decomposition <- tibble(
  predictor = c("vpd_z", "wind_seasonality_z"),
  current_spatial_r_squared = c(
    summary(vpd_spatial_model)$r.sq,
    summary(wind_spatial_model)$r.sq
  )
)

decomposition_validation <- current_decomposition |>
  left_join(archived_decomposition, by = "predictor") |>
  mutate(
    absolute_difference =
      abs(current_spatial_r_squared - archived_spatial_r_squared),
    within_tolerance_0_01 = absolute_difference <= 0.01
  )

write_csv(
  coefficient_validation,
  file.path(output_dir, "12K1_coefficient_validation.csv")
)
write_csv(
  decomposition_validation,
  file.path(output_dir, "12K1_spatial_decomposition_validation.csv")
)

cat("Refitted Spatial+ coefficients:\n")
print(coefficient_validation)
cat("\nSpatial decomposition validation:\n")
print(decomposition_validation)
cat("\n")

if (!all(coefficient_validation$within_tolerance_0_01, na.rm = TRUE)) {
  warning(
    "At least one refitted coefficient differs from the archived Step 12K ",
    "estimate by more than 0.01. Check the mgcv version and input files."
  )
}
if (!all(decomposition_validation$within_tolerance_0_01, na.rm = TRUE)) {
  warning(
    "At least one spatial decomposition R-squared differs from the archived ",
    "Step 12K value by more than 0.01."
  )
}

# ------------------------------------------------------------------------------
# Join to 25-km polygons
# ------------------------------------------------------------------------------

cells_all <- st_read(cells_path, quiet = TRUE)

if (!"grid_cell_id" %in% names(cells_all)) {
  stop("The cell geometry file lacks grid_cell_id.", call. = FALSE)
}

if (is.na(st_crs(cells_all))) {
  st_crs(cells_all) <- 4326
}
cells_all <- st_transform(cells_all, 4326)

map_cells <- cells_all |>
  select(grid_cell_id, geometry) |>
  inner_join(dat, by = "grid_cell_id")

if (nrow(map_cells) != nrow(dat)) {
  stop(
    "Not every modeled cell matched a polygon geometry.",
    call. = FALSE
  )
}

map_centroids <- st_point_on_surface(map_cells)

# Background country outline, provided locally by the maps package.
world_map <- maps::map("world", plot = FALSE, fill = TRUE)
world_sf <- st_as_sf(world_map)
st_crs(world_sf) <- 4326

mexico <- world_sf |>
  filter(grepl("^Mexico", ID)) |>
  st_make_valid()

# ------------------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------------------

max_abs_pp <- max(
  map_centroids$absolute_spatialplus_contribution_percentage_points,
  na.rm = TRUE
)

# Stable size breaks even when effects are extremely small.
if (!is.finite(max_abs_pp) || max_abs_pp <= 0) {
  size_breaks <- c(0)
} else {
  size_breaks <- unique(
    round(pretty(c(0, max_abs_pp), n = 3), 3)
  )
  size_breaks <- size_breaks[
    size_breaks >= 0 & size_breaks <= max_abs_pp
  ]
}

figure <- ggplot() +
  geom_sf(
    data = mexico,
    fill = "grey98",
    colour = "grey40",
    linewidth = 0.35
  ) +
  geom_sf(
    data = cells_all,
    fill = "grey92",
    colour = "white",
    linewidth = 0.12
  ) +
  geom_sf(
    data = map_cells,
    aes(fill = ballooning_proportion_primary),
    colour = "white",
    linewidth = 0.13
  ) +
  geom_hline(
    yintercept = c(24, 26, 28, 30),
    linetype = "22",
    linewidth = 0.25,
    colour = "grey45"
  ) +
  geom_sf(
    data = map_centroids,
    aes(
      colour = spatialplus_contribution_percentage_points,
      size = absolute_spatialplus_contribution_percentage_points
    ),
    shape = 1,
    stroke = 1.05,
    alpha = 0.95
  ) +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 1),
    oob = scales::squish,
    labels = scales::label_percent(accuracy = 1),
    name = "Observed ballooning-\ncapable genera"
  ) +
  scale_colour_gradient2(
    low = "#2166AC",
    mid = "grey80",
    high = "#B2182B",
    midpoint = 0,
    labels = scales::label_number(
      accuracy = 0.01,
      suffix = " pp",
      show_plus = TRUE
    ),
    name = "Signed Spatial+\ncontribution"
  ) +
  scale_size_continuous(
    range = c(0.45, 4.2),
    breaks = size_breaks,
    labels = scales::label_number(
      accuracy = 0.01,
      suffix = " pp"
    ),
    name = "|Spatial+ contribution|"
  ) +
  coord_sf(
    xlim = c(-117.4, -108.4),
    ylim = c(22.6, 32.8),
    expand = FALSE,
    datum = NA
  ) +
  guides(
    fill = guide_colourbar(
      order = 1,
      barheight = grid::unit(38, "mm")
    ),
    colour = guide_colourbar(
      order = 2,
      barheight = grid::unit(38, "mm")
    ),
    size = guide_legend(
      order = 3,
      override.aes = list(colour = "grey30")
    )
  ) +
  labs(
    x = NULL,
    y = NULL,
    caption = paste0(
      "Cell fill: observed proportion of classified genera capable of ",
      "ballooning. Rings: change in modeled ballooning proportion attributable ",
      "to local VPD and wind-seasonality departures after broad spatial trends ",
      "were removed (percentage points, pp). Grey occupied cells were not in ",
      "the complete current Spatial+ dataset (n = ", nrow(dat), ")."
    )
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.major = element_line(
      colour = "grey88",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 9.2),
    legend.text = element_text(size = 8.3),
    plot.caption = element_text(
      size = 8.4,
      hjust = 0,
      margin = margin(t = 8)
    ),
    plot.margin = margin(7, 7, 7, 7)
  )

figure_base <- file.path(
  figure_dir,
  "Figure_4_observed_ballooning_and_SpatialPlus_contribution"
)

ggsave(
  paste0(figure_base, ".png"),
  figure,
  width = 7.5,
  height = 9.2,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  paste0(figure_base, ".tif"),
  figure,
  width = 7.5,
  height = 9.2,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)
ggsave(
  paste0(figure_base, ".pdf"),
  figure,
  width = 7.5,
  height = 9.2,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)
ggsave(
  paste0(figure_base, ".svg"),
  figure,
  width = 7.5,
  height = 9.2,
  units = "in",
  bg = "white"
)

# ------------------------------------------------------------------------------
# Export mapped values, model summary, and caption
# ------------------------------------------------------------------------------

mapped_values <- map_cells |>
  st_drop_geometry() |>
  select(
    grid_cell_id,
    latitude_band,
    easting_km,
    northing_km,
    ballooning_genera_primary,
    non_ballooning_genera_primary,
    classified_genera_primary,
    ballooning_proportion_primary,
    vpd_z,
    wind_seasonality_z,
    vpd_spplus,
    wind_seasonality_spplus,
    spatialplus_contribution_logit,
    spatialplus_contribution_probability,
    spatialplus_contribution_percentage_points,
    predicted_probability_full,
    predicted_probability_without_local_environment
  )

write_csv(
  mapped_values,
  file.path(output_dir, "12K1_mapped_cell_values.csv")
)

effect_summary <- tibble(
  quantity = c(
    "number_of_mapped_cells",
    "beta_vpd_spplus",
    "beta_wind_seasonality_spplus",
    "minimum_combined_contribution_percentage_points",
    "maximum_combined_contribution_percentage_points",
    "median_absolute_contribution_percentage_points",
    "maximum_absolute_contribution_percentage_points",
    "vpd_spatial_r_squared",
    "wind_seasonality_spatial_r_squared"
  ),
  value = c(
    nrow(map_cells),
    beta_vpd,
    beta_wind,
    min(dat$spatialplus_contribution_percentage_points, na.rm = TRUE),
    max(dat$spatialplus_contribution_percentage_points, na.rm = TRUE),
    median(
      dat$absolute_spatialplus_contribution_percentage_points,
      na.rm = TRUE
    ),
    max(
      dat$absolute_spatialplus_contribution_percentage_points,
      na.rm = TRUE
    ),
    summary(vpd_spatial_model)$r.sq,
    summary(wind_spatial_model)$r.sq
  )
)

write_csv(
  effect_summary,
  file.path(output_dir, "12K1_map_effect_summary.csv")
)

caption <- paste0(
  "Figure 4. Observed ballooning composition and Spatial+ environmental ",
  "contributions across occupied 25-km cells of the Baja California Peninsula. ",
  "Cell fill shows the observed proportion of classified genera considered ",
  "capable of ballooning. Centroid rings show the combined contribution of ",
  "local vapor-pressure-deficit and wind-seasonality departures after broad ",
  "spatial trends were removed. Ring colour indicates whether the combined ",
  "Spatial+ terms increased or decreased the fitted ballooning proportion, ",
  "and ring size indicates the absolute magnitude of that change in percentage ",
  "points relative to the same model with the local environmental terms set ",
  "to zero. Grey occupied cells lacked the complete environmental information ",
  "required for the complete Spatial+ analysis (n = ", nrow(dat), "). The mapped contributions are ",
  "within-sample community-composition associations and should not be ",
  "interpreted as habitat suitability, individual ballooning behavior, or ",
  "causal environmental effects."
)

writeLines(
  caption,
  file.path(output_dir, "Figure_4_publication_caption.txt")
)

writeLines(
  c(
    "STEP 12K1 FIGURE DECISION NOTE",
    "================================",
    "",
    "This figure combines two distinct quantities in one map:",
    "  - polygon fill: observed ballooning-capable genus proportion;",
    "  - centroid rings: combined local Spatial+ model contribution.",
    "",
    "The rings are derived from the primary Spatial+ core model and are shown",
    "on the probability scale as percentage-point changes. The actual values",
    "are exported in 12K1_mapped_cell_values.csv.",
    "",
    "Do not rename this figure as a ballooning-suitability, refugia, or",
    "peninsula-wide prediction map."
  ),
  file.path(output_dir, "README_12K1_OUTPUTS.txt")
)

saveRDS(
  list(
    vpd_spatial_model = vpd_spatial_model,
    wind_spatial_model = wind_spatial_model,
    spatial_plus_model = spatial_plus_model,
    mapped_data = mapped_values
  ),
  file.path(output_dir, "12K1_map_model_objects.rds")
)

cat("\nSpatial+ mapped contribution range (percentage points):\n")
cat(
  min(dat$spatialplus_contribution_percentage_points, na.rm = TRUE),
  "to",
  max(dat$spatialplus_contribution_percentage_points, na.rm = TRUE),
  "\n"
)
cat("Figure written to:", figure_dir, "\n")
cat("STEP 12K1 COMPLETED SUCCESSFULLY\n")

session_info <- capture.output(sessionInfo())
writeLines(
  session_info,
  file.path(output_dir, "12K1_session_info.txt")
)
