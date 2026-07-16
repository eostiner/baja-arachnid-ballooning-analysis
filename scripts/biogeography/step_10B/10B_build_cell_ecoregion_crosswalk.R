#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) path.expand(args[[1]]) else path.expand("~/Desktop/Baja_Ballooning_Pipeline")
seed <- if (length(args) >= 2) as.integer(args[[2]]) else 20260715L
set.seed(seed)

required <- c("sf", "ggplot2", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Missing required R packages: ", paste(missing, collapse = ", "),
    "\nInstall them before rerunning, for example:\n",
    "install.packages(c(", paste(sprintf("'%s'", missing), collapse = ", "), "))"
  )
}

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
})

msg <- function(...) cat(sprintf(...), "\n")
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
norm_name <- function(x) tolower(gsub("[^a-z0-9]+", "", x))
clean_chr <- function(x) trimws(as.character(x))

pipeline_root <- file.path(project_root, "04_analysis", "C3_pipeline_rebuild")
step10_root <- file.path(pipeline_root, "09_C3_biogeographic_concordance")
step10a_dir <- file.path(step10_root, "10A_ecoregion_gis_audit")
out_dir <- file.path(step10_root, "10B_cell_ecoregion_crosswalk")
pub_dir <- file.path(out_dir, "publication_outputs")
dir.create(pub_dir, recursive = TRUE, showWarnings = FALSE)

baja_laea <- "+proj=laea +lat_0=27.5 +lon_0=-114 +datum=WGS84 +units=m +no_defs"
cell_side_m <- 25000
cell_area_km2_expected <- (cell_side_m / 1000)^2

# ---------- Required Step 10A input ----------
eco_path <- file.path(step10a_dir, "10A_ecoregions_validated_mainland_only.gpkg")
if (!file.exists(eco_path)) {
  stop("Missing Step 10A mainland ecoregion file:\n", eco_path,
       "\nRun and review Step 10A before Step 10B.")
}
eco <- st_read(eco_path, quiet = TRUE, stringsAsFactors = FALSE)
if (!"ecoregion_label" %in% names(eco)) stop("Step 10A ecoregion file lacks ecoregion_label.")
if (is.na(st_crs(eco))) stop("Step 10A ecoregion CRS is missing; refusing to infer it.")
eco <- st_make_valid(eco)
eco <- eco[!st_is_empty(eco), "ecoregion_label", drop = FALSE]
eco$ecoregion_label <- clean_chr(eco$ecoregion_label)
eco_laea <- st_transform(eco, baja_laea)

# ---------- Find a one-row-per-occupied-cell coordinate table ----------
priority_csvs <- c(
  file.path(pipeline_root, "06_C3_environment_table", "07_C3_cell_environment_table.csv"),
  file.path(pipeline_root, "05_C3_geographic_GLM", "06_cell_level_model_table_common_C1_C4.csv"),
  file.path(pipeline_root, "01_trait_merge", "C1_C4_cell_trait_counts.csv")
)
all_csvs <- list.files(pipeline_root, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
candidate_csvs <- unique(c(priority_csvs[file.exists(priority_csvs)], all_csvs))

pick_named_column <- function(df, exact, regex = NULL) {
  nms <- names(df); nn <- norm_name(nms)
  exn <- norm_name(exact)
  hit <- match(exn, nn, nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit)) return(nms[hit[[1]]])
  if (!is.null(regex)) {
    idx <- grep(regex, nms, ignore.case = TRUE)
    if (length(idx)) return(nms[idx[[1]]])
  }
  NA_character_
}

score_coordinate_column <- function(df, kind = c("lon", "lat")) {
  kind <- match.arg(kind)
  nms <- names(df)
  scores <- rep(-Inf, length(nms))
  for (i in seq_along(nms)) {
    x <- suppressWarnings(as.numeric(df[[i]]))
    good <- is.finite(x)
    if (sum(good) < max(10, floor(0.7 * nrow(df)))) next
    in_range <- if (kind == "lon") mean(x[good] >= -120.5 & x[good] <= -108.0) else mean(x[good] >= 22.0 & x[good] <= 33.5)
    if (in_range < 0.90) next
    nm <- tolower(nms[[i]])
    score <- 10 * in_range
    if (kind == "lon" && grepl("lon|longitude", nm)) score <- score + 20
    if (kind == "lat" && grepl("lat|latitude", nm)) score <- score + 20
    if (grepl("cent|mid|cell|grid", nm)) score <- score + 5
    if (kind == "lon" && grepl("^x$|x_wgs|xcoord", nm)) score <- score + 2
    if (kind == "lat" && grepl("^y$|y_wgs|ycoord", nm)) score <- score + 2
    scores[[i]] <- score
  }
  if (all(!is.finite(scores))) return(NA_character_)
  nms[[which.max(scores)]]
}

pick_id_column <- function(df) {
  exact <- c("grid25km_id", "cell_id", "grid_id", "gridcell_id", "cellid", "site_id")
  hit <- pick_named_column(df, exact)
  if (!is.na(hit)) return(hit)
  idx <- grep("cell|grid", names(df), ignore.case = TRUE)
  if (length(idx)) {
    unq <- vapply(idx, function(i) length(unique(clean_chr(df[[i]])[!is.na(df[[i]])])), integer(1))
    ok <- idx[unq >= 0.80 * nrow(df)]
    if (length(ok)) return(names(df)[ok[[1]]])
  }
  NA_character_
}

cell_table_candidates <- list()
for (f in candidate_csvs) {
  info <- file.info(f)
  if (!is.finite(info$size) || info$size <= 0 || info$size > 2e8) next
  d <- try(read.csv(f, check.names = FALSE, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(d, "try-error") || nrow(d) < 50 || nrow(d) > 5000) next
  lon_col <- score_coordinate_column(d, "lon")
  lat_col <- score_coordinate_column(d, "lat")
  id_col <- pick_id_column(d)
  if (is.na(lon_col) || is.na(lat_col) || lon_col == lat_col) next
  n_unique_points <- length(unique(paste(round(as.numeric(d[[lon_col]]), 7), round(as.numeric(d[[lat_col]]), 7))))
  n_unique_id <- if (!is.na(id_col)) length(unique(clean_chr(d[[id_col]]))) else n_unique_points
  one_row_score <- -abs(nrow(d) - 205) / 10
  if (n_unique_id == nrow(d)) one_row_score <- one_row_score + 30
  if (grepl("07_C3_cell_environment_table", basename(f), fixed = TRUE)) one_row_score <- one_row_score + 100
  if (grepl("cell_level_model_table", basename(f), ignore.case = TRUE)) one_row_score <- one_row_score + 50
  if (grepl("cell_trait_counts", basename(f), ignore.case = TRUE)) one_row_score <- one_row_score + 40
  cell_table_candidates[[length(cell_table_candidates) + 1L]] <- list(
    file = f, data = d, lon_col = lon_col, lat_col = lat_col, id_col = id_col,
    score = one_row_score, nrow = nrow(d), n_unique_id = n_unique_id
  )
}
if (!length(cell_table_candidates)) {
  stop("Could not find a one-row-per-cell CSV with plausible Baja longitude and latitude columns under:\n", pipeline_root)
}
sel_i <- which.max(vapply(cell_table_candidates, `[[`, numeric(1), "score"))
sel <- cell_table_candidates[[sel_i]]
cell_df <- sel$data
lon_col <- sel$lon_col
lat_col <- sel$lat_col
id_col <- sel$id_col
if (is.na(id_col)) {
  id_col <- "step10b_generated_cell_id"
  cell_df[[id_col]] <- sprintf("cell_%03d", seq_len(nrow(cell_df)))
}

cell_df[[id_col]] <- clean_chr(cell_df[[id_col]])
cell_df[[lon_col]] <- suppressWarnings(as.numeric(cell_df[[lon_col]]))
cell_df[[lat_col]] <- suppressWarnings(as.numeric(cell_df[[lat_col]]))
valid_coord <- is.finite(cell_df[[lon_col]]) & is.finite(cell_df[[lat_col]]) &
  cell_df[[lon_col]] >= -120.5 & cell_df[[lon_col]] <= -108 &
  cell_df[[lat_col]] >= 22 & cell_df[[lat_col]] <= 33.5
cell_df <- cell_df[valid_coord, , drop = FALSE]
if (anyDuplicated(cell_df[[id_col]])) {
  # Keep one row only when duplicates carry the same coordinate pair.
  dup_check <- aggregate(
    paste(round(cell_df[[lon_col]], 7), round(cell_df[[lat_col]], 7)),
    by = list(cell_id = cell_df[[id_col]]), FUN = function(z) length(unique(z))
  )
  if (any(dup_check$x > 1)) stop("Selected cell table has duplicated cell IDs with conflicting coordinates.")
  cell_df <- cell_df[!duplicated(cell_df[[id_col]]), , drop = FALSE]
}
if (nrow(cell_df) < 150 || nrow(cell_df) > 500) {
  stop("Selected cell table yielded ", nrow(cell_df), " unique valid cells; expected roughly 205. Review input selection.")
}

points <- st_as_sf(cell_df, coords = c(lon_col, lat_col), crs = 4326, remove = FALSE)
points_laea <- st_transform(points, baja_laea)

# ---------- Prefer an existing polygon grid if it can be matched safely ----------
find_spatial_candidates <- function(root) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
  files[grepl("\\.(gpkg|geojson|shp)$", files, ignore.case = TRUE) & grepl("25|grid|cell", basename(files), ignore.case = TRUE)]
}
spatial_candidates <- find_spatial_candidates(project_root)
spatial_candidates <- spatial_candidates[!grepl("10A_ecoregions|10B_", basename(spatial_candidates), ignore.case = TRUE)]
cell_polys <- NULL
cell_geometry_method <- NA_character_
matched_spatial_file <- NA_character_

for (f in spatial_candidates) {
  g <- try(st_read(f, quiet = TRUE, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(g, "try-error") || !nrow(g) || is.na(st_crs(g))) next
  gt <- unique(as.character(st_geometry_type(g, by_geometry = TRUE)))
  if (!any(grepl("POLYGON", gt))) next
  bbox <- try(st_bbox(st_transform(g, 4326)), silent = TRUE)
  if (inherits(bbox, "try-error") || bbox[["xmin"]] > -108 || bbox[["xmax"]] < -120.5 || bbox[["ymin"]] > 33.5 || bbox[["ymax"]] < 22) next
  gid <- pick_id_column(st_drop_geometry(g))
  if (!is.na(gid)) {
    g_ids <- clean_chr(g[[gid]])
    m <- match(cell_df[[id_col]], g_ids)
    if (mean(!is.na(m)) >= 0.95) {
      keep <- !is.na(m)
      gsub <- g[m[keep], , drop = FALSE]
      gsub$step10b_cell_id <- cell_df[[id_col]][keep]
      # Add any unmatched occupied cells later by reconstruction only if <=5%.
      cell_polys <- st_transform(gsub[, c("step10b_cell_id", "geometry")], baja_laea)
      cell_geometry_method <- "existing_project_polygon_grid_matched_by_id"
      matched_spatial_file <- f
      if (sum(!keep)) {
        xy <- st_coordinates(points_laea[!keep, ])
        half <- cell_side_m / 2
        ps <- lapply(seq_len(nrow(xy)), function(i) {
          x <- xy[i, 1]; y <- xy[i, 2]
          st_polygon(list(matrix(c(
            x-half, y-half, x+half, y-half, x+half, y+half,
            x-half, y+half, x-half, y-half
          ), ncol = 2, byrow = TRUE)))
        })
        add <- st_sf(step10b_cell_id = cell_df[[id_col]][!keep], geometry = st_sfc(ps, crs = baja_laea))
        cell_polys <- rbind(cell_polys, add)
        cell_geometry_method <- paste0(cell_geometry_method, ";", sum(!keep), "_unmatched_cells_reconstructed")
      }
      break
    }
  }
}

if (is.null(cell_polys)) {
  xy <- st_coordinates(points_laea)
  half <- cell_side_m / 2
  ps <- lapply(seq_len(nrow(xy)), function(i) {
    x <- xy[i, 1]; y <- xy[i, 2]
    st_polygon(list(matrix(c(
      x-half, y-half, x+half, y-half, x+half, y+half,
      x-half, y+half, x-half, y-half
    ), ncol = 2, byrow = TRUE)))
  })
  cell_polys <- st_sf(step10b_cell_id = cell_df[[id_col]], geometry = st_sfc(ps, crs = baja_laea))
  cell_geometry_method <- "25km_square_reconstructed_from_cell_centroids_in_baja_equal_area_projection"
}

cell_polys <- st_make_valid(cell_polys)
cell_polys <- cell_polys[match(cell_df[[id_col]], cell_polys$step10b_cell_id), ]
if (any(is.na(cell_polys$step10b_cell_id))) stop("Cell polygon construction failed to preserve all selected cell IDs.")

# Geometry audit.
cell_area_km2 <- as.numeric(st_area(cell_polys)) / 1e6
area_deviation_pct <- 100 * (cell_area_km2 - cell_area_km2_expected) / cell_area_km2_expected
coords <- st_coordinates(points_laea)
nearest_m <- rep(NA_real_, nrow(points_laea))
if (nrow(points_laea) > 1) {
  dd <- as.matrix(st_distance(points_laea))
  diag(dd) <- Inf
  nearest_m <- apply(dd, 1, min)
}

# ---------- Dominant-area assignment ----------
cell_polys$step10b_row <- seq_len(nrow(cell_polys))
inter <- suppressWarnings(st_intersection(
  cell_polys[, c("step10b_cell_id", "step10b_row")],
  eco_laea[, "ecoregion_label"]
))
if (!nrow(inter)) stop("No occupied cell polygons intersected the validated mainland ecoregions.")
inter$overlap_km2 <- as.numeric(st_area(inter)) / 1e6
inter_df <- st_drop_geometry(inter)
inter_df <- inter_df[is.finite(inter_df$overlap_km2) & inter_df$overlap_km2 > 1e-8, ]

cross <- vector("list", nrow(cell_polys))
for (i in seq_len(nrow(cell_polys))) {
  cid <- cell_polys$step10b_cell_id[[i]]
  z <- inter_df[inter_df$step10b_cell_id == cid, c("ecoregion_label", "overlap_km2"), drop = FALSE]
  if (nrow(z)) {
    z <- aggregate(overlap_km2 ~ ecoregion_label, data = z, FUN = sum)
    z <- z[order(-z$overlap_km2, z$ecoregion_label), , drop = FALSE]
  }
  full_area <- cell_area_km2[[i]]
  covered <- if (nrow(z)) sum(z$overlap_km2) else 0
  if (!nrow(z) || covered <= 0) {
    cross[[i]] <- data.frame(
      step10b_cell_id = cid,
      dominant_ecoregion = NA_character_, second_ecoregion = NA_character_,
      dominant_overlap_km2 = 0, second_overlap_km2 = 0,
      ecoregion_covered_km2 = 0, full_cell_area_km2 = full_area,
      dominant_fraction_of_covered_land = NA_real_, second_fraction_of_covered_land = NA_real_,
      dominant_minus_second_fraction = NA_real_, dominant_fraction_of_full_cell = 0,
      ecoregion_coverage_fraction = 0
    )
  } else {
    dom <- z[1, ]
    sec <- if (nrow(z) >= 2) z[2, ] else data.frame(ecoregion_label = NA_character_, overlap_km2 = 0)
    cross[[i]] <- data.frame(
      step10b_cell_id = cid,
      dominant_ecoregion = dom$ecoregion_label,
      second_ecoregion = sec$ecoregion_label,
      dominant_overlap_km2 = dom$overlap_km2,
      second_overlap_km2 = sec$overlap_km2,
      ecoregion_covered_km2 = covered,
      full_cell_area_km2 = full_area,
      dominant_fraction_of_covered_land = dom$overlap_km2 / covered,
      second_fraction_of_covered_land = sec$overlap_km2 / covered,
      dominant_minus_second_fraction = (dom$overlap_km2 - sec$overlap_km2) / covered,
      dominant_fraction_of_full_cell = dom$overlap_km2 / full_area,
      ecoregion_coverage_fraction = covered / full_area
    )
  }
}
cross_df <- do.call(rbind, cross)

# ---------- Centroid assignment and agreement ----------
hits <- st_intersects(points_laea, eco_laea)
centroid_label <- rep(NA_character_, length(hits))
centroid_on_boundary <- lengths(hits) > 1L
for (i in seq_along(hits)) {
  h <- hits[[i]]
  if (length(h) == 1L) centroid_label[[i]] <- eco_laea$ecoregion_label[h]
  if (length(h) > 1L) {
    # Boundary points are assigned to dominant area, while retaining a flag.
    centroid_label[[i]] <- cross_df$dominant_ecoregion[match(cell_df[[id_col]][i], cross_df$step10b_cell_id)]
  }
}
centroid_outside <- lengths(hits) == 0L

cross_df$centroid_ecoregion <- centroid_label[match(cross_df$step10b_cell_id, cell_df[[id_col]])]
cross_df$centroid_on_ecoregion_boundary <- centroid_on_boundary[match(cross_df$step10b_cell_id, cell_df[[id_col]])]
cross_df$centroid_outside_mapped_mainland <- centroid_outside[match(cross_df$step10b_cell_id, cell_df[[id_col]])]
cross_df$centroid_agrees_with_dominant <- !is.na(cross_df$centroid_ecoregion) & cross_df$centroid_ecoregion == cross_df$dominant_ecoregion
cross_df$ambiguous_dominant_assignment <- is.na(cross_df$dominant_ecoregion) |
  cross_df$dominant_fraction_of_covered_land < 0.50 |
  cross_df$dominant_minus_second_fraction < 0.10
cross_df$low_mapped_land_coverage <- cross_df$ecoregion_coverage_fraction < 0.10
cross_df$primary_assignment_eligible <- !is.na(cross_df$dominant_ecoregion) & !cross_df$low_mapped_land_coverage
cross_df$sensitivity_unambiguous_eligible <- cross_df$primary_assignment_eligible & !cross_df$ambiguous_dominant_assignment

# Merge selected biological cell table attributes after renaming the key consistently.
cell_attr <- cell_df
names(cell_attr)[names(cell_attr) == id_col] <- "step10b_cell_id"
# Avoid duplicate coordinate names in output by preserving source names as-is.
out <- merge(cross_df, cell_attr, by = "step10b_cell_id", all.x = TRUE, sort = FALSE)
out <- out[match(cell_polys$step10b_cell_id, out$step10b_cell_id), ]

# ---------- Summaries ----------
count_primary <- aggregate(
  step10b_cell_id ~ dominant_ecoregion,
  data = out[out$primary_assignment_eligible, ], FUN = length
)
names(count_primary)[2] <- "n_cells_primary"
count_unamb <- aggregate(
  step10b_cell_id ~ dominant_ecoregion,
  data = out[out$sensitivity_unambiguous_eligible, ], FUN = length
)
names(count_unamb)[2] <- "n_cells_unambiguous_sensitivity"
eco_counts <- merge(
  data.frame(dominant_ecoregion = sort(unique(eco$ecoregion_label))),
  count_primary, by = "dominant_ecoregion", all.x = TRUE
)
eco_counts <- merge(eco_counts, count_unamb, by = "dominant_ecoregion", all.x = TRUE)
eco_counts[is.na(eco_counts)] <- 0

summary_df <- data.frame(
  metric = c(
    "selected_cell_table", "cell_id_column", "longitude_column", "latitude_column",
    "occupied_cells", "cell_geometry_method", "matched_spatial_grid_file",
    "median_cell_area_km2", "max_abs_cell_area_deviation_percent",
    "median_nearest_centroid_distance_km", "ecoregions_mainland",
    "cells_with_dominant_assignment", "ambiguous_dominant_assignments",
    "centroids_outside_mapped_mainland", "centroid_dominant_disagreements",
    "cells_with_low_mapped_land_coverage", "unambiguous_sensitivity_cells"
  ),
  value = c(
    normalizePath(sel$file, winslash = "/", mustWork = FALSE), id_col, lon_col, lat_col,
    nrow(out), cell_geometry_method, matched_spatial_file,
    sprintf("%.6f", median(cell_area_km2, na.rm = TRUE)), sprintf("%.6f", max(abs(area_deviation_pct), na.rm = TRUE)),
    sprintf("%.6f", median(nearest_m, na.rm = TRUE) / 1000), nrow(eco),
    sum(out$primary_assignment_eligible), sum(out$ambiguous_dominant_assignment, na.rm = TRUE),
    sum(out$centroid_outside_mapped_mainland, na.rm = TRUE),
    sum(!out$centroid_agrees_with_dominant & !out$centroid_outside_mapped_mainland, na.rm = TRUE),
    sum(out$low_mapped_land_coverage, na.rm = TRUE), sum(out$sensitivity_unambiguous_eligible, na.rm = TRUE)
  )
)

input_manifest <- data.frame(
  role = c("validated_mainland_ecoregions", "occupied_cell_coordinate_table", "matched_existing_grid_if_any"),
  file = c(eco_path, sel$file, matched_spatial_file),
  md5 = c(
    unname(tools::md5sum(eco_path)),
    unname(tools::md5sum(sel$file)),
    if (!is.na(matched_spatial_file) && file.exists(matched_spatial_file)) unname(tools::md5sum(matched_spatial_file)) else NA_character_
  ),
  stringsAsFactors = FALSE
)

# ---------- Write tabular and spatial outputs ----------
write.csv(out, file.path(out_dir, "10B_cell_ecoregion_crosswalk.csv"), row.names = FALSE, na = "")
write.csv(eco_counts, file.path(out_dir, "10B_ecoregion_cell_counts.csv"), row.names = FALSE, na = "")
write.csv(summary_df, file.path(out_dir, "10B_assignment_audit_summary.csv"), row.names = FALSE, na = "")
write.csv(input_manifest, file.path(out_dir, "10B_input_manifest.csv"), row.names = FALSE, na = "")
write.csv(inter_df[, c("step10b_cell_id", "ecoregion_label", "overlap_km2")],
          file.path(out_dir, "10B_all_cell_ecoregion_overlaps.csv"), row.names = FALSE, na = "")

cell_out <- merge(cell_polys, out, by = "step10b_cell_id", all.x = TRUE, sort = FALSE)
cell_out <- cell_out[match(cell_polys$step10b_cell_id, cell_out$step10b_cell_id), ]
st_write(st_transform(cell_out, 4326), file.path(out_dir, "10B_occupied_cells_with_ecoregion_assignments.gpkg"), delete_dsn = TRUE, quiet = TRUE)

# ---------- Audit map ----------
map_cells <- st_transform(cell_out, 4326)
map_points <- st_transform(points_laea, 4326)
map_points$ambiguous <- out$ambiguous_dominant_assignment[match(cell_df[[id_col]], out$step10b_cell_id)]

p <- ggplot() +
  geom_sf(data = eco, aes(fill = ecoregion_label), color = "grey35", linewidth = 0.18, alpha = 0.75) +
  geom_sf(data = map_cells, fill = NA, color = "black", linewidth = 0.16, alpha = 0.6) +
  geom_sf(data = map_points[map_points$ambiguous %in% TRUE, ], shape = 4, size = 2.0, stroke = 0.7, color = "black") +
  coord_sf(expand = FALSE) +
  guides(fill = guide_legend(title = "Dominant mainland ecoregion", ncol = 1)) +
  labs(
    title = "Step 10B: occupied 25-km cell to ecoregion crosswalk",
    subtitle = "Cell outlines; X marks assignments flagged as boundary-ambiguous",
    x = NULL, y = NULL,
    caption = paste0("Primary assignment: largest area overlap. Cell geometry: ", cell_geometry_method)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 8),
    plot.caption = element_text(size = 7, hjust = 0)
  )

png_path <- file.path(pub_dir, "10B_cell_ecoregion_crosswalk_audit_map.png")
pdf_path <- file.path(pub_dir, "10B_cell_ecoregion_crosswalk_audit_map.pdf")
ggsave(png_path, p, width = 10.5, height = 8.5, units = "in", dpi = 300, bg = "white")
ggsave(pdf_path, p, width = 10.5, height = 8.5, units = "in", device = "pdf")

status <- "PASS"
if (sum(out$primary_assignment_eligible) != nrow(out)) status <- "PASS_WITH_INELIGIBLE_CELLS_REQUIRING_REVIEW"
if (sum(out$ambiguous_dominant_assignment, na.rm = TRUE) > 0) status <- paste0(status, "_AND_AMBIGUOUS_BOUNDARY_CELLS")

summary_json <- list(
  step = "10B",
  audit_status = status,
  seed = seed,
  project_root = project_root,
  output_dir = out_dir,
  selected_cell_table = normalizePath(sel$file, winslash = "/", mustWork = FALSE),
  cell_id_column = id_col,
  longitude_column = lon_col,
  latitude_column = lat_col,
  occupied_cells = nrow(out),
  cell_geometry_method = cell_geometry_method,
  matched_existing_grid_file = matched_spatial_file,
  mainland_ecoregions = nrow(eco),
  primary_assignment_eligible_cells = sum(out$primary_assignment_eligible),
  ambiguous_cells = sum(out$ambiguous_dominant_assignment, na.rm = TRUE),
  unambiguous_sensitivity_cells = sum(out$sensitivity_unambiguous_eligible, na.rm = TRUE),
  centroid_outside_cells = sum(out$centroid_outside_mapped_mainland, na.rm = TRUE),
  centroid_dominant_disagreements = sum(!out$centroid_agrees_with_dominant & !out$centroid_outside_mapped_mainland, na.rm = TRUE),
  next_step = "Review the crosswalk map and ambiguity table before equal-cell richness and neighboring-cell turnover tests."
)
jsonlite::write_json(summary_json, file.path(out_dir, "10B_audit_summary.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")

readme <- c(
  "STEP 10B — OCCUPIED 25-KM CELL TO ECOREGION CROSSWALK",
  "",
  paste0("AUDIT_STATUS=", status),
  paste0("OUTPUT_DIR=", out_dir),
  paste0("SELECTED_CELL_TABLE=", normalizePath(sel$file, winslash = "/", mustWork = FALSE)),
  paste0("CELL_ID_COLUMN=", id_col),
  paste0("LONGITUDE_COLUMN=", lon_col),
  paste0("LATITUDE_COLUMN=", lat_col),
  paste0("OCCUPIED_CELLS=", nrow(out)),
  paste0("CELL_GEOMETRY_METHOD=", cell_geometry_method),
  paste0("PRIMARY_ASSIGNMENT_ELIGIBLE=", sum(out$primary_assignment_eligible)),
  paste0("AMBIGUOUS_ASSIGNMENTS=", sum(out$ambiguous_dominant_assignment, na.rm = TRUE)),
  paste0("CENTROID_OUTSIDE_MAINLAND_MAP=", sum(out$centroid_outside_mapped_mainland, na.rm = TRUE)),
  paste0("CENTROID_DOMINANT_DISAGREEMENTS=", sum(!out$centroid_agrees_with_dominant & !out$centroid_outside_mapped_mainland, na.rm = TRUE)),
  paste0("UNAMBIGUOUS_SENSITIVITY_CELLS=", sum(out$sensitivity_unambiguous_eligible, na.rm = TRUE)),
  "",
  "PRIMARY RULE:",
  "- Assign each occupied cell to the independently mapped mainland ecoregion with the largest area overlap.",
  "- Fractions are calculated relative to the cell area covered by mapped mainland ecoregions, so ocean area does not make coastal cells artificially ambiguous.",
  "",
  "AMBIGUITY RULE:",
  "- Flag dominant overlap < 0.50 of mapped land, or dominant-minus-second overlap < 0.10.",
  "- Also record centroid assignment, centroid disagreement, and low mapped-land coverage (<0.10 of the full cell).",
  "",
  "REVIEW BEFORE STEP 10C/10D:",
  "1. publication_outputs/10B_cell_ecoregion_crosswalk_audit_map.png",
  "2. 10B_assignment_audit_summary.csv",
  "3. 10B_ecoregion_cell_counts.csv",
  "4. 10B_cell_ecoregion_crosswalk.csv",
  "",
  "This step creates the spatial crosswalk only. It does not test richness or turnover."
)
writeLines(readme, file.path(out_dir, "README_RESULTS_FIRST.txt"), useBytes = TRUE)

msg("STEP 10B COMPLETE")
msg("AUDIT_STATUS=%s", status)
msg("OUTPUT_DIR=%s", out_dir)
msg("SELECTED_CELL_TABLE=%s", sel$file)
msg("OCCUPIED_CELLS=%d", nrow(out))
msg("PRIMARY_ASSIGNMENT_ELIGIBLE=%d", sum(out$primary_assignment_eligible))
msg("AMBIGUOUS_ASSIGNMENTS=%d", sum(out$ambiguous_dominant_assignment, na.rm = TRUE))
msg("CENTROID_OUTSIDE_MAINLAND_MAP=%d", sum(out$centroid_outside_mapped_mainland, na.rm = TRUE))
msg("UNAMBIGUOUS_SENSITIVITY_CELLS=%d", sum(out$sensitivity_unambiguous_eligible, na.rm = TRUE))
