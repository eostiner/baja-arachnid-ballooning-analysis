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

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
msg <- function(...) cat(sprintf(...), "\n")
clean_text <- function(x) {
  x <- as.character(x)
  x <- iconv(x, from = "", to = "UTF-8", sub = "byte")
  trimws(x)
}

pipeline_root <- file.path(project_root, "04_analysis", "C3_pipeline_rebuild")
out_dir <- file.path(pipeline_root, "09_C3_biogeographic_concordance", "10A_ecoregion_gis_audit")
source_dir <- file.path(out_dir, "00_source_download")
unzip_dir <- file.path(source_dir, "unzipped")
pub_dir <- file.path(out_dir, "publication_outputs")
dir.create(unzip_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pub_dir, recursive = TRUE, showWarnings = FALSE)

zip_path <- file.path(source_dir, "ecoreg_final_version_2010.zip")
urls <- c(
  "https://ezcurralab.ucr.edu/baja_california_ecoregions/ecoreg_final_version_2010.zip",
  "http://ezcurralab.ucr.edu/baja_california_ecoregions/ecoreg_final_version_2010.zip"
)
arcgis_item_id <- "e09e39d83c3347879e133cabfabbdcb1"
arcgis_item_page <- paste0("https://www.arcgis.com/home/item.html?id=", arcgis_item_id)
arcgis_rest_base <- paste0("https://www.arcgis.com/sharing/rest/content/items/", arcgis_item_id)
source_method <- NA_character_
source_reference <- NA_character_
source_files <- character()

download_with_fallback <- function(url, dest, quiet = TRUE) {
  unlink(dest, force = TRUE)
  result <- try(
    suppressWarnings(utils::download.file(url, dest, mode = "wb", method = "libcurl", quiet = quiet)),
    silent = TRUE
  )
  ok <- !inherits(result, "try-error") && file.exists(dest) && file.info(dest)$size > 100
  if (!ok) {
    curl <- Sys.which("curl")
    if (nzchar(curl)) {
      status <- suppressWarnings(system2(curl, c("-L", "--fail", "--retry", "3", "--connect-timeout", "30", "-o", shQuote(dest), shQuote(url))))
      ok <- identical(status, 0L) && file.exists(dest) && file.info(dest)$size > 100
    }
  }
  ok
}

fetch_json <- function(url, dest) {
  if (!download_with_fallback(url, dest, quiet = TRUE)) return(NULL)
  parsed <- try(jsonlite::fromJSON(dest, simplifyVector = FALSE), silent = TRUE)
  if (inherits(parsed, "try-error") || !is.list(parsed) || !is.null(parsed$error)) return(NULL)
  parsed
}

extract_service_urls <- function(x) {
  out <- character()
  recurse <- function(z) {
    if (!is.list(z)) return(invisible(NULL))
    if (!is.null(z$url) && is.character(z$url) && length(z$url) == 1L) out <<- c(out, z$url)
    for (el in z) recurse(el)
    invisible(NULL)
  }
  recurse(x)
  unique(out[grepl("/(FeatureServer|MapServer)(/|$)", out, ignore.case = TRUE)])
}

# 1. Prefer the authors' official companion ZIP when it is available.
official_ok <- FALSE
if (file.exists(zip_path) && file.info(zip_path)$size >= 1000) {
  official_ok <- TRUE
} else {
  msg("Downloading the official Gonzalez-Abraham et al. companion GIS...")
  for (u in urls) {
    msg("  trying: %s", u)
    if (download_with_fallback(u, zip_path, quiet = FALSE) && file.info(zip_path)$size >= 1000) {
      official_ok <- TRUE
      break
    }
  }
}

if (official_ok) {
  unzip_result <- try(utils::unzip(zip_path, exdir = unzip_dir, overwrite = TRUE), silent = TRUE)
  if (!inherits(unzip_result, "try-error")) {
    source_files <- list.files(unzip_dir, pattern = "\\.(shp|geojson|gpkg)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  }
  if (length(source_files)) {
    source_method <- "official_UCR_companion_ZIP"
    source_reference <- urls[[1]]
  } else {
    official_ok <- FALSE
    msg("The downloaded UCR file was not a readable spatial archive; trying the ArcGIS fallback.")
  }
}

# 2. The UCR page currently points to a retired path. Fall back to a public
# ArcGIS item carrying the same Gonzalez-Abraham ecoregion map, and discover
# its polygon service dynamically rather than hard-coding a transient server URL.
if (!official_ok) {
  unlink(zip_path, force = TRUE)
  msg("Official UCR archive unavailable. Trying the ArcGIS-hosted Baja Ecoregions item...")
  meta_path <- file.path(source_dir, "arcgis_item_metadata.json")
  data_path <- file.path(source_dir, "arcgis_item_data.json")
  item_meta <- fetch_json(paste0(arcgis_rest_base, "?f=json"), meta_path)
  item_data <- fetch_json(paste0(arcgis_rest_base, "/data?f=json"), data_path)
  if (is.null(item_meta)) stop("Could not read the ArcGIS item metadata after the official UCR URL returned 404.")

  service_urls <- character()
  if (!is.null(item_meta$url) && is.character(item_meta$url)) service_urls <- c(service_urls, item_meta$url)
  if (!is.null(item_data)) service_urls <- c(service_urls, extract_service_urls(item_data))
  service_urls <- unique(service_urls[grepl("/(FeatureServer|MapServer)(/|$)", service_urls, ignore.case = TRUE)])
  if (!length(service_urls)) stop("The ArcGIS fallback item did not expose a FeatureServer or MapServer URL.")

  downloaded_layers <- character()
  layer_manifest <- list()
  manifest_i <- 0L
  for (sidx in seq_along(service_urls)) {
    service_url <- sub("/+$", "", service_urls[[sidx]])
    layer_urls <- character()
    if (grepl("/(FeatureServer|MapServer)/[0-9]+$", service_url, ignore.case = TRUE)) {
      layer_urls <- service_url
    } else {
      service_json_path <- file.path(source_dir, sprintf("arcgis_service_%02d.json", sidx))
      service_meta <- fetch_json(paste0(service_url, "?f=json"), service_json_path)
      if (is.null(service_meta) || is.null(service_meta$layers)) next
      entries <- service_meta$layers
      if (!is.list(entries)) next
      for (entry in entries) {
        if (is.list(entry) && !is.null(entry$id)) layer_urls <- c(layer_urls, paste0(service_url, "/", entry$id))
      }
    }

    for (layer_url in unique(layer_urls)) {
      layer_id <- sub("^.*/", "", layer_url)
      layer_json_path <- file.path(source_dir, sprintf("arcgis_layer_%02d_%s_metadata.json", sidx, layer_id))
      layer_meta <- fetch_json(paste0(layer_url, "?f=json"), layer_json_path)
      if (is.null(layer_meta)) next
      geom_type <- layer_meta$geometryType %||% ""
      if (!grepl("Polygon", geom_type, ignore.case = TRUE)) next

      geojson_path <- file.path(unzip_dir, sprintf("arcgis_baja_ecoregions_service%02d_layer%s.geojson", sidx, layer_id))
      query_url <- paste0(
        layer_url,
        "/query?where=1%3D1&outFields=*&returnGeometry=true&outSR=4326&f=geojson"
      )
      ok <- download_with_fallback(query_url, geojson_path, quiet = TRUE)
      if (!ok) next
      test_obj <- try(st_read(geojson_path, quiet = TRUE, stringsAsFactors = FALSE), silent = TRUE)
      if (inherits(test_obj, "try-error") || !nrow(test_obj) || !any(grepl("POLYGON", as.character(st_geometry_type(test_obj, by_geometry = TRUE))))) {
        unlink(geojson_path, force = TRUE)
        next
      }
      downloaded_layers <- c(downloaded_layers, geojson_path)
      manifest_i <- manifest_i + 1L
      layer_manifest[[manifest_i]] <- data.frame(
        service_url = service_url,
        layer_url = layer_url,
        layer_id = layer_id,
        layer_name = layer_meta$name %||% NA_character_,
        geometry_type = geom_type,
        feature_count_downloaded = nrow(test_obj),
        downloaded_file = normalizePath(geojson_path, winslash = "/", mustWork = FALSE)
      )
    }
  }
  if (!length(downloaded_layers)) stop("The ArcGIS fallback was reachable, but no readable polygon layer could be downloaded.")
  source_files <- unique(downloaded_layers)
  if (length(layer_manifest)) write.csv(do.call(rbind, layer_manifest), file.path(out_dir, "10A_arcgis_fallback_layer_manifest.csv"), row.names = FALSE)
  source_method <- "ArcGIS_public_item_fallback"
  source_reference <- arcgis_item_page
}

spatial_files <- unique(c(
  source_files,
  list.files(unzip_dir, pattern = "\\.(shp|geojson|gpkg)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
))
spatial_files <- spatial_files[file.exists(spatial_files)]
if (!length(spatial_files)) stop("No readable shapefile, GeoJSON, or GeoPackage source files were available.")

source_md5s <- unname(tools::md5sum(spatial_files))
zip_md5 <- if (official_ok && file.exists(zip_path)) unname(tools::md5sum(zip_path)) else paste(source_md5s, collapse = ";")
# Baja-centered equal-area projection for areas and boundary lengths.
baja_laea <- "+proj=laea +lat_0=27.5 +lon_0=-114 +datum=WGS84 +units=m +no_defs"

candidate_field_pattern <- "eco|region|nombre|name|unidad|unit|descrip|class|clase|tipo|type|legend|leyenda"
layer_inventory <- list()
field_inventory <- list()
layer_objects <- list()

msg("Reading and auditing %d spatial layer(s)...", length(spatial_files))
for (i in seq_along(spatial_files)) {
  f <- spatial_files[[i]]
  obj <- try(st_read(f, quiet = TRUE, stringsAsFactors = FALSE, options = "ENCODING=UTF-8"), silent = TRUE)
  if (inherits(obj, "try-error")) obj <- try(st_read(f, quiet = TRUE, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(obj, "try-error")) {
    layer_inventory[[i]] <- data.frame(
      layer_id = i, file = f, readable = FALSE, read_error = as.character(obj),
      n_features = NA, geometry_type = NA, crs_epsg = NA, crs_name = NA,
      bbox_xmin = NA, bbox_ymin = NA, bbox_xmax = NA, bbox_ymax = NA,
      bbox_plausible_baja = FALSE, invalid_geometries = NA, empty_geometries = NA,
      candidate_label_fields = NA, selection_score = -Inf
    )
    next
  }
  layer_objects[[as.character(i)]] <- obj
  geom_type <- paste(sort(unique(as.character(st_geometry_type(obj, by_geometry = TRUE)))), collapse = ";")
  crs <- st_crs(obj)
  bbox_wgs <- try(st_bbox(st_transform(obj, 4326)), silent = TRUE)
  if (inherits(bbox_wgs, "try-error")) bbox_wgs <- c(xmin = NA, ymin = NA, xmax = NA, ymax = NA)
  plausible <- all(is.finite(bbox_wgs)) &&
    bbox_wgs[["xmin"]] < -108 && bbox_wgs[["xmax"]] > -119.5 &&
    bbox_wgs[["ymin"]] < 33.5 && bbox_wgs[["ymax"]] > 22
  candidate_fields <- names(obj)[grepl(candidate_field_pattern, names(obj), ignore.case = TRUE)]
  score <- 0
  is_poly <- any(grepl("POLYGON", geom_type))
  if (is_poly) score <- score + 20
  if (plausible) score <- score + 20
  if (length(candidate_fields)) score <- score + 10
  if (nrow(obj) >= 5 && nrow(obj) <= 100) score <- score + 10
  if (nrow(obj) > 0) score <- score + min(log1p(nrow(obj)), 5)
  layer_inventory[[i]] <- data.frame(
    layer_id = i, file = normalizePath(f, winslash = "/", mustWork = FALSE), readable = TRUE,
    read_error = NA_character_, n_features = nrow(obj), geometry_type = geom_type,
    crs_epsg = crs$epsg %||% NA_integer_, crs_name = crs$Name %||% crs$input %||% NA_character_,
    bbox_xmin = bbox_wgs[["xmin"]], bbox_ymin = bbox_wgs[["ymin"]],
    bbox_xmax = bbox_wgs[["xmax"]], bbox_ymax = bbox_wgs[["ymax"]],
    bbox_plausible_baja = plausible,
    invalid_geometries = sum(!st_is_valid(obj), na.rm = TRUE),
    empty_geometries = sum(st_is_empty(obj), na.rm = TRUE),
    candidate_label_fields = paste(candidate_fields, collapse = ";"),
    selection_score = score
  )
  attrs <- st_drop_geometry(obj)
  if (ncol(attrs)) {
    field_inventory[[length(field_inventory) + 1L]] <- do.call(rbind, lapply(names(attrs), function(nm) {
      vals <- attrs[[nm]]
      data.frame(
        layer_id = i,
        field = nm,
        class = paste(class(vals), collapse = ";"),
        n_unique_nonmissing = length(unique(vals[!is.na(vals)])),
        n_missing = sum(is.na(vals)),
        example_values = paste(head(unique(clean_text(vals[!is.na(vals)])), 8), collapse = " | ")
      )
    }))
  }
}

layer_inventory_df <- do.call(rbind, layer_inventory)
field_inventory_df <- if (length(field_inventory)) do.call(rbind, field_inventory) else data.frame()
write.csv(layer_inventory_df, file.path(out_dir, "10A_layer_inventory.csv"), row.names = FALSE, na = "")
write.csv(field_inventory_df, file.path(out_dir, "10A_field_inventory.csv"), row.names = FALSE, na = "")

eligible <- layer_inventory_df$readable & grepl("POLYGON", layer_inventory_df$geometry_type) & layer_inventory_df$bbox_plausible_baja
if (!any(eligible)) stop("No readable polygon layer had a plausible Baja California extent.")
selected_row <- layer_inventory_df[eligible, ][which.max(layer_inventory_df$selection_score[eligible]), , drop = FALSE]
selected_id <- as.character(selected_row$layer_id[[1]])
eco_raw <- layer_objects[[selected_id]]

if (is.na(st_crs(eco_raw))) stop("Selected ecoregion layer has no CRS. The script will not infer a CRS silently.")

eco <- st_make_valid(eco_raw)
eco <- eco[!st_is_empty(eco), , drop = FALSE]
eco <- st_transform(eco, 4326)
eco$.source_feature_id <- seq_len(nrow(eco))

# Choose a human-readable label field. Prefer likely field names with 3-40 unique values.
attrs <- st_drop_geometry(eco)
label_candidates <- names(attrs)[grepl(candidate_field_pattern, names(attrs), ignore.case = TRUE)]
if (length(label_candidates)) {
  nunique <- vapply(label_candidates, function(nm) length(unique(clean_text(attrs[[nm]][!is.na(attrs[[nm]])]))), integer(1))
  preferred <- label_candidates[nunique >= 3 & nunique <= 40]
  if (!length(preferred)) preferred <- label_candidates
  label_field <- preferred[[which.max(vapply(preferred, function(nm) {
    mean(nchar(clean_text(attrs[[nm]])), na.rm = TRUE)
  }, numeric(1)))]]
} else {
  non_geom <- setdiff(names(eco), attr(eco, "sf_column"))
  usable <- non_geom[vapply(non_geom, function(nm) {
    x <- attrs[[nm]]
    is.character(x) || is.factor(x)
  }, logical(1))]
  label_field <- if (length(usable)) usable[[1]] else ".source_feature_id"
}

eco$.eco_label <- clean_text(eco[[label_field]])
eco$.eco_label[is.na(eco$.eco_label) | !nzchar(eco$.eco_label)] <- paste0("Unlabelled_", eco$.source_feature_id[is.na(eco$.eco_label) | !nzchar(eco$.eco_label)])

# Keep source polygons intact, but also generate a dissolved layer by label.
eco_laea <- st_transform(eco, baja_laea)
area_km2 <- as.numeric(st_area(eco_laea)) / 1e6
geom_text <- st_as_text(st_geometry(eco), digits = 12)
duplicate_geometry <- duplicated(geom_text) | duplicated(geom_text, fromLast = TRUE)

eco_summary <- data.frame(
  source_feature_id = eco$.source_feature_id,
  ecoregion_label = eco$.eco_label,
  source_label_field = label_field,
  geometry_type = as.character(st_geometry_type(eco)),
  valid_after_repair = st_is_valid(eco),
  empty = st_is_empty(eco),
  area_km2 = area_km2,
  duplicate_geometry = duplicate_geometry
)
write.csv(eco_summary, file.path(out_dir, "10A_selected_layer_feature_audit.csv"), row.names = FALSE)

attribute_values <- as.data.frame(table(ecoregion_label = eco$.eco_label), stringsAsFactors = FALSE)
names(attribute_values)[2] <- "source_polygon_count"
attribute_values$area_km2 <- vapply(attribute_values$ecoregion_label, function(lbl) {
  sum(area_km2[eco$.eco_label == lbl])
}, numeric(1))
attribute_values <- attribute_values[order(attribute_values$ecoregion_label), ]
write.csv(attribute_values, file.path(out_dir, "10A_ecoregion_labels_and_areas.csv"), row.names = FALSE)

# Dissolve by label to simplify later cell assignment and prevent multipart labels from acting as separate classes.
split_idx <- split(seq_len(nrow(eco)), eco$.eco_label)
dissolved_list <- lapply(names(split_idx), function(lbl) {
  g <- st_union(st_geometry(eco[split_idx[[lbl]], ]))
  st_sf(ecoregion_label = lbl, source_polygon_count = length(split_idx[[lbl]]), geometry = g)
})
eco_dissolved <- do.call(rbind, dissolved_list)
eco_dissolved <- st_make_valid(eco_dissolved)

# Pairwise overlap audit in an equal-area CRS.
diss_laea <- st_transform(eco_dissolved, baja_laea)
overlap_rows <- list()
if (nrow(diss_laea) > 1) {
  k <- 0L
  for (i in seq_len(nrow(diss_laea) - 1L)) {
    for (j in (i + 1L):nrow(diss_laea)) {
      if (length(st_intersects(diss_laea[i, ], diss_laea[j, ], sparse = TRUE)[[1]]) > 0) {
        inter <- suppressWarnings(st_intersection(st_geometry(diss_laea[i, ]), st_geometry(diss_laea[j, ])))
        a <- if (length(inter)) sum(as.numeric(st_area(inter))) / 1e6 else 0
        if (is.finite(a) && a > 1e-6) {
          k <- k + 1L
          overlap_rows[[k]] <- data.frame(
            ecoregion_1 = diss_laea$ecoregion_label[[i]],
            ecoregion_2 = diss_laea$ecoregion_label[[j]],
            overlap_km2 = a
          )
        }
      }
    }
  }
}
overlap_df <- if (length(overlap_rows)) do.call(rbind, overlap_rows) else data.frame(ecoregion_1 = character(), ecoregion_2 = character(), overlap_km2 = numeric())
write.csv(overlap_df, file.path(out_dir, "10A_pairwise_polygon_overlaps.csv"), row.names = FALSE)

sum_area <- sum(as.numeric(st_area(diss_laea))) / 1e6
union_geom <- st_union(st_geometry(diss_laea))
union_area <- as.numeric(st_area(union_geom)) / 1e6
overlap_total <- max(0, sum_area - union_area)

# Connected components identify the mainland versus islands without using the arachnid data.
components <- suppressWarnings(st_cast(union_geom, "POLYGON"))
component_sf <- st_sf(component_id = seq_along(components), geometry = components)
component_sf$area_km2 <- as.numeric(st_area(component_sf)) / 1e6
component_sf <- component_sf[order(component_sf$area_km2, decreasing = TRUE), ]
component_sf$area_rank <- seq_len(nrow(component_sf))
component_sf$is_largest_mainland_candidate <- component_sf$area_rank == 1L
cent <- suppressWarnings(st_point_on_surface(component_sf))
cent_wgs <- st_transform(cent, 4326)
xy <- st_coordinates(cent_wgs)
component_table <- st_drop_geometry(component_sf)
component_table$centroid_lon <- xy[, 1]
component_table$centroid_lat <- xy[, 2]
write.csv(component_table, file.path(out_dir, "10A_connected_component_audit.csv"), row.names = FALSE)

mainland_outline <- st_transform(component_sf[component_sf$area_rank == 1L, ], 4326)
st_write(mainland_outline, file.path(out_dir, "10A_mainland_outline_largest_component.gpkg"), delete_dsn = TRUE, quiet = TRUE)

# Adjacency/shared-boundary inventory for later boundary-specific tests.
adj_rows <- list()
k <- 0L
for (i in seq_len(nrow(diss_laea))) {
  for (j in seq_len(nrow(diss_laea))) {
    if (j <= i) next
    touches <- lengths(st_touches(diss_laea[i, ], diss_laea[j, ])) > 0
    if (touches) {
      shared <- suppressWarnings(st_intersection(st_boundary(st_geometry(diss_laea[i, ])), st_boundary(st_geometry(diss_laea[j, ]))))
      shared_km <- if (length(shared)) sum(as.numeric(st_length(shared))) / 1000 else 0
      k <- k + 1L
      adj_rows[[k]] <- data.frame(
        ecoregion_1 = diss_laea$ecoregion_label[[i]],
        ecoregion_2 = diss_laea$ecoregion_label[[j]],
        shared_boundary_km = shared_km
      )
    }
  }
}
adj_df <- if (length(adj_rows)) do.call(rbind, adj_rows) else data.frame(ecoregion_1 = character(), ecoregion_2 = character(), shared_boundary_km = numeric())
write.csv(adj_df, file.path(out_dir, "10A_ecoregion_adjacency.csv"), row.names = FALSE)

# Save validated layers. The mainland-only version is clipped to the largest connected land component.
st_write(eco_dissolved, file.path(out_dir, "10A_ecoregions_validated_all_components.gpkg"), delete_dsn = TRUE, quiet = TRUE)
mainland_wgs <- st_transform(mainland_outline, 4326)
eco_mainland <- suppressWarnings(st_intersection(eco_dissolved, st_geometry(mainland_wgs)))
eco_mainland <- eco_mainland[!st_is_empty(eco_mainland), ]
st_write(eco_mainland, file.path(out_dir, "10A_ecoregions_validated_mainland_only.gpkg"), delete_dsn = TRUE, quiet = TRUE)

# Locate candidate biological inputs but do not select one silently.
scan_root <- pipeline_root
candidate_inputs <- data.frame()
if (dir.exists(scan_root)) {
  all_files <- list.files(scan_root, recursive = TRUE, full.names = TRUE)
  keep_ext <- grepl("\\.(csv|tsv|rds|rda|gpkg|geojson|shp|parquet)$", all_files, ignore.case = TRUE)
  keep_name <- grepl("cell|25km|incidence|richness|turnover|jaccard|simpson|trait|C3|spatial", basename(all_files), ignore.case = TRUE)
  cand <- all_files[keep_ext & keep_name & !grepl("09_C3_biogeographic_concordance", all_files, fixed = TRUE)]
  if (length(cand)) {
    info <- file.info(cand)
    candidate_inputs <- data.frame(
      file = normalizePath(cand, winslash = "/", mustWork = FALSE),
      filename = basename(cand),
      extension = tolower(tools::file_ext(cand)),
      size_bytes = info$size,
      modified = format(info$mtime, "%Y-%m-%d %H:%M:%S")
    )
    candidate_inputs <- candidate_inputs[order(candidate_inputs$modified, decreasing = TRUE), ]
  }
}
write.csv(candidate_inputs, file.path(out_dir, "10A_candidate_project_inputs_for_step10B.csv"), row.names = FALSE)

# Publication-quality audit map (not yet a biological result figure).
p <- ggplot() +
  geom_sf(data = eco_mainland, aes(fill = ecoregion_label), color = "white", linewidth = 0.18) +
  coord_sf(datum = NA) +
  labs(
    title = "González-Abraham et al. (2010) Baja California ecoregions",
    subtitle = "Validated mainland component; audit map only",
    fill = "Ecoregion"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title.position = "plot"
  )

ggsave(file.path(pub_dir, "10A_validated_ecoregion_audit_map.png"), p, width = 10, height = 8, dpi = 400)
ggsave(file.path(pub_dir, "10A_validated_ecoregion_audit_map.pdf"), p, width = 10, height = 8, device = "pdf")

bbox <- st_bbox(eco_dissolved)
selected_summary <- data.frame(
  selected_source_file = selected_row$file,
  source_method = source_method,
  source_reference = source_reference,
  selected_layer_id = selected_row$layer_id,
  source_label_field = label_field,
  source_feature_count = nrow(eco),
  dissolved_ecoregion_count = nrow(eco_dissolved),
  invalid_before_repair = selected_row$invalid_geometries,
  invalid_after_repair = sum(!st_is_valid(eco), na.rm = TRUE),
  empty_after_repair = sum(st_is_empty(eco), na.rm = TRUE),
  exact_duplicate_geometry_features = sum(duplicate_geometry),
  pairwise_overlap_total_km2 = overlap_total,
  pairwise_overlap_percent_of_union = if (union_area > 0) 100 * overlap_total / union_area else NA_real_,
  connected_components = nrow(component_sf),
  mainland_candidate_area_km2 = component_sf$area_km2[[1]],
  other_components_area_km2 = sum(component_sf$area_km2[-1]),
  bbox_xmin_wgs84 = bbox[["xmin"]], bbox_ymin_wgs84 = bbox[["ymin"]],
  bbox_xmax_wgs84 = bbox[["xmax"]], bbox_ymax_wgs84 = bbox[["ymax"]],
  source_zip_md5 = zip_md5,
  seed = seed
)
write.csv(selected_summary, file.path(out_dir, "10A_selected_layer_summary.csv"), row.names = FALSE)

warnings <- character()
if (selected_summary$invalid_after_repair > 0) warnings <- c(warnings, "Some geometries remain invalid after st_make_valid().")
if (selected_summary$empty_after_repair > 0) warnings <- c(warnings, "Some geometries are empty after repair.")
if (selected_summary$exact_duplicate_geometry_features > 0) warnings <- c(warnings, "Exact duplicate geometries were detected.")
if (selected_summary$pairwise_overlap_percent_of_union > 0.01) warnings <- c(warnings, "Non-trivial overlap exists among dissolved ecoregion polygons.")
if (nrow(component_sf) > 1) warnings <- c(warnings, "The GIS contains disconnected components; the largest is retained as the mainland candidate for primary analysis.")
if (!nrow(candidate_inputs)) warnings <- c(warnings, "No candidate Step 10B biological input files were discovered automatically.")

audit_status <- if (any(grepl("remain invalid|Non-trivial overlap", warnings))) "REVIEW_REQUIRED" else "PASS_WITH_DOCUMENTED_COMPONENT_FILTER"

jsonlite::write_json(
  list(
    step = "10A",
    audit_status = audit_status,
    project_root = project_root,
    output_dir = out_dir,
    official_dataset_urls = urls,
    arcgis_fallback_item = arcgis_item_page,
    source_method = source_method,
    source_reference = source_reference,
    citation = "Gonzalez-Abraham CE, Garcillan PP, Ezcurra E, and Ecoregions Working Group. 2010. Ecorregiones de la peninsula de Baja California: una sintesis. Boletin de la Sociedad Botanica de Mexico 87:69-82.",
    selected_summary = as.list(selected_summary[1, ]),
    warnings = warnings,
    next_step = "Review the audit map, labels, disconnected components, and candidate biological inputs before running cell-to-ecoregion assignment and concordance tests."
  ),
  file.path(out_dir, "10A_audit_summary.json"),
  pretty = TRUE, auto_unbox = TRUE, na = "null"
)

readme_lines <- c(
  "STEP 10A — GONZALEZ-ABRAHAM ECOREGION GIS AUDIT",
  "",
  paste0("AUDIT_STATUS=", audit_status),
  paste0("OUTPUT_DIR=", out_dir),
  paste0("SELECTED_LAYER=", selected_row$file),
  paste0("LABEL_FIELD=", label_field),
  paste0("SOURCE_FEATURES=", nrow(eco)),
  paste0("DISSOLVED_ECOREGIONS=", nrow(eco_dissolved)),
  paste0("CONNECTED_COMPONENTS=", nrow(component_sf)),
  paste0("MAINLAND_CANDIDATE_AREA_KM2=", round(component_sf$area_km2[[1]], 2)),
  paste0("OTHER_COMPONENTS_AREA_KM2=", round(sum(component_sf$area_km2[-1]), 2)),
  paste0("OVERLAP_PERCENT=", signif(selected_summary$pairwise_overlap_percent_of_union, 5)),
  paste0("SOURCE_ZIP_MD5=", zip_md5),
  "",
  "WARNINGS:",
  if (length(warnings)) paste0("- ", warnings) else "- None",
  "",
  "DO NOT RUN THE BIOLOGICAL CONCORDANCE TESTS UNTIL THESE ARE REVIEWED:",
  "1. publication_outputs/10A_validated_ecoregion_audit_map.png",
  "2. 10A_ecoregion_labels_and_areas.csv",
  "3. 10A_connected_component_audit.csv",
  "4. 10A_pairwise_polygon_overlaps.csv",
  "5. 10A_candidate_project_inputs_for_step10B.csv",
  "",
  paste0("SOURCE_METHOD=", source_method),
  paste0("SOURCE_REFERENCE=", source_reference),
  "The script first attempts the official companion ZIP. If that retired URL returns 404, it uses the documented ArcGIS public-item fallback.",
  "Confirm redistribution terms before placing the downloaded shapefile itself in a public GitHub release."
)
writeLines(readme_lines, file.path(out_dir, "README_RESULTS_FIRST.txt"))

msg("")
msg("STEP 10A COMPLETE")
msg("AUDIT_STATUS=%s", audit_status)
msg("OUTPUT_DIR=%s", out_dir)
msg("SELECTED_LAYER=%s", selected_row$file)
msg("LABEL_FIELD=%s", label_field)
msg("DISSOLVED_ECOREGIONS=%d", nrow(eco_dissolved))
msg("CONNECTED_COMPONENTS=%d", nrow(component_sf))
msg("PAIRWISE_OVERLAP_PERCENT=%.6f", selected_summary$pairwise_overlap_percent_of_union)
msg("SOURCE_ZIP_MD5=%s", zip_md5)
