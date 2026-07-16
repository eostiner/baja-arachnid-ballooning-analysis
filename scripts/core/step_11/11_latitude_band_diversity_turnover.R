#!/usr/bin/env Rscript

# ============================================================================
# STEP 11 v3 — LATITUDE-BAND DIVERSITY, COMPLETENESS, TURNOVER, AND C3/N0 TRAITS
# Baja Ballooning Publication
#
# Default project:
#   ~/Desktop/Baja_Ballooning_Pipeline
#
# Run:
#   Rscript 11_latitude_band_diversity_turnover.R
#
# Optional:
#   Rscript 11_latitude_band_diversity_turnover.R \
#     ~/Desktop/Baja_Ballooning_Pipeline 200
#
# The second argument is the number of iNEXT bootstrap replicates.
# ============================================================================

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  warn = 1
)

# ------------------------------ arguments ----------------------------------

args <- commandArgs(
  trailingOnly = TRUE
)

project_root <- if (
  length(args) >= 1 &&
  nzchar(args[[1]])
) {
  path.expand(
    args[[1]]
  )
} else {
  path.expand(
    "~/Desktop/Baja_Ballooning_Pipeline"
  )
}

nboot <- if (
  length(args) >= 2 &&
  nzchar(args[[2]])
) {
  as.integer(
    args[[2]]
  )
} else {
  200L
}

if (
  is.na(nboot) ||
  nboot < 0
) {
  stop(
    "The iNEXT bootstrap count must be a non-negative integer."
  )
}

set.seed(
  20260713
)

band_order <- c(
  "23-24N",
  "24-26N",
  "26-28N",
  "28-30N",
  "30-32N"
)

band_labels <- c(
  "23–24°N",
  "24–26°N",
  "26–28°N",
  "28–30°N",
  "30–32°N"
)

names(
  band_labels
) <- band_order

# ------------------------------ packages -----------------------------------

required_packages <- c(
  "iNEXT",
  "ggplot2",
  "jsonlite"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(
      1
    )
  )
]

if (
  length(
    missing_packages
  ) > 0
) {
  message(
    "Installing missing R packages: ",
    paste(
      missing_packages,
      collapse = ", "
    )
  )

  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org",
    dependencies = TRUE
  )
}

suppressPackageStartupMessages(
  library(
    iNEXT
  )
)

suppressPackageStartupMessages(
  library(
    ggplot2
  )
)

suppressPackageStartupMessages(
  library(
    jsonlite
  )
)

# ------------------------------ folders ------------------------------------

if (
  !dir.exists(
    project_root
  )
) {
  stop(
    "Project folder not found:\n",
    project_root
  )
}

analysis_ready <- file.path(
  project_root,
  "ANALYSIS_READY_INPUTS"
)

grid_fallback <- file.path(
  project_root,
  "02_data_clean",
  "08_grid25km_incidence"
)

trait_fallback <- file.path(
  project_root,
  "02_data_clean",
  "07_final_trait_merge"
)

output_dir <- file.path(
  project_root,
  "04_analysis",
  "11_latitude_band_diversity_turnover"
)

archive_root <- file.path(
  project_root,
  "08_archive"
)

timestamp <- format(
  Sys.time(),
  "%Y%m%dT%H%M%S"
)

if (
  dir.exists(
    output_dir
  )
) {
  existing_files <- list.files(
    output_dir,
    all.files = TRUE,
    no.. = TRUE
  )

  if (
    length(
      existing_files
    ) > 0
  ) {
    archive_dir <- file.path(
      archive_root,
      paste0(
        "11_latitude_band_diversity_turnover_",
        timestamp
      )
    )

    dir.create(
      archive_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    copied <- file.copy(
      from = file.path(
        output_dir,
        existing_files
      ),
      to = archive_dir,
      recursive = TRUE,
      overwrite = TRUE
    )

    if (
      !all(
        copied
      )
    ) {
      stop(
        "Could not archive all existing Step 11 outputs."
      )
    }

    unlink(
      output_dir,
      recursive = TRUE,
      force = TRUE
    )
  }
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

figure_dir <- file.path(
  output_dir,
  "figures"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

log_path <- file.path(
  output_dir,
  "11_analysis_log.txt"
)

log_connection <- file(
  log_path,
  open = "wt"
)

sink(
  log_connection,
  split = TRUE
)

sink(
  log_connection,
  type = "message",
  append = TRUE
)

on.exit(
  {
    while (
      sink.number(
        type = "message"
      ) > 0
    ) {
      sink(
        type = "message"
      )
    }

    while (
      sink.number() > 0
    ) {
      sink()
    }

    close(
      log_connection
    )
  },
  add = TRUE
)

message(
  "STEP 11 STARTED"
)

message(
  "Project: ",
  project_root
)

message(
  "iNEXT bootstrap replicates: ",
  nboot
)

# ------------------------------ helpers ------------------------------------

first_existing <- function(
  paths,
  required = TRUE,
  label = "file"
) {
  paths <- path.expand(
    paths
  )

  existing <- paths[
    file.exists(
      paths
    )
  ]

  if (
    length(
      existing
    ) > 0
  ) {
    return(
      existing[[1]]
    )
  }

  if (
    required
  ) {
    stop(
      "Could not find ",
      label,
      ". Tried:\n",
      paste(
        paths,
        collapse = "\n"
      )
    )
  }

  return(
    NA_character_
  )
}


find_field <- function(
  fields,
  candidates,
  required = TRUE,
  label = "field"
) {
  for (
    candidate in candidates
  ) {
    if (
      candidate %in% fields
    ) {
      return(
        candidate
      )
    }
  }

  if (
    required
  ) {
    stop(
      "Could not identify ",
      label,
      ". Tried: ",
      paste(
        candidates,
        collapse = ", "
      )
    )
  }

  return(
    NA_character_
  )
}


trim_character <- function(
  value
) {
  trimws(
    as.character(
      value
    )
  )
}


read_incidence_matrix <- function(
  path
) {
  data <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (
    ncol(
      data
    ) < 2
  ) {
    stop(
      "Incidence matrix has fewer than two columns:\n",
      path
    )
  }

  genus_field <- find_field(
    names(
      data
    ),
    c(
      "genus",
      "analysis_genus"
    ),
    label = "genus column"
  )

  genera <- trim_character(
    data[[genus_field]]
  )

  if (
    any(
      !nzchar(
        genera
      )
    )
  ) {
    stop(
      "Blank genus names found in:\n",
      path
    )
  }

  if (
    anyDuplicated(
      tolower(
        genera
      )
    )
  ) {
    stop(
      "Duplicate genus names found in:\n",
      path
    )
  }

  matrix_columns <- setdiff(
    names(
      data
    ),
    genus_field
  )

  matrix_data <- as.matrix(
    data[
      matrix_columns
    ]
  )

  suppressWarnings(
    storage.mode(
      matrix_data
    ) <- "numeric"
  )

  if (
    any(
      is.na(
        matrix_data
      )
    )
  ) {
    stop(
      "Non-numeric or missing matrix entries found in:\n",
      path
    )
  }

  if (
    any(
      !matrix_data %in% c(
        0,
        1
      )
    )
  ) {
    stop(
      "Matrix contains entries other than 0 and 1:\n",
      path
    )
  }

  rownames(
    matrix_data
  ) <- genera

  list(
    path = path,
    matrix = matrix_data,
    genera = genera,
    cells = matrix_columns
  )
}


write_matrix_csv <- function(
  matrix_object,
  path,
  row_name = "genus"
) {
  output <- data.frame(
    row_name_value = rownames(
      matrix_object
    ),
    matrix_object,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  names(
    output
  )[
    1
  ] <- row_name

  write.csv(
    output,
    path,
    row.names = FALSE,
    na = ""
  )
}


safe_numeric <- function(
  value
) {
  suppressWarnings(
    as.numeric(
      trim_character(
        value
      )
    )
  )
}


normalize_status_to_binary <- function(
  status
) {
  status <- tolower(
    trim_character(
      status
    )
  )

  output <- rep(
    NA_integer_,
    length(
      status
    )
  )

  output[
    status %in% c(
      "ballooning",
      "ballooner",
      "yes",
      "true",
      "1"
    )
  ] <- 1L

  output[
    status %in% c(
      "non_ballooning",
      "non-ballooning",
      "nonballooning",
      "no",
      "false",
      "0"
    )
  ] <- 0L

  output
}


load_traits <- function(
  path,
  required_genera
) {
  traits <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  genus_field <- find_field(
    names(traits),
    c("genus", "analysis_genus"),
    label = "trait genus column"
  )

  confidence_field <- find_field(
    names(traits),
    c(
      "final_confidence",
      "trait_final_confidence",
      "trait_confidence",
      "trait_ballooning_confidence"
    ),
    required = FALSE
  )
  order_field <- find_field(names(traits), c("order", "taxon_order"), required = FALSE)
  family_field <- find_field(names(traits), c("family", "taxon_family"), required = FALSE)

  parse_evidence_class <- function(value) {
    text <- toupper(trimws(as.character(value)))
    if (is.na(text) || !nzchar(text)) return(NA_character_)
    hits <- regmatches(
      text,
      gregexpr(
        "(?<![A-Z0-9])(D1|D2|D3|D4|N0|C3)(?![A-Z0-9])",
        text,
        perl = TRUE
      )
    )[[1]]
    hits <- unique(hits[hits != ""])
    if (length(hits) == 1L && hits %in% c("D1", "D2", "D3", "D4", "N0", "C3")) {
      return(hits)
    }
    normalized <- tolower(gsub("[^a-z0-9]+", "", text))
    if (normalized %in% c(
      "nonballooning", "fixednonballooning", "referencenonballooning",
      "noballooning", "nonballooningreference"
    )) return("N0")
    if (normalized %in% c("c3", "primaryc3", "d1d2d3", "d1tod3")) return("C3")
    if (normalized %in% c("d4excluded", "excludedd4")) return("D4")
    NA_character_
  }

  preferred_fields <- c(
    "evidence_class", "final_evidence_class", "final_evidence_category",
    "evidence_category", "evidence_level", "trait_evidence_level",
    "d_level", "dlevel", "trait_class", "primary_class", "analysis_class",
    "ballooning_evidence_tier", "ballooning_evidence_category",
    "final_designation", "designation"
  )
  candidate_scores <- vapply(
    names(traits),
    function(field) {
      parsed <- vapply(traits[[field]], parse_evidence_class, character(1))
      n_ok <- sum(!is.na(parsed))
      if (n_ok == 0L) return(-Inf)
      fraction <- n_ok / max(1L, nrow(traits))
      classes <- length(unique(parsed[!is.na(parsed)]))
      field_clean <- tolower(gsub("[^a-z0-9]+", "", field))
      bonus <- 0
      if (tolower(field) %in% tolower(preferred_fields)) bonus <- bonus + 100
      if (grepl("evidence|tier|class|designation|decision", field_clean)) bonus <- bonus + 20
      if (classes < 2L || fraction < 0.25) return(-Inf)
      bonus + 100 * fraction + 5 * classes
    },
    numeric(1)
  )
  if (!any(is.finite(candidate_scores))) {
    stop(
      "The trait table lacks an explicit D1/D2/D3/D4/N0 or C3/N0 evidence field. ",
      "Legacy binary ballooning fields are intentionally rejected because they cannot ",
      "distinguish D4 from fixed N0."
    )
  }
  evidence_field <- names(which.max(candidate_scores))

  genus <- trim_character(traits[[genus_field]])
  if (anyDuplicated(tolower(genus))) {
    duplicated_genera <- unique(genus[
      duplicated(tolower(genus)) | duplicated(tolower(genus), fromLast = TRUE)
    ])
    stop("Duplicate genera in trait table:\n", paste(duplicated_genera, collapse = "\n"))
  }

  evidence <- vapply(traits[[evidence_field]], parse_evidence_class, character(1))
  analysis_class <- ifelse(
    evidence %in% c("D1", "D2", "D3", "C3"),
    "C3",
    ifelse(evidence == "N0", "N0", ifelse(evidence == "D4", "D4_excluded", NA_character_))
  )
  confidence <- if (is.na(confidence_field)) {
    rep("UNSPECIFIED", nrow(traits))
  } else {
    toupper(trim_character(traits[[confidence_field]]))
  }
  confidence[!nzchar(confidence)] <- "UNSPECIFIED"
  order_value <- if (is.na(order_field)) rep("", nrow(traits)) else trim_character(traits[[order_field]])
  family_value <- if (is.na(family_field)) rep("", nrow(traits)) else trim_character(traits[[family_field]])

  # Compatibility field: 1 = C3, 0 = fixed N0, -1 = D4 excluded.
  binary <- ifelse(analysis_class == "C3", 1L, ifelse(analysis_class == "N0", 0L, -1L))
  normalized <- data.frame(
    genus = genus,
    order = order_value,
    family = family_value,
    ballooning_binary = binary,
    ballooning_status = ifelse(
      analysis_class == "C3", "ballooning_C3",
      ifelse(analysis_class == "N0", "non_ballooning_N0", "excluded_D4")
    ),
    analysis_class = analysis_class,
    final_confidence = confidence,
    evidence_level = evidence,
    evidence_field = evidence_field,
    stringsAsFactors = FALSE
  )

  index <- match(tolower(required_genera), tolower(normalized$genus))
  if (any(is.na(index))) {
    stop(
      "The following matrix genera are missing from the final trait table:\n",
      paste(required_genera[is.na(index)], collapse = "\n")
    )
  }
  normalized <- normalized[index, , drop = FALSE]
  if (any(is.na(normalized$analysis_class))) {
    stop(
      "The following incidence-matrix genera lack an explicit D1/D2/D3/D4/N0 or C3/N0 designation:\n",
      paste(normalized$genus[is.na(normalized$analysis_class)], collapse = "\n")
    )
  }
  if (any(!normalized$ballooning_binary %in% c(-1L, 0L, 1L))) {
    stop("Invalid C3/N0/D4 values remain in the normalized trait table.")
  }
  normalized
}

make_band_matrices <- function(
  incidence_matrix,
  cell_lookup
) {
  output <- vector(
    "list",
    length(
      band_order
    )
  )

  names(
    output
  ) <- band_order

  for (
    band in band_order
  ) {
    cells <- cell_lookup$grid_cell_id[
      cell_lookup$latitude_band == band
    ]

    output[[band]] <- incidence_matrix[
      ,
      cells,
      drop = FALSE
    ]
  }

  output
}


make_genus_by_band <- function(
  band_matrices
) {
  output <- sapply(
    band_order,
    function(
      band
    ) {
      as.integer(
        rowSums(
          band_matrices[[band]]
        ) > 0
      )
    }
  )

  if (
    is.null(
      dim(
        output
      )
    )
  ) {
    output <- matrix(
      output,
      ncol = length(
        band_order
      )
    )
  }

  rownames(
    output
  ) <- rownames(
    band_matrices[[1]]
  )

  colnames(
    output
  ) <- band_order

  output
}


prepare_inext_list <- function(
  band_matrices
) {
  output <- lapply(
    band_order,
    function(
      band
    ) {
      matrix_object <- band_matrices[[band]]

      observed <- rowSums(
        matrix_object
      ) > 0

      matrix_object[
        observed,
        ,
        drop = FALSE
      ]
    }
  )

  names(
    output
  ) <- band_order

  output
}


clean_inext_table <- function(
  object
) {
  output <- as.data.frame(
    object,
    stringsAsFactors = FALSE
  )

  if (
    "Assemblage" %in% names(
      output
    )
  ) {
    output$Assemblage <- factor(
      output$Assemblage,
      levels = band_order
    )

    output <- output[
      order(
        output$Assemblage
      ),
      ,
      drop = FALSE
    ]

    output$Assemblage <- as.character(
      output$Assemblage
    )
  }

  rownames(
    output
  ) <- NULL

  output
}


extract_inext_curves <- function(
  inext_object
) {
  size_based <- NULL
  coverage_based <- NULL

  if (
    is.list(
      inext_object$iNextEst
    ) &&
    "size_based" %in% names(
      inext_object$iNextEst
    )
  ) {
    size_based <- inext_object$iNextEst$size_based
    coverage_based <- inext_object$iNextEst$coverage_based
  } else {
    size_based <- inext_object$iNextEst
  }

  list(
    size_based = clean_inext_table(
      size_based
    ),
    coverage_based = if (
      is.null(
        coverage_based
      )
    ) {
      data.frame()
    } else {
      clean_inext_table(
        coverage_based
      )
    }
  )
}


run_inext_analysis <- function(
  dataset_key,
  band_matrices,
  output_dir,
  nboot
) {
  assemblages <- prepare_inext_list(
    band_matrices
  )

  sample_sizes <- vapply(
    assemblages,
    ncol,
    FUN.VALUE = integer(
      1
    )
  )

  minimum_cells <- min(
    sample_sizes
  )

  data_info <- iNEXT::DataInfo(
    assemblages,
    datatype = "incidence_raw"
  )

  rarefaction <- iNEXT::iNEXT(
    assemblages,
    q = 0,
    datatype = "incidence_raw",
    endpoint = NULL,
    knots = 40,
    se = TRUE,
    conf = 0.95,
    nboot = nboot
  )

  curves <- extract_inext_curves(
    rarefaction
  )

  equal_cell <- iNEXT::estimateD(
    assemblages,
    q = 0,
    datatype = "incidence_raw",
    base = "size",
    level = minimum_cells,
    nboot = nboot,
    conf = 0.95
  )

  coverage_standardized <- iNEXT::estimateD(
    assemblages,
    q = 0,
    datatype = "incidence_raw",
    base = "coverage",
    level = NULL,
    nboot = nboot,
    conf = 0.95
  )

  asymptotic <- rarefaction$AsyEst

  write.csv(
    clean_inext_table(
      data_info
    ),
    file.path(
      output_dir,
      paste0(
        "11_",
        dataset_key,
        "_inext_data_info.csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )

  write.csv(
    curves$size_based,
    file.path(
      output_dir,
      paste0(
        "11_",
        dataset_key,
        "_inext_size_based_curve.csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )

  if (
    nrow(
      curves$coverage_based
    ) > 0
  ) {
    write.csv(
      curves$coverage_based,
      file.path(
        output_dir,
        paste0(
          "11_",
          dataset_key,
          "_inext_coverage_based_curve.csv"
        )
      ),
      row.names = FALSE,
      na = ""
    )
  }

  write.csv(
    clean_inext_table(
      equal_cell
    ),
    file.path(
      output_dir,
      paste0(
        "11_",
        dataset_key,
        "_equal_cell_standardized_richness.csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )

  write.csv(
    clean_inext_table(
      coverage_standardized
    ),
    file.path(
      output_dir,
      paste0(
        "11_",
        dataset_key,
        "_coverage_standardized_richness.csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )

  write.csv(
    clean_inext_table(
      asymptotic
    ),
    file.path(
      output_dir,
      paste0(
        "11_",
        dataset_key,
        "_asymptotic_richness.csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )

  size_plot <- iNEXT::ggiNEXT(
    rarefaction,
    type = 1,
    se = TRUE,
    facet.var = "None",
    color.var = "Assemblage"
  ) +
    ggplot2::theme_bw(
      base_size = 12
    ) +
    ggplot2::labs(
      title = paste0(
        dataset_key,
        ": incidence-based richness by sampling effort"
      ),
      x = "Number of occupied 25-km grid cells",
      y = "Estimated genus richness",
      color = "Latitude band",
      fill = "Latitude band"
    )

  completeness_plot <- iNEXT::ggiNEXT(
    rarefaction,
    type = 2,
    se = TRUE,
    facet.var = "None",
    color.var = "Assemblage"
  ) +
    ggplot2::theme_bw(
      base_size = 12
    ) +
    ggplot2::labs(
      title = paste0(
        dataset_key,
        ": sample completeness"
      ),
      x = "Number of occupied 25-km grid cells",
      y = "Sample coverage",
      color = "Latitude band",
      fill = "Latitude band"
    )

  coverage_plot <- iNEXT::ggiNEXT(
    rarefaction,
    type = 3,
    se = TRUE,
    facet.var = "None",
    color.var = "Assemblage"
  ) +
    ggplot2::theme_bw(
      base_size = 12
    ) +
    ggplot2::labs(
      title = paste0(
        dataset_key,
        ": coverage-standardized richness"
      ),
      x = "Sample coverage",
      y = "Estimated genus richness",
      color = "Latitude band",
      fill = "Latitude band"
    )

  ggplot2::ggsave(
    file.path(
      figure_dir,
      paste0(
        "11_",
        dataset_key,
        "_rarefaction_by_cells.png"
      )
    ),
    size_plot,
    width = 8.5,
    height = 6.2,
    dpi = 320
  )

  ggplot2::ggsave(
    file.path(
      figure_dir,
      paste0(
        "11_",
        dataset_key,
        "_sample_completeness.png"
      )
    ),
    completeness_plot,
    width = 8.5,
    height = 6.2,
    dpi = 320
  )

  ggplot2::ggsave(
    file.path(
      figure_dir,
      paste0(
        "11_",
        dataset_key,
        "_coverage_standardized_curve.png"
      )
    ),
    coverage_plot,
    width = 8.5,
    height = 6.2,
    dpi = 320
  )

  list(
    assemblages = assemblages,
    sample_sizes = sample_sizes,
    minimum_cells = minimum_cells,
    data_info = clean_inext_table(
      data_info
    ),
    equal_cell = clean_inext_table(
      equal_cell
    ),
    coverage_standardized = clean_inext_table(
      coverage_standardized
    ),
    asymptotic = clean_inext_table(
      asymptotic
    ),
    rarefaction = rarefaction
  )
}


pairwise_beta <- function(
  genus_by_band
) {
  band_by_genus <- t(
    genus_by_band
  )

  output <- list()
  index <- 1L

  for (
    first_index in seq_len(
      length(
        band_order
      ) - 1L
    )
  ) {
    for (
      second_index in seq.int(
        first_index + 1L,
        length(
          band_order
        )
      )
    ) {
      first_band <- band_order[
        first_index
      ]

      second_band <- band_order[
        second_index
      ]

      first_presence <- band_by_genus[
        first_band,
      ] > 0

      second_presence <- band_by_genus[
        second_band,
      ] > 0

      shared <- sum(
        first_presence &
          second_presence
      )

      unique_first <- sum(
        first_presence &
          !second_presence
      )

      unique_second <- sum(
        !first_presence &
          second_presence
      )

      jaccard_denominator <- shared +
        unique_first +
        unique_second

      sorensen_denominator <- 2 * shared +
        unique_first +
        unique_second

      minimum_unique <- min(
        unique_first,
        unique_second
      )

      simpson_denominator <- shared +
        minimum_unique

      jaccard <- if (
        jaccard_denominator == 0
      ) {
        NA_real_
      } else {
        (
          unique_first +
            unique_second
        ) /
          jaccard_denominator
      }

      sorensen <- if (
        sorensen_denominator == 0
      ) {
        NA_real_
      } else {
        (
          unique_first +
            unique_second
        ) /
          sorensen_denominator
      }

      simpson_turnover <- if (
        simpson_denominator == 0
      ) {
        NA_real_
      } else {
        minimum_unique /
          simpson_denominator
      }

      nestedness_resultant <- if (
        is.na(
          sorensen
        ) ||
        is.na(
          simpson_turnover
        )
      ) {
        NA_real_
      } else {
        max(
          0,
          sorensen -
            simpson_turnover
        )
      }

      output[[index]] <- data.frame(
        band_1 = first_band,
        band_2 = second_band,
        shared_genera = shared,
        unique_to_band_1 = unique_first,
        unique_to_band_2 = unique_second,
        jaccard_dissimilarity = jaccard,
        sorensen_dissimilarity = sorensen,
        simpson_turnover = simpson_turnover,
        sorensen_nestedness_resultant = nestedness_resultant,
        adjacent_bands = second_index == first_index + 1L,
        stringsAsFactors = FALSE
      )

      index <- index + 1L
    }
  }

  do.call(
    rbind,
    output
  )
}


pairwise_to_matrix <- function(
  pairwise_table,
  value_field
) {
  output <- matrix(
    0,
    nrow = length(
      band_order
    ),
    ncol = length(
      band_order
    ),
    dimnames = list(
      band_order,
      band_order
    )
  )

  for (
    row_index in seq_len(
      nrow(
        pairwise_table
      )
    )
  ) {
    first_band <- pairwise_table$band_1[
      row_index
    ]

    second_band <- pairwise_table$band_2[
      row_index
    ]

    value <- pairwise_table[[value_field]][
      row_index
    ]

    output[
      first_band,
      second_band
    ] <- value

    output[
      second_band,
      first_band
    ] <- value
  }

  output
}


write_beta_outputs <- function(
  dataset_key,
  genus_by_band,
  output_dir
) {
  pairwise <- pairwise_beta(
    genus_by_band
  )

  write.csv(
    pairwise,
    file.path(
      output_dir,
      paste0(
        "11_",
        dataset_key,
        "_pairwise_beta_diversity.csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )

  metrics <- c(
    "jaccard_dissimilarity",
    "sorensen_dissimilarity",
    "simpson_turnover",
    "sorensen_nestedness_resultant"
  )

  matrix_outputs <- list()

  for (
    metric in metrics
  ) {
    metric_matrix <- pairwise_to_matrix(
      pairwise,
      metric
    )

    matrix_outputs[[metric]] <- metric_matrix

    write.csv(
      data.frame(
        latitude_band = rownames(
          metric_matrix
        ),
        metric_matrix,
        check.names = FALSE
      ),
      file.path(
        output_dir,
        paste0(
          "11_",
          dataset_key,
          "_",
          metric,
          "_matrix.csv"
        )
      ),
      row.names = FALSE,
      na = ""
    )
  }

  adjacent <- pairwise[
    pairwise$adjacent_bands,
    ,
    drop = FALSE
  ]

  write.csv(
    adjacent,
    file.path(
      output_dir,
      paste0(
        "11_",
        dataset_key,
        "_adjacent_band_turnover.csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )

  heatmap_data <- pairwise

  heatmap_data$band_1 <- factor(
    heatmap_data$band_1,
    levels = rev(
      band_order
    )
  )

  heatmap_data$band_2 <- factor(
    heatmap_data$band_2,
    levels = band_order
  )

  heatmap_plot <- ggplot2::ggplot(
    heatmap_data,
    ggplot2::aes(
      x = band_2,
      y = band_1,
      fill = jaccard_dissimilarity
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(
      ggplot2::aes(
        label = sprintf(
          "%.2f",
          jaccard_dissimilarity
        )
      ),
      size = 3.5
    ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "black",
      limits = c(
        0,
        1
      )
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_bw(
      base_size = 12
    ) +
    ggplot2::labs(
      title = paste0(
        dataset_key,
        ": pairwise Jaccard dissimilarity"
      ),
      x = "Latitude band",
      y = "Latitude band",
      fill = "Jaccard"
    )

  adjacent_plot <- ggplot2::ggplot(
    adjacent,
    ggplot2::aes(
      x = factor(
        paste(
          band_1,
          band_2,
          sep = " to "
        ),
        levels = paste(
          adjacent$band_1,
          adjacent$band_2,
          sep = " to "
        )
      ),
      y = simpson_turnover
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_point(
      ggplot2::aes(
        y = sorensen_nestedness_resultant
      ),
      size = 2.5
    ) +
    ggplot2::theme_bw(
      base_size = 12
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 35,
        hjust = 1
      )
    ) +
    ggplot2::labs(
      title = paste0(
        dataset_key,
        ": adjacent-band beta-diversity components"
      ),
      x = "Adjacent latitude bands",
      y = "Dissimilarity",
      caption = "Bars: Simpson turnover; points: Sørensen nestedness-resultant component"
    )

  ggplot2::ggsave(
    file.path(
      figure_dir,
      paste0(
        "11_",
        dataset_key,
        "_jaccard_heatmap.png"
      )
    ),
    heatmap_plot,
    width = 7.5,
    height = 6.5,
    dpi = 320
  )

  ggplot2::ggsave(
    file.path(
      figure_dir,
      paste0(
        "11_",
        dataset_key,
        "_adjacent_turnover_components.png"
      )
    ),
    adjacent_plot,
    width = 8.5,
    height = 6,
    dpi = 320
  )

  list(
    pairwise = pairwise,
    adjacent = adjacent,
    matrices = matrix_outputs
  )
}


band_sampling_summary <- function(
  band_matrices,
  cell_lookup,
  data_info = NULL,
  dataset_key = "biodiversity_final"
) {
  output <- lapply(
    band_order,
    function(
      band
    ) {
      matrix_object <- band_matrices[[band]]

      incidence_frequency <- rowSums(
        matrix_object
      )

      cell_richness <- colSums(
        matrix_object
      )

      band_cells <- cell_lookup[
        cell_lookup$latitude_band == band,
        ,
        drop = FALSE
      ]

      record_count <- if (
        "biodiversity_record_count" %in% names(
          band_cells
        ) &&
        dataset_key == "biodiversity_final"
      ) {
        sum(
          safe_numeric(
            band_cells$biodiversity_record_count
          ),
          na.rm = TRUE
        )
      } else {
        NA_real_
      }

      data.frame(
        dataset = dataset_key,
        latitude_band = band,
        latitude_band_label = unname(
          band_labels[
            band
          ]
        ),
        occupied_25km_cells = ncol(
          matrix_object
        ),
        occurrence_records = record_count,
        genus_cell_incidences = sum(
          matrix_object
        ),
        observed_genus_richness = sum(
          incidence_frequency > 0
        ),
        genera_in_one_cell = sum(
          incidence_frequency == 1
        ),
        genera_in_two_cells = sum(
          incidence_frequency == 2
        ),
        mean_genus_richness_per_cell = mean(
          cell_richness
        ),
        median_genus_richness_per_cell = stats::median(
          cell_richness
        ),
        minimum_genus_richness_per_cell = min(
          cell_richness
        ),
        maximum_genus_richness_per_cell = max(
          cell_richness
        ),
        stringsAsFactors = FALSE
      )
    }
  )

  output <- do.call(
    rbind,
    output
  )

  if (
    !is.null(
      data_info
    ) &&
    "Assemblage" %in% names(
      data_info
    )
  ) {
    coverage_field <- find_field(
      names(
        data_info
      ),
      c(
        "SC",
        "SampleCoverage",
        "sample_coverage"
      ),
      required = FALSE
    )

    if (
      !is.na(
        coverage_field
      )
    ) {
      coverage_lookup <- data.frame(
        latitude_band = as.character(
          data_info$Assemblage
        ),
        estimated_sample_coverage = safe_numeric(
          data_info[[coverage_field]]
        ),
        stringsAsFactors = FALSE
      )

      output <- merge(
        output,
        coverage_lookup,
        by = "latitude_band",
        all.x = TRUE,
        sort = FALSE
      )

      output <- output[
        match(
          band_order,
          output$latitude_band
        ),
        ,
        drop = FALSE
      ]
    }
  }

  rownames(
    output
  ) <- NULL

  output
}


cell_richness_table <- function(
  incidence_matrix,
  cell_lookup,
  traits
) {
  ballooning_index <- traits$ballooning_binary == 1L
  non_ballooning_index <- traits$ballooning_binary == 0L
  d4_index <- traits$ballooning_binary == -1L
  classified_index <- ballooning_index | non_ballooning_index
  confidence_restricted_index <- traits$final_confidence != "LOW" & classified_index

  ballooning_matrix <- incidence_matrix[ballooning_index, , drop = FALSE]
  non_ballooning_matrix <- incidence_matrix[non_ballooning_index, , drop = FALSE]
  d4_matrix <- incidence_matrix[d4_index, , drop = FALSE]
  classified_matrix <- incidence_matrix[classified_index, , drop = FALSE]
  restricted_matrix <- incidence_matrix[confidence_restricted_index, , drop = FALSE]
  restricted_traits <- traits[confidence_restricted_index, , drop = FALSE]
  restricted_ballooning <- restricted_traits$ballooning_binary == 1L
  restricted_non_ballooning <- restricted_traits$ballooning_binary == 0L

  output <- cell_lookup
  cells <- output$grid_cell_id
  output$total_genus_richness <- colSums(incidence_matrix[, cells, drop = FALSE])
  output$classified_C3_N0_genus_richness <- colSums(classified_matrix[, cells, drop = FALSE])
  output$ballooning_genus_richness <- colSums(ballooning_matrix[, cells, drop = FALSE])
  output$non_ballooning_genus_richness <- colSums(non_ballooning_matrix[, cells, drop = FALSE])
  output$excluded_D4_genus_richness <- colSums(d4_matrix[, cells, drop = FALSE])
  output$ballooning_genus_proportion <- ifelse(
    output$classified_C3_N0_genus_richness > 0,
    output$ballooning_genus_richness / output$classified_C3_N0_genus_richness,
    NA_real_
  )

  output$confidence_restricted_total_richness <- colSums(restricted_matrix[, cells, drop = FALSE])
  output$confidence_restricted_ballooning_richness <- colSums(
    restricted_matrix[restricted_ballooning, cells, drop = FALSE]
  )
  output$confidence_restricted_non_ballooning_richness <- colSums(
    restricted_matrix[restricted_non_ballooning, cells, drop = FALSE]
  )
  output$confidence_restricted_ballooning_proportion <- ifelse(
    output$confidence_restricted_total_richness > 0,
    output$confidence_restricted_ballooning_richness /
      output$confidence_restricted_total_richness,
    NA_real_
  )
  output
}


ballooning_band_summary <- function(
  genus_by_band,
  band_matrices,
  traits
) {
  output <- lapply(
    band_order,
    function(band) {
      present <- genus_by_band[, band] > 0
      ballooning <- present & traits$ballooning_binary == 1L
      non_ballooning <- present & traits$ballooning_binary == 0L
      excluded_d4 <- present & traits$ballooning_binary == -1L
      classified_present <- ballooning | non_ballooning

      retained_confidence <- traits$final_confidence != "LOW"
      restricted_ballooning <- ballooning & retained_confidence
      restricted_non_ballooning <- non_ballooning & retained_confidence
      restricted_present <- restricted_ballooning | restricted_non_ballooning

      matrix_object <- band_matrices[[band]]
      ballooning_incidence <- sum(matrix_object[traits$ballooning_binary == 1L, , drop = FALSE])
      non_ballooning_incidence <- sum(matrix_object[traits$ballooning_binary == 0L, , drop = FALSE])
      d4_incidence <- sum(matrix_object[traits$ballooning_binary == -1L, , drop = FALSE])
      classified_incidence <- ballooning_incidence + non_ballooning_incidence

      data.frame(
        latitude_band = band,
        latitude_band_label = unname(band_labels[band]),
        total_observed_genera = sum(present),
        classified_C3_N0_genera = sum(classified_present),
        ballooning_genera = sum(ballooning),
        non_ballooning_genera = sum(non_ballooning),
        excluded_D4_genera = sum(excluded_d4),
        ballooning_genus_proportion = if (sum(classified_present) == 0) {
          NA_real_
        } else {
          sum(ballooning) / sum(classified_present)
        },
        genus_cell_incidences = sum(matrix_object),
        classified_C3_N0_genus_cell_incidences = classified_incidence,
        ballooning_genus_cell_incidences = ballooning_incidence,
        non_ballooning_genus_cell_incidences = non_ballooning_incidence,
        excluded_D4_genus_cell_incidences = d4_incidence,
        ballooning_incidence_proportion = if (classified_incidence == 0) {
          NA_real_
        } else {
          ballooning_incidence / classified_incidence
        },
        low_confidence_genera_present = sum(present & traits$final_confidence == "LOW"),
        confidence_restricted_total_genera = sum(restricted_present),
        confidence_restricted_ballooning_genera = sum(restricted_ballooning),
        confidence_restricted_non_ballooning_genera = sum(restricted_non_ballooning),
        confidence_restricted_ballooning_proportion = if (sum(restricted_present) == 0) {
          NA_real_
        } else {
          sum(restricted_ballooning) / sum(restricted_present)
        },
        stringsAsFactors = FALSE
      )
    }
  )
  do.call(rbind, output)
}

standardized_field <- function(
  table,
  candidates,
  label
) {
  field <- find_field(
    names(
      table
    ),
    candidates,
    required = FALSE
  )

  if (
    is.na(
      field
    )
  ) {
    stop(
      "Could not identify ",
      label,
      " in standardized-richness output."
    )
  }

  field
}


plot_standardized_richness <- function(
  table,
  title,
  filename
) {
  assemblage_field <- standardized_field(
    table,
    c(
      "Assemblage",
      "assemblage"
    ),
    "assemblage field"
  )

  estimate_field <- standardized_field(
    table,
    c(
      "qD",
      "qD.obs",
      "Estimate"
    ),
    "richness estimate"
  )

  lower_field <- find_field(
    names(
      table
    ),
    c(
      "qD.LCL",
      "LCL"
    ),
    required = FALSE
  )

  upper_field <- find_field(
    names(
      table
    ),
    c(
      "qD.UCL",
      "UCL"
    ),
    required = FALSE
  )

  plot_data <- table

  plot_data$latitude_band <- factor(
    plot_data[[assemblage_field]],
    levels = band_order
  )

  plot_data$estimate_value <- safe_numeric(
    plot_data[[estimate_field]]
  )

  plot_object <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = latitude_band,
      y = estimate_value
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::theme_bw(
      base_size = 12
    ) +
    ggplot2::labs(
      title = title,
      x = "Latitude band",
      y = "Estimated genus richness"
    )

  if (
    !is.na(
      lower_field
    ) &&
    !is.na(
      upper_field
    )
  ) {
    plot_data$lower_value <- safe_numeric(
      plot_data[[lower_field]]
    )

    plot_data$upper_value <- safe_numeric(
      plot_data[[upper_field]]
    )

    plot_object <- plot_object +
      ggplot2::geom_errorbar(
        data = plot_data,
        ggplot2::aes(
          ymin = lower_value,
          ymax = upper_value
        ),
        width = 0.2
      )
  }

  ggplot2::ggsave(
    filename,
    plot_object,
    width = 8.2,
    height = 6,
    dpi = 320
  )
}


sha256_file <- function(
  path
) {
  unname(
    tools::md5sum(
      path
    )
  )
}

# Note: field name retained for compatibility with earlier pipeline manifests.
# R base tools::md5sum produces MD5, so the provenance key is explicitly md5.

# ------------------------------ inputs -------------------------------------

primary_matrix_path <- first_existing(
  c(
    file.path(
      analysis_ready,
      "02_incidence_matrices_25km",
      "10_biodiversity_final_genus_by_grid25km_incidence.csv"
    ),
    file.path(
      grid_fallback,
      "10_biodiversity_final_genus_by_grid25km_incidence.csv"
    )
  ),
  label = "primary biodiversity incidence matrix"
)

ballooning_matrix_path <- first_existing(
  c(
    file.path(
      analysis_ready,
      "02_incidence_matrices_25km",
      "10_ballooning_final_genus_by_grid25km_incidence.csv"
    ),
    file.path(
      grid_fallback,
      "10_ballooning_final_genus_by_grid25km_incidence.csv"
    )
  ),
  required = FALSE,
  label = "ballooning incidence matrix"
)

taxonomy_strict_matrix_path <- first_existing(
  c(
    file.path(
      analysis_ready,
      "02_incidence_matrices_25km",
      "10_biodiversity_taxonomy_strict_genus_by_grid25km_incidence.csv"
    ),
    file.path(
      grid_fallback,
      "10_biodiversity_taxonomy_strict_genus_by_grid25km_incidence.csv"
    )
  ),
  required = FALSE,
  label = "taxonomy-strict biodiversity matrix"
)

cell_lookup_path <- first_existing(
  c(
    file.path(
      analysis_ready,
      "04_spatial_reference",
      "10_common_grid25km_cell_lookup.csv"
    ),
    file.path(
      grid_fallback,
      "10_common_grid25km_cell_lookup.csv"
    )
  ),
  label = "25-km grid-cell lookup"
)

trait_path <- first_existing(
  c(
    file.path(
      analysis_ready,
      "03_trait_tables",
      "07_reviewed_genus_trait_lookup_normalized.csv"
    ),
    file.path(
      analysis_ready,
      "03_trait_tables",
      "07_reviewed_genus_trait_lookup_final.csv"
    ),
    file.path(
      trait_fallback,
      "07_reviewed_genus_trait_lookup_final.csv"
    )
  ),
  label = "final reviewed genus trait lookup"
)

message(
  "Primary matrix: ",
  primary_matrix_path
)

message(
  "Cell lookup: ",
  cell_lookup_path
)

message(
  "Trait lookup: ",
  trait_path
)

# ------------------------------ load ---------------------------------------

primary <- read_incidence_matrix(
  primary_matrix_path
)

if (
  "Fesa" %in% primary$genera
) {
  stop(
    "Fesa remains in the primary incidence matrix."
  )
}

if (
  !is.na(
    ballooning_matrix_path
  )
) {
  ballooning_matrix <- read_incidence_matrix(
    ballooning_matrix_path
  )

  if (
    !identical(
      primary$genera,
      ballooning_matrix$genera
    ) ||
    !identical(
      primary$cells,
      ballooning_matrix$cells
    ) ||
    !identical(
      unname(
        primary$matrix
      ),
      unname(
        ballooning_matrix$matrix
      )
    )
  ) {
    stop(
      "Primary biodiversity and ballooning matrices are not identical."
    )
  }
}

cell_lookup_raw <- read.csv(
  cell_lookup_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cell_id_field <- find_field(
  names(
    cell_lookup_raw
  ),
  c(
    "grid_cell_id",
    "cell_id"
  ),
  label = "grid-cell ID"
)

band_field <- find_field(
  names(
    cell_lookup_raw
  ),
  c(
    "centroid_latitude_band",
    "latitude_band"
  ),
  label = "centroid latitude band"
)

centroid_latitude_field <- find_field(
  names(
    cell_lookup_raw
  ),
  c(
    "centroid_latitude",
    "latitude"
  ),
  label = "cell centroid latitude"
)

centroid_longitude_field <- find_field(
  names(
    cell_lookup_raw
  ),
  c(
    "centroid_longitude",
    "longitude"
  ),
  label = "cell centroid longitude"
)

cell_lookup_raw$grid_cell_id <- trim_character(
  cell_lookup_raw[[cell_id_field]]
)

cell_lookup_raw$latitude_band <- trim_character(
  cell_lookup_raw[[band_field]]
)

cell_lookup_raw$centroid_latitude_standardized <- safe_numeric(
  cell_lookup_raw[[centroid_latitude_field]]
)

cell_lookup_raw$centroid_longitude_standardized <- safe_numeric(
  cell_lookup_raw[[centroid_longitude_field]]
)

missing_cells <- setdiff(
  primary$cells,
  cell_lookup_raw$grid_cell_id
)

if (
  length(
    missing_cells
  ) > 0
) {
  stop(
    "Primary matrix cells missing from the cell lookup:\n",
    paste(
      missing_cells,
      collapse = "\n"
    )
  )
}

cell_lookup <- cell_lookup_raw[
  match(
    primary$cells,
    cell_lookup_raw$grid_cell_id
  ),
  ,
  drop = FALSE
]

if (
  any(
    !cell_lookup$latitude_band %in% band_order
  )
) {
  stop(
    "Cells assigned outside the five manuscript latitude bands:\n",
    paste(
      unique(
        cell_lookup$latitude_band[
          !cell_lookup$latitude_band %in% band_order
        ]
      ),
      collapse = "\n"
    )
  )
}

cell_lookup$latitude_band <- factor(
  cell_lookup$latitude_band,
  levels = band_order,
  ordered = TRUE
)

traits <- load_traits(
  trait_path,
  primary$genera
)

if (
  !identical(
    tolower(
      primary$genera
    ),
    tolower(
      traits$genus
    )
  )
) {
  stop(
    "Trait rows do not align with the primary incidence-matrix rows."
  )
}

# ----------------------- cell-to-band and matrices --------------------------

cell_band_output <- data.frame(
  grid_cell_id = cell_lookup$grid_cell_id,
  latitude_band = as.character(
    cell_lookup$latitude_band
  ),
  latitude_band_label = unname(
    band_labels[
      as.character(
        cell_lookup$latitude_band
      )
    ]
  ),
  centroid_latitude = cell_lookup$centroid_latitude_standardized,
  centroid_longitude = cell_lookup$centroid_longitude_standardized,
  cell_lookup,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

duplicate_names <- duplicated(
  names(
    cell_band_output
  )
)

cell_band_output <- cell_band_output[
  !duplicate_names
]

write.csv(
  cell_band_output,
  file.path(
    output_dir,
    "11_grid25km_cell_to_latitude_band.csv"
  ),
  row.names = FALSE,
  na = ""
)

primary_band_matrices <- make_band_matrices(
  primary$matrix,
  cell_lookup
)

primary_genus_by_band <- make_genus_by_band(
  primary_band_matrices
)

primary_genus_by_band_output <- data.frame(
  genus = rownames(
    primary_genus_by_band
  ),
  traits[
    ,
    c(
      "order",
      "family",
      "ballooning_status",
      "ballooning_binary",
      "analysis_class",
      "final_confidence",
      "evidence_level"
    ),
    drop = FALSE
  ],
  primary_genus_by_band,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.csv(
  primary_genus_by_band_output,
  file.path(
    output_dir,
    "11_biodiversity_final_genus_by_latitude_band_incidence.csv"
  ),
  row.names = FALSE,
  na = ""
)

for (
  band in band_order
) {
  write_matrix_csv(
    primary_band_matrices[[band]],
    file.path(
      output_dir,
      paste0(
        "11_biodiversity_final_",
        band,
        "_genus_by_grid25km_incidence.csv"
      )
    )
  )
}

# ------------------------- iNEXT: primary ----------------------------------

primary_inext <- run_inext_analysis(
  dataset_key = "biodiversity_final",
  band_matrices = primary_band_matrices,
  output_dir = output_dir,
  nboot = nboot
)

primary_sampling <- band_sampling_summary(
  band_matrices = primary_band_matrices,
  cell_lookup = cell_lookup,
  data_info = primary_inext$data_info,
  dataset_key = "biodiversity_final"
)

write.csv(
  primary_sampling,
  file.path(
    output_dir,
    "11_biodiversity_final_sampling_summary_by_band.csv"
  ),
  row.names = FALSE,
  na = ""
)

plot_standardized_richness(
  primary_inext$equal_cell,
  "Equal-cell standardized genus richness",
  file.path(
    figure_dir,
    "11_biodiversity_final_equal_cell_richness.png"
  )
)

plot_standardized_richness(
  primary_inext$coverage_standardized,
  "Coverage-standardized genus richness",
  file.path(
    figure_dir,
    "11_biodiversity_final_coverage_standardized_richness.png"
  )
)

# ------------------------- turnover: primary -------------------------------

primary_beta <- write_beta_outputs(
  dataset_key = "biodiversity_final",
  genus_by_band = primary_genus_by_band,
  output_dir = output_dir
)

# ------------------------- ballooning composition ---------------------------

cell_richness <- cell_richness_table(
  incidence_matrix = primary$matrix,
  cell_lookup = cell_lookup,
  traits = traits
)

write.csv(
  cell_richness,
  file.path(
    output_dir,
    "11_grid25km_cell_richness_and_ballooning_composition.csv"
  ),
  row.names = FALSE,
  na = ""
)

ballooning_summary <- ballooning_band_summary(
  genus_by_band = primary_genus_by_band,
  band_matrices = primary_band_matrices,
  traits = traits
)

write.csv(
  ballooning_summary,
  file.path(
    output_dir,
    "11_ballooning_composition_by_latitude_band.csv"
  ),
  row.names = FALSE,
  na = ""
)

cell_boxplot <- ggplot2::ggplot(
  cell_richness,
  ggplot2::aes(
    x = factor(
      latitude_band,
      levels = band_order
    ),
    y = total_genus_richness
  )
) +
  ggplot2::geom_boxplot() +
  ggplot2::geom_jitter(
    width = 0.15,
    alpha = 0.45,
    size = 1.3
  ) +
  ggplot2::theme_bw(
    base_size = 12
  ) +
  ggplot2::labs(
    title = "Observed genus richness among occupied 25-km cells",
    x = "Latitude band",
    y = "Genus richness per occupied cell"
  )

ballooning_plot <- ggplot2::ggplot(
  ballooning_summary,
  ggplot2::aes(
    x = factor(
      latitude_band,
      levels = band_order
    ),
    y = ballooning_genus_proportion
  )
) +
  ggplot2::geom_col() +
  ggplot2::geom_point(
    ggplot2::aes(
      y = confidence_restricted_ballooning_proportion
    ),
    size = 2.7
  ) +
  ggplot2::scale_y_continuous(
    limits = c(
      0,
      1
    )
  ) +
  ggplot2::theme_bw(
    base_size = 12
  ) +
  ggplot2::labs(
    title = "Ballooning-capable genera by latitude band",
    x = "Latitude band",
    y = "Proportion of observed genera",
    caption = "Bars: all completed traits; points: LOW-confidence genera excluded"
  )

ggplot2::ggsave(
  file.path(
    figure_dir,
    "11_cell_genus_richness_by_latitude_band.png"
  ),
  cell_boxplot,
  width = 8.2,
  height = 6,
  dpi = 320
)

ggplot2::ggsave(
  file.path(
    figure_dir,
    "11_ballooning_proportion_by_latitude_band.png"
  ),
  ballooning_plot,
  width = 8.2,
  height = 6,
  dpi = 320
)

# --------------------- taxonomy-strict sensitivity --------------------------

taxonomy_results <- NULL

if (
  !is.na(
    taxonomy_strict_matrix_path
  )
) {
  message(
    "Running taxonomy-strict sensitivity analysis."
  )

  taxonomy <- read_incidence_matrix(
    taxonomy_strict_matrix_path
  )

  if (
    !identical(
      primary$genera,
      taxonomy$genera
    ) ||
    !identical(
      primary$cells,
      taxonomy$cells
    )
  ) {
    stop(
      "Taxonomy-strict matrix is not aligned to the primary genus/cell universe."
    )
  }

  taxonomy_band_matrices <- make_band_matrices(
    taxonomy$matrix,
    cell_lookup
  )

  taxonomy_genus_by_band <- make_genus_by_band(
    taxonomy_band_matrices
  )

  taxonomy_genus_by_band_output <- data.frame(
    genus = rownames(
      taxonomy_genus_by_band
    ),
    taxonomy_genus_by_band,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  write.csv(
    taxonomy_genus_by_band_output,
    file.path(
      output_dir,
      "11_taxonomy_strict_genus_by_latitude_band_incidence.csv"
    ),
    row.names = FALSE,
    na = ""
  )

  taxonomy_inext <- run_inext_analysis(
    dataset_key = "taxonomy_strict",
    band_matrices = taxonomy_band_matrices,
    output_dir = output_dir,
    nboot = nboot
  )

  taxonomy_sampling <- band_sampling_summary(
    band_matrices = taxonomy_band_matrices,
    cell_lookup = cell_lookup,
    data_info = taxonomy_inext$data_info,
    dataset_key = "taxonomy_strict"
  )

  write.csv(
    taxonomy_sampling,
    file.path(
      output_dir,
      "11_taxonomy_strict_sampling_summary_by_band.csv"
    ),
    row.names = FALSE,
    na = ""
  )

  taxonomy_beta <- write_beta_outputs(
    dataset_key = "taxonomy_strict",
    genus_by_band = taxonomy_genus_by_band,
    output_dir = output_dir
  )

  richness_comparison <- merge(
    primary_sampling[
      ,
      c(
        "latitude_band",
        "observed_genus_richness"
      )
    ],
    taxonomy_sampling[
      ,
      c(
        "latitude_band",
        "observed_genus_richness"
      )
    ],
    by = "latitude_band",
    suffixes = c(
      "_primary",
      "_taxonomy_strict"
    ),
    all = TRUE,
    sort = FALSE
  )

  richness_comparison$richness_difference <- (
    richness_comparison$observed_genus_richness_taxonomy_strict -
      richness_comparison$observed_genus_richness_primary
  )

  write.csv(
    richness_comparison,
    file.path(
      output_dir,
      "11_primary_vs_taxonomy_strict_richness_comparison.csv"
    ),
    row.names = FALSE,
    na = ""
  )

  beta_comparison <- merge(
    primary_beta$pairwise[
      ,
      c(
        "band_1",
        "band_2",
        "jaccard_dissimilarity",
        "simpson_turnover",
        "sorensen_nestedness_resultant"
      )
    ],
    taxonomy_beta$pairwise[
      ,
      c(
        "band_1",
        "band_2",
        "jaccard_dissimilarity",
        "simpson_turnover",
        "sorensen_nestedness_resultant"
      )
    ],
    by = c(
      "band_1",
      "band_2"
    ),
    suffixes = c(
      "_primary",
      "_taxonomy_strict"
    ),
    all = TRUE,
    sort = FALSE
  )

  write.csv(
    beta_comparison,
    file.path(
      output_dir,
      "11_primary_vs_taxonomy_strict_beta_comparison.csv"
    ),
    row.names = FALSE,
    na = ""
  )

  taxonomy_results <- list(
    matrix_path = taxonomy_strict_matrix_path,
    sampling = taxonomy_sampling,
    inext = taxonomy_inext,
    beta = taxonomy_beta,
    genus_by_band = taxonomy_genus_by_band
  )
} else {
  message(
    "Taxonomy-strict matrix was not found; skipping that sensitivity analysis."
  )
}

# ------------------------------ validation ---------------------------------

cell_counts_by_band <- table(
  factor(
    cell_lookup$latitude_band,
    levels = band_order
  )
)

validation <- data.frame(
  check = c(
    "primary_matrix_binary",
    "primary_genera_unique",
    "primary_cells_unique",
    "all_primary_cells_in_lookup",
    "all_cells_in_manuscript_bands",
    "all_primary_genera_have_traits",
    "all_primary_genera_have_explicit_C3_N0_D4_trait",
    "fesa_absent",
    "band_cell_counts_sum_to_total",
    "primary_genus_band_rows_match_primary_genera",
    "primary_genus_band_columns_match_five_bands"
  ),
  passed = c(
    all(
      primary$matrix %in% c(
        0,
        1
      )
    ),
    !anyDuplicated(
      tolower(
        primary$genera
      )
    ),
    !anyDuplicated(
      primary$cells
    ),
    length(
      missing_cells
    ) == 0,
    all(
      as.character(
        cell_lookup$latitude_band
      ) %in% band_order
    ),
    nrow(
      traits
    ) == nrow(
      primary$matrix
    ),
    all(
      traits$ballooning_binary %in% c(
        -1,
        0,
        1
      )
    ),
    !"Fesa" %in% primary$genera,
    sum(
      cell_counts_by_band
    ) == ncol(
      primary$matrix
    ),
    nrow(
      primary_genus_by_band
    ) == nrow(
      primary$matrix
    ),
    ncol(
      primary_genus_by_band
    ) == length(
      band_order
    )
  ),
  stringsAsFactors = FALSE
)

if (
  any(
    !validation$passed
  )
) {
  write.csv(
    validation,
    file.path(
      output_dir,
      "11_validation.csv"
    ),
    row.names = FALSE
  )

  stop(
    "One or more Step 11 validation checks failed."
  )
}

write.csv(
  validation,
  file.path(
    output_dir,
    "11_validation.csv"
  ),
  row.names = FALSE
)

trait_confidence_summary <- as.data.frame(
  table(
    traits$final_confidence,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)

names(
  trait_confidence_summary
) <- c(
  "trait_confidence",
  "genera"
)

write.csv(
  trait_confidence_summary,
  file.path(
    output_dir,
    "11_trait_confidence_summary.csv"
  ),
  row.names = FALSE,
  na = ""
)

# ------------------------------ provenance ---------------------------------

input_files <- c(
  primary_matrix_path,
  cell_lookup_path,
  trait_path
)

if (
  !is.na(
    ballooning_matrix_path
  )
) {
  input_files <- c(
    input_files,
    ballooning_matrix_path
  )
}

if (
  !is.na(
    taxonomy_strict_matrix_path
  )
) {
  input_files <- c(
    input_files,
    taxonomy_strict_matrix_path
  )
}

input_manifest <- data.frame(
  path = input_files,
  md5 = unname(
    tools::md5sum(
      input_files
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  input_manifest,
  file.path(
    output_dir,
    "11_input_file_manifest.csv"
  ),
  row.names = FALSE
)

output_files <- list.files(
  output_dir,
  recursive = TRUE,
  full.names = TRUE
)

output_files <- output_files[
  file.info(
    output_files
  )$isdir == FALSE
]

output_manifest <- data.frame(
  path = output_files,
  relative_path = substring(
    output_files,
    nchar(
      output_dir
    ) + 2
  ),
  bytes = file.info(
    output_files
  )$size,
  md5 = unname(
    tools::md5sum(
      output_files
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  output_manifest,
  file.path(
    output_dir,
    "11_output_file_manifest.csv"
  ),
  row.names = FALSE
)

provenance <- list(
  created_utc = format(
    Sys.time(),
    tz = "UTC",
    usetz = TRUE
  ),
  script = "11_latitude_band_diversity_turnover.R",
  project_root = normalizePath(
    project_root,
    winslash = "/",
    mustWork = TRUE
  ),
  random_seed = 20260713,
  inext_bootstrap_replicates = nboot,
  latitude_band_order = band_order,
  latitude_band_labels = unname(
    band_labels
  ),
  sampling_unit = "occupied 25-km equal-area grid cell",
  primary_trait_definition = "C3 = D1 + D2 + D3 versus fixed N0; D4 excluded",
  trait_evidence_field = unique(traits$evidence_field),
  trait_class_counts = as.list(table(traits$analysis_class)),
  incidence_definition = (
    "A genus is present in a latitude band when it occurs in at least one "
  ),
  incidence_definition_continued = (
    "occupied 25-km cell assigned to that band by cell centroid."
  ),
  primary_dimensions = list(
    genera = nrow(
      primary$matrix
    ),
    occupied_cells = ncol(
      primary$matrix
    ),
    genus_cell_presences = sum(
      primary$matrix
    )
  ),
  occupied_cells_by_band = as.list(
    as.integer(
      cell_counts_by_band
    )
  ),
  minimum_cells_for_equal_cell_standardization = primary_inext$minimum_cells,
  confidence_sensitivity_rule = (
    "Exclude genera explicitly classified LOW confidence; retain HIGH, MEDIUM, "
  ),
  confidence_sensitivity_rule_continued = (
    "and legacy resolved traits whose confidence field is unspecified."
  ),
  beta_diversity = list(
    jaccard = "(b+c)/(a+b+c)",
    sorensen = "(b+c)/(2a+b+c)",
    simpson_turnover = "min(b,c)/(a+min(b,c))",
    sorensen_nestedness_resultant = "sorensen - simpson_turnover"
  ),
  taxonomy_strict_analysis_completed = !is.null(
    taxonomy_results
  ),
  fesa_excluded = TRUE,
  validation_passed = all(
    validation$passed
  ),
  input_manifest = input_manifest
)

jsonlite::write_json(
  provenance,
  file.path(
    output_dir,
    "11_provenance.json"
  ),
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "11_R_session_info.txt"
  )
)

# Copy the running script into the output folder when possible.
script_argument <- grep(
  "^--file=",
  commandArgs(
    trailingOnly = FALSE
  ),
  value = TRUE
)

if (
  length(
    script_argument
  ) == 1
) {
  running_script <- sub(
    "^--file=",
    "",
    script_argument
  )

  if (
    file.exists(
      running_script
    )
  ) {
    file.copy(
      running_script,
      file.path(
        output_dir,
        "11_latitude_band_diversity_turnover.R"
      ),
      overwrite = TRUE
    )
  }
}

# ------------------------------ finish -------------------------------------

message(
  "\n",
  strrep(
    "=",
    78
  )
)

message(
  "STEP 11 COMPLETED SUCCESSFULLY"
)

message(
  strrep(
    "=",
    78
  )
)

message(
  "Primary genera: ",
  nrow(
    primary$matrix
  )
)

message(
  "Occupied 25-km cells: ",
  ncol(
    primary$matrix
  )
)

message(
  "Genus-cell incidences: ",
  sum(
    primary$matrix
  )
)

message(
  "Cells by latitude band:"
)

for (
  band in band_order
) {
  message(
    "  ",
    band,
    ": ",
    cell_counts_by_band[[band]]
  )
}

message(
  "Equal-cell standardization level: ",
  primary_inext$minimum_cells,
  " occupied cells per band"
)

message(
  "Taxonomy-strict sensitivity completed: ",
  !is.null(
    taxonomy_results
  )
)

message(
  "Outputs:\n",
  output_dir
)
