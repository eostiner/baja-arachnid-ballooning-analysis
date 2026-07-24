#!/usr/bin/env Rscript

# Coverage-standardized iNEXT diagnostic for the Baja ballooning publication QC.
# This is intentionally separate from the Baselga replacement-nestedness analysis.
# It consumes 07_iNEXT_INCIDENCE_FREQUENCY_INPUT.csv produced by the Python core.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 02_run_iNEXT_hill.R <QC output directory> [nboot] [seed]")
}
out_dir <- normalizePath(path.expand(args[[1]]), mustWork = TRUE)
nboot <- if (length(args) >= 2) as.integer(args[[2]]) else 200L
seed <- if (length(args) >= 3) as.integer(args[[3]]) else 20260723L
if (is.na(nboot) || nboot < 20L) stop("nboot must be at least 20")
if (is.na(seed)) stop("seed must be an integer")
set.seed(seed)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

needed <- c("iNEXT", "ggplot2")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  message("Installing missing R packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(iNEXT)
  library(ggplot2)
})

infile <- file.path(out_dir, "07_iNEXT_INCIDENCE_FREQUENCY_INPUT.csv")
if (!file.exists(infile)) stop("Missing input: ", infile)
dat <- read.csv(infile, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("scope", "assemblage", "assemblage_label", "trait_class", "trait_label", "n_sampling_units", "genus", "incidence_frequency")
missing_cols <- setdiff(required, names(dat))
if (length(missing_cols) > 0) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

# Build incidence-frequency vectors. First entry is number of sampling units.
dat$key <- paste(dat$scope, dat$assemblage, dat$trait_class, sep = "||")
keys <- unique(dat$key)
assemblage_list <- lapply(keys, function(k) {
  d <- dat[dat$key == k, , drop = FALSE]
  c(as.integer(d$n_sampling_units[[1]]), as.integer(d$incidence_frequency))
})
names(assemblage_list) <- keys

# Keep key metadata for parsing results.
meta <- unique(dat[c("key", "scope", "assemblage", "assemblage_label", "trait_class", "trait_label", "n_sampling_units")])

message("Running iNEXT for ", length(assemblage_list), " assemblages with nboot=", nboot, "; seed=", seed)
fit_rds <- file.path(out_dir, paste0("08_iNEXT_FIT_nboot", nboot, "_seed", seed, ".rds"))
if (file.exists(fit_rds)) {
  message("Resuming from saved iNEXT fit: ", fit_rds)
  fit <- readRDS(fit_rds)
} else {
  fit <- iNEXT(
    assemblage_list,
    q = c(0, 1, 2),
    datatype = "incidence_freq",
    knots = 40,
    se = TRUE,
    conf = 0.95,
    nboot = nboot
  )
  saveRDS(fit, fit_rds)
}

# DataInfo gives observed sample coverage for each assemblage.
info <- DataInfo(assemblage_list, datatype = "incidence_freq")
write.csv(info, file.path(out_dir, "08_iNEXT_DATA_INFO.csv"), row.names = FALSE)

coverage_candidates <- grep("^(SC|Coverage|Sample.coverage)$", names(info), value = TRUE, ignore.case = TRUE)
if (length(coverage_candidates) == 0) {
  coverage_candidates <- grep("cover", names(info), value = TRUE, ignore.case = TRUE)
}
if (length(coverage_candidates) == 0) {
  warning("Could not identify observed sample-coverage column in DataInfo; using 0.90 as conservative target.")
  target_coverage <- 0.90
} else {
  target_coverage <- min(as.numeric(info[[coverage_candidates[[1]]]]), na.rm = TRUE)
  target_coverage <- max(0.50, min(0.999, target_coverage))
}
message("Shared coverage target: ", format(target_coverage, digits = 5))

est <- estimateD(
  assemblage_list,
  q = c(0, 1, 2),
  datatype = "incidence_freq",
  base = "coverage",
  level = target_coverage,
  conf = 0.95,
  nboot = nboot
)

# Normalize the assemblage-name column across iNEXT versions.
name_col <- intersect(c("Assemblage", "assemblage", "site", "Site"), names(est))
if (length(name_col) == 0) stop("Could not identify assemblage column in estimateD output")
names(est)[names(est) == name_col[[1]]] <- "key"
est <- merge(est, meta, by = "key", all.x = TRUE, sort = FALSE)
est$target_coverage <- target_coverage
write.csv(est, file.path(out_dir, "09_iNEXT_COVERAGE_STANDARDIZED_HILL.csv"), row.names = FALSE)

# Extract a stable diversity-estimate column name.
d_col <- intersect(c("qD", "qD.base", "Diversity"), names(est))
q_col <- intersect(c("Order.q", "q", "Order"), names(est))
if (length(d_col) > 0 && length(q_col) > 0) {
  names(est)[names(est) == d_col[[1]]] <- "qD_value"
  names(est)[names(est) == q_col[[1]]] <- "q_order"

  # Derive a transparent adjacent-pair Hill beta diagnostic from equally weighted,
  # coverage-standardized band estimates and the corresponding pooled-pair gamma.
  band <- est[est$scope == "band", , drop = FALSE]
  pool <- est[est$scope == "adjacent_pair_pool", , drop = FALSE]
  pair_names <- c("23-24N__24-26N", "24-26N__26-28N", "26-28N__28-30N", "28-30N__30-32N")
  beta_rows <- list()
  z <- 1L
  for (pair in pair_names) {
    parts <- strsplit(pair, "__", fixed = TRUE)[[1]]
    for (trait in unique(dat$trait_class)) {
      for (qq in sort(unique(est$q_order))) {
        d1 <- band$qD_value[band$assemblage == parts[[1]] & band$trait_class == trait & band$q_order == qq]
        d2 <- band$qD_value[band$assemblage == parts[[2]] & band$trait_class == trait & band$q_order == qq]
        gamma <- pool$qD_value[pool$assemblage == pair & pool$trait_class == trait & pool$q_order == qq]
        if (length(d1) != 1 || length(d2) != 1 || length(gamma) != 1) next
        if (qq == 0) {
          alpha <- mean(c(d1, d2))
        } else if (qq == 1) {
          alpha <- sqrt(d1 * d2)
        } else {
          alpha <- (mean(c(d1^(1 - qq), d2^(1 - qq))))^(1 / (1 - qq))
        }
        beta_rows[[z]] <- data.frame(
          pair = pair,
          band_1 = parts[[1]],
          band_2 = parts[[2]],
          trait_class = trait,
          q_order = qq,
          alpha_equal_weight = alpha,
          gamma_pooled = gamma,
          hill_beta_gamma_over_alpha = gamma / alpha,
          target_coverage = target_coverage,
          stringsAsFactors = FALSE
        )
        z <- z + 1L
      }
    }
  }
  if (length(beta_rows) > 0) {
    beta <- do.call(rbind, beta_rows)
    write.csv(beta, file.path(out_dir, "10_iNEXT_DERIVED_ADJACENT_HILL_BETA.csv"), row.names = FALSE)
  }

  caption_lines <- c(
    "iNEXT caption addendum (values inserted from the current run)",
    "",
    paste0("Incidence-based Hill diversity was estimated for q = 0, 1, and 2 and standardized to a shared sample coverage of ", round(target_coverage, 3), ". These estimates provide a coverage-standardized diversity diagnostic and are reported separately from the Baselga replacement-nestedness partition."),
    "",
    "Coverage-standardized band estimates:"
  )
  band_sorted <- band[order(match(band$assemblage, c("23-24N", "24-26N", "26-28N", "28-30N", "30-32N")), band$trait_class, band$q_order), , drop = FALSE]
  for (i in seq_len(nrow(band_sorted))) {
    caption_lines <- c(caption_lines, paste0(
      band_sorted$assemblage_label[[i]], "; ", band_sorted$trait_label[[i]],
      "; q=", band_sorted$q_order[[i]], ": ", round(band_sorted$qD_value[[i]], 3)
    ))
  }
  writeLines(caption_lines, file.path(out_dir, "PUBLICATION_CAPTION_iNEXT_ADDENDUM.txt"))
}

# Export raw iNEXT estimates for audit.
if (!is.null(fit$iNextEst)) {
  for (nm in names(fit$iNextEst)) {
    obj <- fit$iNextEst[[nm]]
    if (is.data.frame(obj)) {
      safe <- gsub("[^A-Za-z0-9]+", "_", nm)
      write.csv(obj, file.path(out_dir, paste0("08_iNEXT_RAW_", safe, ".csv")), row.names = FALSE)
    }
  }
}

# Curves. iNEXT 3.x uses the output-column names "Order.q" and
# "Assemblage" for faceting and colour. Older versions accepted "order" and
# "site". Try the current API first, then the legacy API, and finally a
# no-facet fallback so a plotting-version difference cannot invalidate the
# completed diversity calculations.
make_ggiNEXT <- function(fit_object, plot_type) {
  errors <- character(0)

  current <- tryCatch(
    ggiNEXT(fit_object, type = plot_type, facet.var = "Order.q", color.var = "Assemblage"),
    error = function(e) {
      errors <<- c(errors, paste0("current API: ", conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(current)) return(current)

  legacy <- tryCatch(
    ggiNEXT(fit_object, type = plot_type, facet.var = "order", color.var = "site"),
    error = function(e) {
      errors <<- c(errors, paste0("legacy API: ", conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(legacy)) return(legacy)

  fallback <- tryCatch(
    ggiNEXT(fit_object, type = plot_type),
    error = function(e) {
      errors <<- c(errors, paste0("default API: ", conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(fallback)) return(fallback)

  stop("Unable to create ggiNEXT plot. ", paste(errors, collapse = " | "))
}

# Use ASCII-safe titles to avoid font warnings.
p1 <- make_ggiNEXT(fit, 1) +
  labs(title = "iNEXT sample-size-based Hill diversity", subtitle = "Band and adjacent-pair incidence-frequency assemblages", x = "Number of 25-km sampling units", y = "Hill diversity") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
p2 <- make_ggiNEXT(fit, 3) +
  labs(title = "iNEXT coverage-based Hill diversity", subtitle = paste0("Shared coverage target = ", round(target_coverage, 3)), x = "Sample coverage", y = "Hill diversity") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

ggsave(file.path(out_dir, "figures", "Figure_QC_04_iNEXT_sample_size_curves.png"), p1, width = 14, height = 9, dpi = 300)
ggsave(file.path(out_dir, "figures", "Figure_QC_04_iNEXT_sample_size_curves.pdf"), p1, width = 14, height = 9)
ggsave(file.path(out_dir, "figures", "Figure_QC_05_iNEXT_coverage_curves.png"), p2, width = 14, height = 9, dpi = 300)
ggsave(file.path(out_dir, "figures", "Figure_QC_05_iNEXT_coverage_curves.pdf"), p2, width = 14, height = 9)

writeLines(c(
  "iNEXT stage completed.",
  paste0("nboot=", nboot),
  paste0("seed=", seed),
  paste0("shared_coverage_target=", target_coverage),
  "Important: these Hill-diversity outputs are coverage-standardized diagnostics and are not equivalent to the Baselga replacement-nestedness partition."
), file.path(out_dir, "iNEXT_RUN_SUMMARY.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "iNEXT_R_SESSION_INFO.txt"))
message("iNEXT outputs written to: ", out_dir)
