#!/usr/bin/env Rscript

# Base-R publication figures for the Baja Ballooning QC package.
# This script deliberately uses no add-on packages, so figures are still
# produced when Python matplotlib is unavailable.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 03_make_publication_figures.R <QC output directory>")
out_dir <- normalizePath(path.expand(args[[1]]), mustWork = TRUE)
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_required <- function(name) {
  path <- file.path(out_dir, name)
  if (!file.exists(path)) stop("Missing required file: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

rich <- read_required("04_EQUAL_CELL_RICHNESS_SUMMARY.csv")
beta <- read_required("05_BASELGA_PAIRWISE_SUMMARY.csv")
contrast <- read_required("06_TRAIT_CONTRAST_SUMMARY.csv")

as_adjacent <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}
beta$adjacent <- as_adjacent(beta$adjacent)
contrast$adjacent <- as_adjacent(contrast$adjacent)

band_order <- c("23-24N", "24-26N", "26-28N", "28-30N", "30-32N")
band_labels <- c("23–24°N", "24–26°N", "26–28°N", "28–30°N", "30–32°N")
pair_order <- data.frame(
  band_1 = c("23-24N", "24-26N", "26-28N", "28-30N"),
  band_2 = c("24-26N", "26-28N", "28-30N", "30-32N"),
  label = c("23–24 / 24–26°N", "24–26 / 26–28°N", "26–28 / 28–30°N", "28–30 / 30–32°N"),
  stringsAsFactors = FALSE
)
trait_order <- c("ballooning", "non_ballooning")
trait_labels <- c(ballooning = "Ballooning-capable", non_ballooning = "Non-ballooning")

find_row <- function(dat, ...) {
  filters <- list(...)
  keep <- rep(TRUE, nrow(dat))
  for (nm in names(filters)) keep <- keep & dat[[nm]] == filters[[nm]]
  rows <- dat[keep, , drop = FALSE]
  if (nrow(rows) != 1) stop("Expected one row but found ", nrow(rows), " for ", paste(names(filters), unlist(filters), collapse = "; "))
  rows
}

error_bar <- function(x, y, lo, hi, length = 0.04, ...) {
  arrows(x, lo, x, hi, angle = 90, code = 3, length = length, ...)
  points(x, y, ...)
}

save_triplet <- function(stem, width, height, plot_fun) {
  png(file.path(fig_dir, paste0(stem, ".png")), width = width, height = height, units = "in", res = 300)
  plot_fun(); dev.off()
  pdf(file.path(fig_dir, paste0(stem, ".pdf")), width = width, height = height, useDingbats = FALSE)
  plot_fun(); dev.off()
  svg(file.path(fig_dir, paste0(stem, ".svg")), width = width, height = height)
  plot_fun(); dev.off()
}

plot_main <- function() {
  old <- par(mar = c(7.2, 5.1, 2.0, 1.2), las = 1, xpd = NA)
  on.exit(par(old))
  metrics <- c("jaccard_total", "jaccard_turnover", "jaccard_nestedness")
  labels <- c("Total Jaccard", "Replacement", "Nestedness-resultant")
  pch <- c(16, 15, 17)
  offsets <- c(-0.18, 0, 0.18)
  ylim <- range(c(contrast$p025[contrast$adjacent & contrast$metric %in% metrics], contrast$p975[contrast$adjacent & contrast$metric %in% metrics]), na.rm = TRUE)
  ylim <- c(min(-0.65, ylim[1] - 0.03), max(0.65, ylim[2] + 0.03))
  plot(NA, xlim = c(0.55, 4.45), ylim = ylim, xaxt = "n",
       xlab = "Adjacent latitude-band comparison",
       ylab = "Ballooning-capable minus non-ballooning")
  abline(h = 0, lty = 2, lwd = 1.1)
  axis(1, at = 1:4, labels = pair_order$label, las = 2)
  for (m in seq_along(metrics)) {
    for (i in seq_len(nrow(pair_order))) {
      r <- find_row(contrast, band_1 = pair_order$band_1[i], band_2 = pair_order$band_2[i], metric = metrics[m])
      error_bar(i + offsets[m], r$median, r$p025, r$p975, pch = pch[m], cex = 1.05, lwd = 1.25)
    }
  }
  legend("top", inset = c(0, -0.12), legend = labels, pch = pch, horiz = TRUE, bty = "n", xpd = NA)
  mtext("Trait contrast in Jaccard dissimilarity and its Baselga partition", side = 3, line = 0.5, font = 2)
}

plot_partition <- function() {
  old <- par(mfrow = c(2, 1), mar = c(5.5, 4.8, 2.7, 1.2), las = 1)
  on.exit(par(old))
  for (trait in trait_order) {
    replacement <- nested <- numeric(4)
    for (i in seq_len(nrow(pair_order))) {
      replacement[i] <- find_row(beta, band_1 = pair_order$band_1[i], band_2 = pair_order$band_2[i], trait_class = trait, metric = "jaccard_turnover")$median
      nested[i] <- find_row(beta, band_1 = pair_order$band_1[i], band_2 = pair_order$band_2[i], trait_class = trait, metric = "jaccard_nestedness")$median
    }
    mat <- rbind(replacement, nested)
    barplot(mat, beside = FALSE, names.arg = pair_order$label, las = 2, ylim = c(0, 1),
            ylab = "Jaccard dissimilarity", main = trait_labels[[trait]], border = NA)
    legend("topleft", legend = c("Replacement", "Nestedness-resultant"), fill = gray.colors(2, start = 0.35, end = 0.75), bty = "n")
  }
}

plot_richness_jaccard <- function() {
  old <- par(mfrow = c(2, 1), mar = c(5.7, 4.8, 3.0, 1.2), las = 1)
  on.exit(par(old))
  pchs <- c(ballooning = 16, non_ballooning = 15)
  all_y <- rich$median[rich$metric == "richness"]
  plot(1:5, rep(NA, 5), xlim = c(0.8, 5.2), ylim = range(c(rich$p025, rich$p975), na.rm = TRUE), xaxt = "n",
       xlab = "Latitude band", ylab = "Genus richness", main = "A. Equal-cell genus richness")
  axis(1, at = 1:5, labels = band_labels)
  for (trait in trait_order) {
    med <- lo <- hi <- numeric(5)
    for (i in seq_along(band_order)) {
      r <- find_row(rich, latitude_band = band_order[i], trait_class = trait, metric = "richness")
      med[i] <- r$median; lo[i] <- r$p025; hi[i] <- r$p975
    }
    lines(1:5, med, type = "b", pch = pchs[[trait]], lwd = 1.2)
    arrows(1:5, lo, 1:5, hi, angle = 90, code = 3, length = 0.04)
  }
  legend("topright", legend = unname(trait_labels[trait_order]), pch = unname(pchs[trait_order]), lty = 1, bty = "n")

  plot(1:4, rep(NA, 4), xlim = c(0.8, 4.2), ylim = c(0, 1), xaxt = "n",
       xlab = "Adjacent latitude bands", ylab = "Total Jaccard dissimilarity",
       main = "B. Adjacent-band total compositional dissimilarity")
  axis(1, at = 1:4, labels = pair_order$label, las = 2)
  for (trait in trait_order) {
    med <- lo <- hi <- numeric(4)
    for (i in seq_len(nrow(pair_order))) {
      r <- find_row(beta, band_1 = pair_order$band_1[i], band_2 = pair_order$band_2[i], trait_class = trait, metric = "jaccard_total")
      med[i] <- r$median; lo[i] <- r$p025; hi[i] <- r$p975
    }
    lines(1:4, med, type = "b", pch = pchs[[trait]], lwd = 1.2)
    arrows(1:4, lo, 1:4, hi, angle = 90, code = 3, length = 0.04)
  }
  legend("topleft", legend = unname(trait_labels[trait_order]), pch = unname(pchs[trait_order]), lty = 1, bty = "n")
}

plot_all_contrasts <- function() {
  old <- par(mar = c(7.2, 5.2, 3.0, 1.2), las = 1, xpd = NA)
  on.exit(par(old))
  metrics <- c("jaccard_total", "jaccard_turnover", "jaccard_nestedness", "simpson_replacement", "sorensen_nestedness")
  labels <- c("Jaccard total", "Jaccard replacement", "Jaccard nestedness", "Simpson replacement", "Sørensen nestedness")
  pch <- c(16, 15, 17, 1, 2)
  offsets <- seq(-0.28, 0.28, length.out = length(metrics))
  d <- contrast[contrast$adjacent & contrast$metric %in% metrics, , drop = FALSE]
  ylim <- range(c(d$p025, d$p975), na.rm = TRUE)
  plot(NA, xlim = c(0.5, 4.5), ylim = c(min(-0.7, ylim[1]), max(0.7, ylim[2])), xaxt = "n",
       xlab = "Adjacent latitude-band comparison", ylab = "Ballooning-capable minus non-ballooning",
       main = "Paired trait contrasts from identical cell draws")
  axis(1, at = 1:4, labels = pair_order$label, las = 2)
  abline(h = 0, lty = 2)
  for (m in seq_along(metrics)) {
    for (i in seq_len(nrow(pair_order))) {
      r <- find_row(contrast, band_1 = pair_order$band_1[i], band_2 = pair_order$band_2[i], metric = metrics[m])
      error_bar(i + offsets[m], r$median, r$p025, r$p975, pch = pch[m], cex = 0.9, lwd = 1.1)
    }
  }
  legend("top", inset = c(0, -0.13), legend = labels, pch = pch, ncol = 3, bty = "n", xpd = NA)
}

save_triplet("Figure_MAIN_Jaccard_Baselga_trait_contrasts", 10.2, 5.8, plot_main)
save_triplet("Figure_QC_01_richness_and_Jaccard", 11.5, 10.0, plot_richness_jaccard)
save_triplet("Figure_QC_02_Baselga_Jaccard_partition", 12.5, 10.0, plot_partition)
save_triplet("Figure_QC_03_trait_contrasts", 12.5, 7.0, plot_all_contrasts)

marker <- file.path(out_dir, "FIGURES_NOT_CREATED.txt")
if (file.exists(marker)) unlink(marker)
writeLines(c(
  "Base-R publication figures completed.",
  paste0("output_directory=", fig_dir),
  "formats=PNG,PDF,SVG",
  "Figure_MAIN_Jaccard_Baselga_trait_contrasts is the recommended main beta-diversity figure."
), file.path(out_dir, "FIGURE_RUN_SUMMARY.txt"))
message("Publication figures written to: ", fig_dir)
