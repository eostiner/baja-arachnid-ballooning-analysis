#!/usr/bin/env Rscript
packages <- c(
  "gdm", "vegan", "dplyr", "tidyr", "purrr", "readr", "tibble",
  "ggplot2", "patchwork", "scales", "svglite"
)
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (!length(missing)) {
  cat("All Step 12J packages are already installed.\n")
  quit(status = 0)
}
cat("Installing:", paste(missing, collapse = ", "), "\n")
install.packages(missing, repos = "https://cloud.r-project.org", dependencies = TRUE)
still_missing <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing)) stop("Packages still missing after installation: ", paste(still_missing, collapse = ", "))
cat("All Step 12J packages installed successfully.\n")
