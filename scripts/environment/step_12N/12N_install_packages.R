#!/usr/bin/env Rscript

required <- c("ade4", "ggplot2", "sf", "patchwork", "jsonlite", "svglite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (!length(missing)) {
  cat("All Step 12N packages are already installed.\n")
  quit(status = 0)
}

cat("Installing missing Step 12N packages:", paste(missing, collapse = ", "), "\n")
install.packages(missing, repos = "https://cloud.r-project.org")

still_missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing)) {
  stop("Packages still missing after installation: ", paste(still_missing, collapse = ", "))
}
cat("Step 12N package installation complete.\n")
