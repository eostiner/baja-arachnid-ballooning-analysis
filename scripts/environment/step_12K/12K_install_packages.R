#!/usr/bin/env Rscript

packages <- c("mgcv", "ggplot2", "patchwork")
repos <- "https://cloud.r-project.org"

installed <- rownames(installed.packages())
needed <- setdiff(packages, installed)

if (length(needed)) {
  install.packages(needed, repos = repos, dependencies = TRUE)
} else {
  cat("All Step 12K packages are already installed.\n")
}

for (pkg in packages) {
  cat(pkg, ": ", as.character(packageVersion(pkg)), "\n", sep = "")
}
