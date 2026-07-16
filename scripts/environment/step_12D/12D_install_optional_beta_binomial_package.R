#!/usr/bin/env Rscript
if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  install.packages("glmmTMB", repos = "https://cloud.r-project.org")
}
cat("glmmTMB available:", requireNamespace("glmmTMB", quietly = TRUE), "\n")
