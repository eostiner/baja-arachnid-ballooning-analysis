# R packages

The retained workflow uses:

- spatial/data: `sf`, `terra`, `readxl`, `readr`, `dplyr`, `tidyr`, `purrr`, `tibble`, `jsonlite`;
- richness/turnover: `iNEXT`, `vegan`;
- models: `mgcv`, `glmmTMB`;
- figures: `ggplot2`, `patchwork`, `maps`, `scales`, `viridis`, `svglite`.

The optional exploratory Step 12J additionally requires `gdm`.

The provided `environment.yml` requests these packages from conda-forge. After the final rerun, save `sessionInfo()` with the archived analysis outputs and use it when preparing the immutable GitHub/Zenodo release.

Step 12N additionally requires `ade4` for Outlying Mean Index environmental niche analysis.

The optional exploratory Step 12N OMI environmental-niche analysis additionally requires `ade4`.
