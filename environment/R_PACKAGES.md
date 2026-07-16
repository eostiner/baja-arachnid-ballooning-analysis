# R packages

The retained workflow uses:

- spatial/data: `sf`, `terra`, `readxl`, `readr`, `dplyr`, `tidyr`, `purrr`, `tibble`, `jsonlite`;
- richness/turnover: `iNEXT`, `vegan`;
- models: `mgcv`, `glmmTMB`;
- figures: `ggplot2`, `patchwork`, `maps`, `scales`, `viridis`, `svglite`.

The optional exploratory Step 12J additionally requires `gdm`.

The provided `environment.yml` requests these packages from conda-forge. After the final rerun, save `sessionInfo()` with the archived analysis outputs and use it when preparing the immutable GitHub/Zenodo release.
