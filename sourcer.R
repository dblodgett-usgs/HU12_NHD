cran_packages <- c("sf", "dplyr", "drake", "snow", "xml2", 
                   "readr", "rmapshaper", "R.utils", "tidyr",
                   "igraph", "readr", "pbapply", "remotes")

usgs_r_github_packages <- c("nhdplusTools", "HUCAgg")

local_package <-  "mainstems"

packages <- c(cran_packages, usgs_r_github_packages, local_package)

inst <- installed.packages()

if(any(!(inst_check <- cran_packages %in% inst[, 1]))) {
  install.packages(cran_packages[inst_check])
}

if(any(!(ur_check <- usgs_r_github_packages %in% inst[, 1]))) {
  lapply(usgs_r_github_packages[ur_check], 
         function(x) remotes::install_github(paste0("usgs-r/", x)))
}


if(!local_package %in% inst[, 1]) {
  install.packages(local_package)
}

lapply(packages, require, character.only = TRUE)

source("R/2_fixes.R")
source("R/3_setup.R")
source("R/4_find_match.R")
source("R/5_find_outlets.R")
source("R/6_visualize.R")
source("R/10_build_mainstems_table.R")
