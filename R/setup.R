cran_packages <- c("sf", "dplyr", "drake", "snow", "xml2", 
                   "readr", "R.utils", "tidyr", # "rmapshaper"
                   "readr", "pbapply", "remotes", "future")

usgs_r_github_packages <- c("nhdplusTools", "HUCAgg")

local_package <-  "mainstems"

igraph <- "igraph"

packages <- c(cran_packages, usgs_r_github_packages, local_package, igraph)

inst <- installed.packages()

if(!"devtools" %in% inst[, 1])
  install.packages("devtools", repos = "https://cloud.r-project.org")

if(any(!(ur_check <- usgs_r_github_packages %in% inst[, 1]))) {
  lapply(usgs_r_github_packages[!ur_check], 
         function(x) remotes::install_github(paste0("usgs-r/", x)))
}

if(!"units" %in% inst[, 1]) 
  install.packages('udunits2', configure.args='--with-udunits2-lib=/cxfs/projects/root/opt/udunits/2.2.16/gnu/lib --with-udunits2-include=/cxfs/projects/root/opt/udunits/2.2.16/gnu/include')

if(!igraph %in% inst[, 1]) {
  devtools::install_version('igraph', '1.2.2')
}
  
if(any(!(inst_check <- cran_packages %in% inst[, 1]))) {
  install.packages(cran_packages[!inst_check])
}

if(!local_package %in% inst[, 1]) {
  install.packages(local_package, repos = NULL)
}

lapply(packages, require, character.only = TRUE)

