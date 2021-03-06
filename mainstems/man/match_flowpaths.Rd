% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match_flowpath.R
\name{match_flowpaths}
\alias{match_flowpaths}
\title{Match Flowpaths}
\usage{
match_flowpaths(source_flowline, target_flowline, hw_pair, cores = NA)
}
\arguments{
\item{source_flowline}{sf data.frame with COMID and LevelPathI attributes.}

\item{target_flowline}{sf data.frame either NHDPlusHR with NHDPlusID, LENGTHKM, DnHydroseq, 
Hydroseq, and LevelPathI or an ad-hoc network with ID and toID attributes.}

\item{hw_pair}{data.frame with headwater pairs}

\item{cores}{integer number of cores to use in NHDPlus downstream network navigation.}
}
\value{
data.frame with the headwater id from the source_flowline input and 
levelpath from the source_flowline input for each identifier matched from the 
target dataset.
}
\description{
Implements a flowpath-matching algorithm that traces downstream along
the target flowline network and determines which levelpath from the source flowlines
best matches the resulting downstream traces. The algorithm starts from the outlet
location of the upstream most catchment in the source flowlines to stay away from
complexity that occurs near drainage divides.
}
\examples{
source(system.file("extdata/nhdplushr_data.R", package = "mainstems"))
source(system.file("extdata/new_hope_data.R", package = "nhdplusTools"))

hr_catchment <- nhdplusTools::align_nhdplus_names(hr_catchment)
hw_pair <- sf::st_set_geometry(sf::st_join(get_hw_points(new_hope_flowline),
                              dplyr::select(hr_catchment, FEATUREID),
                              join = sf::st_within), NULL)

lp_df_df <- match_flowpaths(new_hope_flowline, hr_flowline, hw_pair)
matched <- dplyr::left_join(dplyr::select(hr_flowline, COMID),
                            dplyr::select(lp_df_df, member_NHDPlusID,
                                          MR_LevelPathI = mr_LevelPathI), 
                                          by = c("COMID" = "member_NHDPlusID"))

lp <- min(matched$MR_LevelPathI, na.rm = TRUE)
mr_lp <- dplyr::filter(new_hope_flowline, LevelPathI <= lp)
hr_lp <- dplyr::filter(matched, MR_LevelPathI <= lp)
plot(sf::st_geometry(matched), col = "blue", lwd = 0.5)
plot(sf::st_geometry(mr_lp), col = "red", lwd = 3, add = TRUE)
plot(sf::st_geometry(hr_lp), col = "black", add = TRUE)

}
