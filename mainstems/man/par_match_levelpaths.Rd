% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match_levelpaths.R
\name{par_match_levelpaths}
\alias{par_match_levelpaths}
\title{Parallel match levelpaths}
\usage{
par_match_levelpaths(
  net,
  wbd,
  simp,
  cores,
  temp_dir = "temp/",
  out_file = "temp.csv",
  net_int = NULL,
  purge_temp = TRUE
)
}
\arguments{
\item{net}{NHDPlus network}

\item{wbd}{WBD HU12 polygons}

\item{simp}{simplification in units of WBD and NHDPlus to limit memory usage}

\item{cores}{Number of parallel processes to run. (1 for now process splitting)}

\item{temp_dir}{directory to write temprary files to}

\item{out_file}{file to write results to.}

\item{net_int}{intermediate process artifact to be used in drake plans output 
of get_process_data}

\item{purge_temp}{should temp directory be purged?}
}
\description{
Parallel match levelpaths
}
\details{
Calls match_levelpaths in a parallel mode for use in large process workflows
}
