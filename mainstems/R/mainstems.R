.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste(strwrap(
    "USGS Support Package:
    https://owi.usgs.gov/R/packages.html#support"),
    collapse = "\n"))
  options("pboptions" = list(
    type = if (interactive()) "timer" else "none",
    char = "-",
    txt.width = 50,
    gui.width = 300,
    style = 3,
    initial = 0,
    title = "R progress bar",
    label = "",
    nout = 100L,
    min_time = 0))
  invisible(NULL)
}