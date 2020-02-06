run_match <- function(comid) {
  f <- paste0("data/test_", comid, ".gpkg")
  
  net <- read_sf(f, "NHDFlowline_Network") %>%
    st_zm()
  hu <- read_sf(f, "HUC12") %>%
    rename(HUC12 = HUC_12, TOHUC = HU_12_DS)
  
  net_int <- mainstems:::get_process_data(net, hu, 0)
  
  matched <- match_levelpaths(net_int, comid)
  
  linked_points <- get_linked_points(matched, net, hu, exclude = c())
  
  unlink("temp.csv")
  
  return(linked_points)
}

comid <- 21412883

bad_net <- readRDS("tests/testthat/data/match_levelpaths_21412883.rds")

matched <- match_levelpaths(bad_net, comid)

comid <- 626220

matched <- par_match_levelpaths(net, hu, 0, 1, tempdir(check = TRUE))
