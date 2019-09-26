par_fun <- function(start_comid, net_atts, net_prep, wbd_atts, temp_dir) {
  library(nhdplusTools)
  library(hyRefactor)
  library(sf)
  library(dplyr)
  
  out_file <- paste0(temp_dir, start_comid, ".rds")
  
  if(!file.exists(out_file)) {
    
    message(paste(Sys.time(), out_file, "\n"))
    
    all_comid <- get_UT(net_atts, start_comid)
    
    sub_net <- filter(net_prep, COMID %in% all_comid)
    
    out <- list(NULL)
    
    if(nrow(sub_net) > 0) {
      sub_net <- bind_rows(sub_net,
                           select(filter(wbd_atts,
                                         !HUC12 %in% sub_net$HUC12),
                                  HUC12, TOHUC))
      
      out <- list(match_levelpaths(sub_net, start_comid, add_checks = TRUE))
    }
    
    saveRDS(out, out_file)
  }
  return(out_file)
}

get_hu_joiner <- function(net, wbd, simp, cores, temp_dir = "temp/", out_dir = "") {
  
  out_file <- file.path(out_dir, "map_joiner.csv")
  
  if(file.exists(out_file)) {
    all <- readr::read_csv(out_file)
  } else {
  
  unlink(temp_dir, recursive = TRUE)
  
  net_int <- get_process_data(net, wbd, simp)
  
  net <- st_set_geometry(net, NULL)[ ,1:40]
  wbd <- st_set_geometry(wbd, NULL)
  
  terminals <- net %>%
    select(TerminalPa) %>%
    distinct() %>%
    left_join(select(net, COMID, LevelPathI, Hydroseq), 
              by = c("TerminalPa" = "LevelPathI")) %>%
    filter(COMID %in% net_int$COMID) %>%
    group_by(TerminalPa) %>%
    filter(Hydroseq == min(Hydroseq))
  
  cl <- parallel::makeCluster(rep("localhost", cores), 
                              type = "SOCK", 
                              outfile = "hu_joiner.log")
  
  to_run <- terminals$COMID
  
  dir.create(temp_dir, showWarnings = FALSE)
  already_run <- list.files(temp_dir, pattern = "*.rds")
  already_run <- as.numeric(gsub(".rds", "", already_run))
  
  to_run <- to_run[!to_run %in% already_run]
  
  all_outlets <- parLapply(cl, to_run, par_fun,
                           net_atts = net,
                           net_prep = net_int,
                           wbd_atts = wbd,
                           temp_dir = temp_dir)
  
  stopCluster(cl)
  
  out_files <- list.files(temp_dir, full.names = TRUE)
  
  all <- sapply(out_files, readRDS, USE.NAMES = TRUE)
  
  names(all) <- gsub(temp_dir, "", names(all))
  names(all) <- gsub(".rds", "", names(all))
  
  all <- bind_rows(all)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  readr::write_csv(all, out_file)
  
  }
  
  return(all)
}

get_hu_outlets <- function(hu12, linked_points, add_to_file = NA) {
  hu12 <- st_set_geometry(hu12, NULL)
  
  hu12_sort <- names(igraph::topo_sort(igraph::graph_from_data_frame(
    hu12, directed = TRUE), mode = "out"))
  
  
  hu12_sort <- hu12_sort[hu12_sort %in% hu12$HUC12]
  
  
  hu12 <- left_join(tibble(HUC12 = hu12_sort),
                    hu12, by = "HUC12")
  
  hu12[["sort"]] <- seq_len(nrow(hu12))
  
  hu10 <- get_hu_outlet(hu12, 10, "HUC10")
  hu08 <- get_hu_outlet(hu12, 8, "HUC8")
  hu06 <- get_hu_outlet(hu12, 6, "HUC6")
  hu04 <- get_hu_outlet(hu12, 4, "HUC4")
  hu02 <- get_hu_outlet(hu12, 2, "HUC2")
  
  linked_points <- linked_points %>%
    left_join(hu02, by = c("HUC12" = "outlet_HUC12")) %>%
    left_join(hu04, by = c("HUC12" = "outlet_HUC12")) %>%
    left_join(hu06, by = c("HUC12" = "outlet_HUC12")) %>%
    left_join(hu08, by = c("HUC12" = "outlet_HUC12")) %>%
    left_join(hu10, by = c("HUC12" = "outlet_HUC12"))
  
  if(!is.na(add_to_file)) {
    write_sf(linked_points, add_to_file, "hu_outlets")
  }
  
  linked_points
}

get_hu_outlet <- function(hu12, hu_size, hu_name) {
  hu_outlet <- hu12 %>%
    mutate(hu = substr(.$HUC12, 1, hu_size)) %>%
    group_by(hu) %>%
    filter(sort == max(sort)) %>%
    ungroup() %>%
    mutate(tohu = ifelse(grepl("^[0-9].*", TOHUC), substr(TOHUC, 1, hu_size), TOHUC)) %>%
    select(hu, tohu, outlet_HUC12 = HUC12)
  
  names(hu_outlet) <- c(hu_name, paste0("to", hu_name), "outlet_HUC12")
  
  hu_outlet
}

par_hr_pairs <- function(x, prj, nhdplus_hw_outlets) {
  cats <- sf::read_sf(x, "NHDPlusCatchment")
  cats <- sf::st_transform(cats, prj)
  cats <- nhdplusTools:::rename_nhdplus(cats)
  
  dplyr::filter(mainstems::get_hr_pair(nhdplus_hw_outlets, cats), 
                !is.na(FEATUREID))
}

#' @export
get_hr_pairs <- function(nhdhr_path, nhdplus_hw_outlets, prj, cores) {
  gdb_files <- list.files(nhdhr_path, pattern = ".*.gdb$",
               full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
  
  cl <- parallel::makeCluster(rep("localhost", cores), type = "SOCK")
  hr_pairs <- parallel::parLapply(cl, gdb_files, par_hr_pairs, prj = prj,
                                  nhdplus_hw_outlets = nhdplus_hw_outlets)
  dplyr::bind_rows(hr_pairs)
}

get_hw_pairs <- function(hw_points, nhdp_cats) {
  sf::st_join(hw_points,
              dplyr::select(nhdp_cats, FEATUREID),
              join = sf::st_within) %>%
    st_set_geometry(NULL)
}
