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
  cats <- nhdplusTools::align_nhdplus_names(cats)
  
  pairs <- sf::st_set_geometry(
    sf::st_join(nhdplus_hw_outlets,
            dplyr::select(cats, FEATUREID),
            join = sf::st_within),
    NULL)
  
  dplyr::filter(pairs, 
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
