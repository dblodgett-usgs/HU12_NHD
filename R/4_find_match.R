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
  
  linked_points <- st_transform(linked_points, 4326)
  
  if(!is.na(add_to_file)) {
    write_sf(linked_points, add_to_file, "hu_points")
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

filter_dominant_length <- function(hu_joiner, lp_lengths, net, 
                                   factor_corrected = 2,
                                   factor_intersected = 10,
                                   remove_length_m = 200) {
  lp_lengths <- left_join(lp_lengths, select(net, COMID, LevelPathI), by = "COMID")
  
  lp_lengths$length <- as.numeric(lp_lengths$length)
  
  lp_lengths <- lp_lengths %>%
    select(LevelPathI, HUC12, length) %>%
    group_by(LevelPathI, HUC12) %>%
    summarize(sum_length = sum(length))
  
  dominant <- lp_lengths %>%
    group_by(HUC12) %>%
    filter(sum_length == max(sum_length)) %>%
    rename(dominant_LevelPathI = LevelPathI) %>%
    select(-sum_length) %>% ungroup()
  
  hu_joiner2 <- hu_joiner %>%
    rename(LevelPathI = corrected_LevelPathI) %>%
    left_join(dominant, by = "HUC12") %>%
    left_join(rename(ungroup(lp_lengths), lp_length = sum_length), 
              by = c("HUC12", "LevelPathI")) %>%
    left_join(rename(ungroup(lp_lengths), 
                     dominant_LevelPathI = LevelPathI, 
                     dlp_length = sum_length), 
              by = c("HUC12", "dominant_LevelPathI"))
  
  # Now we have length per corrected levelpath and dominant levelpath.
  # Can use these as metrics.
  
  switched <- filter(hu_joiner2, (dlp_length > factor_corrected * lp_length & LevelPathI != intersected_LevelPathI) | 
                       (dlp_length > factor_intersected * lp_length & LevelPathI == intersected_LevelPathI)) %>%
    filter() %>%
    select(-LevelPathI, corrected_LevelPathI = dominant_LevelPathI, -lp_length, -dlp_length) %>%
    mutate(dominant_override = TRUE)
  
  remove <- filter(hu_joiner2, dlp_length < remove_length_m & lp_length < remove_length_m)
  
  filter(hu_joiner, !HUC12 %in% switched$HUC12) %>%
    mutate(dominant_override = FALSE) %>%
    bind_rows(switched) %>%
    filter(!HUC12 %in% remove$HUC12)
}

add_name_paths <- function(net, heads) {
  
  g <- hy_make_graph(net)
  
  for(head in names(heads)) {
    
    path <- hy_dfs(heads[[head]], "down", g = g, data = FALSE)
    
    path <- path[path %in% net$comid]
    
    replace <- net$nameID == "unknown" & net$comid %in% path
    
    net$nameID[replace] <- rep(head, sum(replace))
    
  }
  
  return(net)
  
}

hy_make_graph <- function(net) {
  suppressWarnings(
    dplyr::select(net, comid, tocomid) %>%
      igraph::graph_from_data_frame(directed = TRUE))
}

hy_dfs <- function(start, mode = "up", net = NULL, g = NULL, data = TRUE) {
  
  map_mode <- list(up = "in", down = "out")
  
  if(is.null(g)) {
    
    g <- hy_make_graph(net)
    
  }
  
  start_node <- which(names(V(g)) == as.character(start))
  
  sub_net <- igraph::dfs(g, 
                         root = start_node, 
                         neimode = map_mode[[mode]], 
                         unreachable = FALSE)$order
  
  sub_net <- as.integer(names(sub_net[!is.na(sub_net)]))
  
  if(data) {
    return(dplyr::filter(net, comid %in% sub_net))
  } else {
    return(sub_net)
  }
}

add_names <- function(merit_natearth_heads, names, merit_atts) {
  
  names <- filter(names, !is.na(name)) %>%
    distinct() %>%
    group_by(rivernum) %>%
    filter(row_number() == 1) # this is a hack
  
  merit_natearth_heads <- arrange(merit_natearth_heads, desc(hydroseq))

  merit_natearth_heads <- left_join(merit_natearth_heads, 
                                    names, by = "rivernum")
    
  merit_natearth_heads <- filter(merit_natearth_heads, !is.na(name))
  
  heads <- setNames(as.list(merit_natearth_heads$COMID), 
                    merit_natearth_heads$name)
  
  merit_atts$nameID <- rep("unknown", nrow(merit_atts))
  
  add_name_paths(merit_atts, heads)
  
}
