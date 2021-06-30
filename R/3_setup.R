# For splitting nhdplusV2 in the Mississippi
split_terminals <- data.frame(outlet = c("ohio", "uppermiss", "missouri", "middlemiss", "lowermiss"), 
                              COMID = c(1844789, 880678, 6018266, 5093446, 22811611),
                              stop = c(0, 0, 0, 880678, 5093446))

rename_hr_fl = function(hr_vpus) {
  new_names <- gsub("NHDPLUS_H_", "", gsub("_HU4_GDB", "", hr_vpus))
  
  if(!all(new_names == hr_vpus)) {
    for(i in seq_along(hr_vpus)) {
      file.rename(hr_vpus[i], new_names[i])
    }
  }
  
  return(new_names)
}

prep_nhdplushr <- function(hr_fline) {
  hr_fline <- filter(hr_fline, !is.na(.data$ToNode))
  
  hr_fline <- left_join(st_set_geometry(hr_fline, NULL), 
                        select(st_set_geometry(hr_fline, NULL),
                               toCOMID = .data$COMID,
                               .data$FromNode),
                        by = c("ToNode" = "FromNode"))
  
  hr_fline$TerminalFl[which(is.na(hr_fline$toCOMID))] <- 1
  
  hr_fline <- prepare_nhdplus(hr_fline, 
                              min_network_size = 0, 
                              min_path_length = 0, 
                              min_path_size = 0, 
                              purge_non_dendritic = FALSE) %>%
    select(ID = COMID, toID = toCOMID, length = LENGTHKM) %>%
    left_join(select(hr_fline,
                     ID = COMID, area = AreaSqKM, nameID = GNIS_ID),
              by = "ID") %>%
    distinct()
  
  hr_fline["weight"] <- nhdplusTools::calculate_arbolate_sum(select(hr_fline, ID, toID, length))
  
  hr_fline <- left_join(hr_fline, 
                        nhdplusTools::get_levelpaths(hr_fline), 
                        by = "ID")
  
  hr_fline <- left_join(hr_fline, 
                        select(hr_fline, ID, 
                               down_topo_sort = topo_sort,
                               down_levelpath = levelpath), 
                        by = c("toID" = "ID"))
  
  
  select(hr_fline, NHDPlusID = ID, 
         LevelPathI = levelpath, 
         DnLevelPat = down_levelpath, 
         DnHydroSeq = down_topo_sort, 
         HydroSeq = topo_sort)
}


get_net <- function(net, prj, nhdplus_update = NULL) {
  
  if("NHDPlusID" %in% names(net)) {
    net <- rename(net, COMID = NHDPlusID, LENGTHKM = LengthKM, FTYPE = FType, 
                  TotDASqKM = TotDASqKm, Hydroseq = HydroSeq, Pathlength = PathLength,
                  AreaSqKM = AreaSqKm, DnHydroseq = DnHydroSeq)
    net$TerminalFl[which(!net$ToNode %in% net$FromNode)] <- 1
  } else {
    new_atts <- data.table::fread(nhdplus_update, 
                                  integer64 = "character")
    
    net <- select(net, COMID) %>%
      right_join(new_atts, by = c("COMID" = "comid")) %>%
      align_nhdplus_names()
    
    # Required in the dataset but not used.
    net$StreamOrde <- 1
    net$StreamCalc <- 1
  }
  return(net %>%
           st_zm() %>%
           st_transform(prj)
  )
}

get_wbd <- function(wbd_gdb, fixes, prj) {
  wbd <- tryCatch(read_sf(wbd_gdb, "HUC12"), 
                  error = function(e) read_sf(wbd_gdb, "WBDHU12"))
  
  if("HUC_12" %in% names(wbd)) {
    wbd <- select(wbd, HUC12 = HUC_12, TOHUC = HU_12_DS)
  } else {
    wbd <- select(wbd, HUC12 = huc12, TOHUC = tohuc)
  }
  
  wbd <- filter(wbd, !grepl("^20.*|^19.*|^21.*|^22.*", wbd$HUC12))
  
  if(!is.null(fixes)) {
    # # Check if what we have is a DAG.
    for(fix in 1:nrow(fixes)) {
      fix_huc <- fixes$HUC12[fix]
      fix_tohuc <- fixes$TOHUC[fix]
      wbd_tohuc <- wbd$TOHUC[wbd$HUC12 == fix_huc]
      if(fix_tohuc != wbd_tohuc) {
        wbd$TOHUC[wbd$HUC12 == fix_huc] <- fix_tohuc
      } else {
        print(paste("unneeded fix", fix_huc))
      }
    }
  }
  
  if(any(wbd$HUC12 == wbd$TOHUC)) {
    error(paste("Some HUCs go to themselves.",  paste(wbd$HUC12[wbd$HUC12 == wbd$TOHUC], collapse = ", ")))
  }
  
  if(!igraph::is.dag(igraph::graph_from_data_frame(st_set_geometry(wbd, NULL)))) {
    g <- igraph::graph_from_data_frame(st_set_geometry(wbd, NULL))
    warning(paste("need to make HUCs into a DAG. Running accumulation to find them."))
    
    library(HUCAgg)
    fromHUC <-sapply(wbd$HUC12, fromHUC_finder,
                     hucs=wbd$HUC12, tohucs=wbd$TOHUC)
    
    aggrHUC <- sapply(wbd$HUC12, HUC_aggregator, fromHUC = fromHUC)
    message("investigate those and try again.")
  }
  
  wbd <- st_transform(wbd, prj)
  
  return(wbd)
}

clean_rf1 <- function(rf1) {
  from_to <- st_set_geometry(rf1, NULL) %>%
    filter(!TYPE %in% c("C", "G", "I", "L", "W", "X", "Z", "N")) %>%
    select(ID = ERF1_2., fnode_temp = FNODE_, tnode= TNODE_, name = PNAME, div_fraction = FRAC) %>%
    group_by(fnode_temp) %>%
    mutate(fnode = ifelse(div_fraction == max(div_fraction), fnode_temp, NA)) %>%
    ungroup() %>%
    select(-fnode_temp, -div_fraction)
  
  left_join(from_to, select(from_to, toID = ID, fnode), by = c("tnode" = "fnode")) %>%
    select(-fnode, -tnode) %>%
    left_join(select(rf1, ID = ERF1_2.), by = "ID") %>%
    st_sf() 
}

map_nhdpv1 <- function(ms, nhdpv1, nhdpv2, xwalk) {
  find <- select(ms, LevelPathI, outlet_nhdpv2_COMID, head_nhdpv2_COMID) %>%
    left_join(select(xwalk, V2_ComID, outlet_V1_ComID = V1_ComID, outlet_XWalkType = XWalkType), 
              by = c("outlet_nhdpv2_COMID" = "V2_ComID")) %>%
    left_join(select(xwalk, V2_ComID, head_V1_ComID = V1_ComID, head_XWalkType = XWalkType),
              by = c("head_nhdpv2_COMID" = "V2_ComID"))
  
  missing_outlet <- find %>%
    filter(is.na(outlet_V1_ComID))
  
  missing_head <- find %>%
    filter(is.na(head_V1_ComID))
  
  find <- find %>%
    filter(!is.na(outlet_V1_ComID) & !is.na(head_V1_ComID))
  
  find <- find %>%
    left_join(select(nhdpv1, COMID, HYDROSEQ), by = c("outlet_V1_ComID" = "COMID")) %>%
    group_by(LevelPathI) %>%
    filter(HYDROSEQ == min(HYDROSEQ) & HYDROSEQ != 0) %>%
    select(-HYDROSEQ) %>%
    ungroup() %>%
    left_join(select(nhdpv1, COMID, HYDROSEQ), by = c("head_V1_ComID" = "COMID")) %>%
    group_by(LevelPathI) %>%
    filter(HYDROSEQ == max(HYDROSEQ) & HYDROSEQ != 0) %>%
    select(-HYDROSEQ) %>%
    ungroup()
  
  dup <- group_by(find, LevelPathI) %>%
    filter(n() > 1)
  
  missing_head <- find$LevelPathI[is.na(find$head_V1_ComID)]
  missing_outlet <- find$LevelPathI[is.na(find$outlet_V1_ComID)]
  
  return(list(mapped = find, missing_head = missing_head, missing_outlet = missing_outlet))
}

get_hu02 <- function(wbd_path, prj, national_viz_simp) {
  wbd <- st_transform(read_sf(wbd_path, "WBDHU2"), prj)
  
  valid <- st_is_valid(st_geometry(wbd), NA_on_exception = TRUE)
  valid[is.na(valid)] <- FALSE
  
  wbd <- wbd[valid,]
  
  st_simplify(wbd, dTolerance = national_viz_simp)
}

get_merit <- function(mdir) {
  merit_files <- list.files(mdir, pattern = "riv_.*shp", full.names = TRUE)
  
  merit_cats <- list.files(mdir, pattern = "cat_.*shp", full.names = TRUE)
  
  riv <- lapply(merit_files, sf::read_sf)
  
  cat <- lapply(merit_cats, sf::read_sf)
  
  cat <- lapply(cat, sf::st_drop_geometry)
  
  cat <- bind_rows(cat)
  
  riv <- do.call(rbind, riv)
  
  riv <- left_join(riv, cat, by = "COMID")
  
  mutate(riv, nameID = "constant") %>%
    select(comid = COMID, 
           tocomid = NextDownID, 
           nameID, 
           lengthkm, 
           areasqkm = unitarea)
}

get_merit_cats <- function(mdir) {
  merit_cats <- list.files(mdir, pattern = "cat_.*shp", full.names = TRUE)
  do.call(rbind, lapply(merit_cats, sf::read_sf))
}

get_merit_atts <- function(riv, cores, merit_cache = "") {
  
  if(file.exists(merit_cache)) {
    message("returning cached merit atts")
    
    return(st_drop_geometry(read_sf(merit_cache, "merit_plus")))
  }
  
  riv_atts <- st_drop_geometry(riv)
  
  add_plus_network_attributes(riv_atts, 
                              cores = cores, 
                              status = TRUE)
  
}

write_merit <- function(riv, riv_atts, out_gpkg, simp = FALSE) {
  
  if(file.exists(out_gpkg)) {
    message("returning cached merit plus data")
    
    return(sf::read_sf(out_gpkg, "merit_plus"))
  }
  
  riv <- left_join(select(riv, comid), riv_atts, by = "comid")
  
  write_sf(riv, out_gpkg, "merit_plus")
  
  riv
}

write_merit_share <- function(merit, merit_w_names, out_gpkg) {
  merit <- select(merit, comid) %>%
    left_join(merit_w_names)
  
  sf::sf_use_s2(FALSE)
  
  merit <- mutate(merit, 
                  comid = as.integer(comid), 
                  tocomid = as.integer(tocomid), 
                  terminalpa = as.integer(terminalpa), 
                  hydroseq = as.integer(hydroseq), 
                  levelpathi = as.integer(levelpathi),
                  dnlevelpat = as.integer(dnlevelpat), 
                  dnhydroseq = as.integer(dnhydroseq),
                  terminalfl = as.integer(terminalfl)) 
  
  merit <- select(merit, -weight)
  
  merit <- mutate(merit, down_levelpaths = gsub(",0$|^0$", "", down_levelpaths))
  
  merit <- st_simplify(merit, dTolerance = 0.005)
  
  write_sf(merit, out_gpkg, "merit_plus")
  
  add_merit_indexes(out_gpkg)
  
  merit
}

add_merit_indexes <- function(x) {
  library(DBI)
  
  db <- dbConnect(RSQLite::SQLite(), x)
  
  dbExecute(db, "CREATE UNIQUE INDEX idx_comid ON merit_plus (comid);")
  
  add_idx <- c("nameID", "terminalpa", "levelpathi")
  
  sapply(add_idx, function(x, db) {
    dbExecute(db, paste0("CREATE INDEX idx_", x, " ON merit_plus (", x, ");"))
  }, db = db)
  
  dbDisconnect(db)
}

trim_line <- function(x, dist = 0.2) {
  
  coords <- st_coordinates(x) %>%
    as.data.frame() %>%
    nhdplusTools:::add_len() %>%
    filter(len > dist & len < max(len) - dist) %>%
    as.matrix()
  
  if(nrow(coords) < 2) {
    
    coords <- st_coordinates(x)
    
    mid <- round((nrow(coords) / 2))
    
    coords <- coords[(mid - 1):(mid + 1), 1:2]
    
  } else {
    
    coords <- coords[, 1:2]
    
  }
  
  return(st_linestring(coords))
  
}

get_trimmed <- function(naturalearth_data) {
  
  naturalearth_data_trim <- do.call(
    st_sfc, lapply(naturalearth_data$geometry, trim_line)
  )
  
  naturalearth_data_trim <- st_sf(st_drop_geometry(naturalearth_data), 
                                  geom = naturalearth_data_trim)
  
  sf::st_crs(naturalearth_data_trim) <- sf::st_crs(naturalearth_data)
  
  naturalearth_data_trim
  
}

get_ends <- function(r, nd1, nd2, naturalearth_data) {
  
  all_r_b <- naturalearth_data$rivernum == r
  
  all_r <- naturalearth_data[all_r_b, ]
  
  all_nd <- bind_rows(nd1[all_r_b, ], nd2[all_r_b, ])
  
  all_dist <- st_distance(all_nd, by_element = FALSE) %>%
    unclass()
  
  all_dist[lower.tri(all_dist, diag = TRUE)] <- NA
  
  all_dist <- tibble::as_tibble(all_dist) %>%
    tibble::rowid_to_column() %>%
    tidyr::gather(colid, distance, tidyr::starts_with("V")) %>%
    arrange(desc(distance))
  
  furthest <- all_dist[1, ]
  
  rbind(all_nd[furthest$rowid, ], 
        all_nd[gsub("V", "", furthest$colid), ])
  
}

get_fp_ends <- function(naturalearth_data_trim) {
  
  nd1 <- get_node(naturalearth_data_trim, "start")
  
  nd2 <- get_node(naturalearth_data_trim, "end")
  
  
  ends <- lapply(unique(naturalearth_data_trim$rivernum), get_ends, 
                 naturalearth_data = naturalearth_data_trim, 
                 nd1 = nd1, nd2 = nd2)
  
  ends <- do.call(rbind, ends)
  
  rivers <- unique(naturalearth_data_trim$rivernum)
  
  rivers <- rivers[rep(1:length(rivers), each = 2)]
  
  st_sf(rivernum = rivers, geom = sf::st_geometry(ends))
}

join_naturalearth_ends <- function(rivers, merit) {
  
  sf::sf_use_s2(FALSE)
  
  merit_simp <- 
    select(merit, COMID = comid) %>%
    mutate(REACHCODE = COMID, ToMeas = 100, FromMeas = 0)
  
  indexes <- merit_simp %>%
    nhdplusTools::get_flowline_index(rivers, search_radius = 0.05, 
                                     max_matches = 10)
  
  indexes %>%
    left_join(select(st_drop_geometry(merit), comid, totdasqkm), 
              by = c("COMID" = "comid")) %>%
    left_join(mutate(rivers, id = seq(1:nrow(rivers))), by = "id") %>%
    group_by(id) %>%
    filter(totdasqkm == max(totdasqkm))
  
}

get_naturalearth_heads <- function(naturalearth_data, merit) {
  naturalearth_data <- 
    dplyr::left_join(naturalearth_data, 
                     select(sf::st_drop_geometry(merit), 
                            comid, hydroseq, levelpathi, 
                            terminalpa, terminalfl), 
                     by = c("COMID" = "comid"))
  
  # Drop very small networks
  small_nets <- sf::st_drop_geometry(merit) %>%
    select(totdasqkm, terminalpa, hydroseq) %>%
    filter(terminalpa == hydroseq & totdasqkm < 1000)
  
  out <- naturalearth_data %>%
    filter(!terminalpa %in% small_nets$terminalpa &
             terminalfl != 1) %>%
    group_by(rivernum) %>%
    filter(!is.na(hydroseq) & hydroseq == max(hydroseq)) %>%
    ungroup()
  
}

