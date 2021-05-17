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
    wbd <- select(wbd, HUC12, TOHUC)
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