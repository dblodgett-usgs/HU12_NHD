nhdhr_mod <- function(nhdhr_path, out_gpkg, min_size, simp, proj, force_terminal, fix_terminals = FALSE) {
  
  if(!file.exists(out_gpkg)){
    gdbs <- do.call(c, as.list(sapply(nhdhr_path, list.files, 
                              pattern = ".*[0-9][0-9][0-9][0-9].*.gdb$",
                              full.names = TRUE, recursive = TRUE, 
                              include.dirs = TRUE, USE.NAMES = FALSE)))
    
    get_hr_data_fun <- function(hr_gdb, min_size, simp, proj) {
      print(hr_gdb)
      
      hr_data <- nhdplusTools:::get_hr_data(hr_gdb, "NHDFlowline")
      
      hr_data <- nhdplusTools:::rename_nhdplus(hr_data)
      
      hr_data <- st_zm(hr_data)
      hr_data <- st_transform(hr_data, proj)
      hr_data <- st_simplify(hr_data, dTolerance = simp)
      
      filter_data <- select(st_set_geometry(hr_data, NULL), LevelPathI, TotDASqKM) %>%
        group_by(LevelPathI) %>%
        filter(TotDASqKM == max(TotDASqKM)) %>%
        ungroup() %>%
        distinct() %>%
        filter(TotDASqKM > min_size)
      
      hr_data <- hr_data[hr_data$LevelPathI %in% filter_data$LevelPathI, ]
      
      hr_data <- select(hr_data, -Permanent_Identifier, -FDate, -Resolution, GNIS_ID, 
                        -GNIS_Name, LENGTHKM, REACHCODE, -FlowDir, -WBArea_Permanent_Identifier,
                        FTYPE, -FCode, -MainPath, -InNetwork, VisibilityFilter, -Shape_Length, 
                        COMID, VPUID, -Enabled, -Shape, StreamLeve, StreamOrde, StreamCalc, 
                        FromNode, ToNode, Hydroseq, LevelPathI, Pathlength, TerminalPa, ArbolateSu, 
                        Divergence, StartFlag, TerminalFl, UpLevelPat, UpHydroSeq, DnLevel, 
                        DnLevelPat, DnHydroseq, DnMinorHyd, -DnDrainCou, FromMeas, ToMeas, 
                        -RtnDiv, -Thinner, -VPUIn, -VPUOut, AreaSqKM, TotDASqKM, -DivDASqKm, 
                        -MaxElevRaw, -MinElevRaw, -MaxElevSmo, -MinElevSmo, -Slope, -SlopeLenKm, 
                        -ElevFixed, HWType,-HWNodeSqKm, -StatusFlag)
    }
    
    hr_data <- pblapply(gdbs, get_hr_data_fun, min_size = 6, simp = simp, proj = proj)
    hr_data <- do.call(rbind, hr_data)
    hr_data <- st_sf(hr_data)
    
    terminals <- hr_data[hr_data$TerminalFl == 1, ]
    
    terminal_test <- hr_data$TerminalPa %in% terminals$TerminalPa
    
    if(fix_terminals) {
      hr_data_null <- hr_data[!terminal_test, ]
      
      outlets <- st_set_geometry(hr_data_null, NULL) %>%
        group_by(TerminalPa) %>%
        filter(Hydroseq == min(Hydroseq)) %>%
        select(Hydroseq, TerminalPa)
      
      for(term in unique(hr_data_null$TerminalPa)) {
        hr_data$TerminalFl[hr_data$Hydroseq == term] <- 1
        hr_data$TerminalPa[hr_data$TerminalPa == term] <- outlets$Hydroseq[outlets$TerminalPa == term]
      }
      
    } else {
      
      warning(paste("Removing", sum(!terminal_test), 
                    "flowlines that are missing terminal paths."))
      
      hr_data <- hr_data[terminal_test, ]
    }
    
    hr_data <- hr_data[(hr_data$FTYPE != 566 & hr_data$TerminalFl != 1), ]
    
    if(force_terminal) {
      t_atts <- select(st_set_geometry(hr_data, NULL), COMID, ToNode, FromNode, TerminalFl)
      
      t_atts <- left_join(t_atts, select(t_atts,
                                         toCOMID = COMID,
                                         FromNode),
                          by = c("ToNode" = "FromNode"))
      
      na_t_atts <- filter(t_atts, is.na(t_atts$toCOMID) & TerminalFl == 0) 
      
      warning(paste("Found", nrow(na_t_atts), "broken outlets where no toNode and not terminal. Fixing."))
      
      hr_data$TerminalFl[which(hr_data$COMID %in% na_t_atts$COMID)] <- 1
      
      # out <- filter(hr_data, COMID %in% na_t_atts$COMID)
      # 
      # write_sf(out, "./bad.gpkg")
    }
    
    write_sf(hr_data, layer = "NHDFlowline", dsn = out_gpkg)
  } else {
    hr_data <- read_sf(out_gpkg, layer = "NHDFlowline")
  }
  return(hr_data)
}

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
  hr_fline <- left_join(st_set_geometry(hr_fline, NULL), 
                        select(st_set_geometry(hr_fline, NULL),
                               toCOMID = .data$NHDPlusID,
                               .data$FromNode),
                        by = c("ToNode" = "FromNode"))
  
  hr_fline$TerminalFl[which(is.na(hr_fline$toCOMID))] <- 1
  
  hr_fline <- prepare_nhdplus(hr_fline, 
                              min_network_size = 0, 
                              min_path_length = 0, 
                              min_path_size = 0, 
                              purge_non_dendritic = TRUE) %>%
    select(ID = COMID, toID = toCOMID, length = LENGTHKM) %>%
    left_join(select(hr_fline,
                     ID = NHDPlusID, area = AreaSqKm, nameID = GNIS_ID),
              by = "ID") %>%
    distinct()
  
  hr_fline["weight"] <- nhdplusTools::calculate_arbolate_sum(select(hr_fline, ID, toID, length))
  
  hr_fline <- left_join(hr_fline, 
                        nhdplusTools::calculate_levelpaths(hr_fline), 
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


get_net <- function(net, prj) {
  
  if("NHDPlusID" %in% names(net)) {
    net <- rename(net, COMID = NHDPlusID, LENGTHKM = LengthKM, FTYPE = FType, 
                  TotDASqKM = TotDASqKm, Hydroseq = HydroSeq, Pathlength = PathLength,
                  AreaSqKM = AreaSqKm, DnHydroseq = DnHydroSeq)
    net$TerminalFl[which(!net$ToNode %in% net$FromNode)] <- 1
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

