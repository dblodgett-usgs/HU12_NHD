build_mainstem_table <- function(nhdplus_net, 
                                 olp, 
                                 nlp,
                                 rfp,
                                 rf1,
                                 nhdp_wbd,
                                 wbd,
                                 nhd1) {
  
  nhdplus_net <- left_join(nhdplus_net, 
                           rename(nhdplusTools::get_vaa(atts = c("gnis_id", "gnis_name")),
                                  COMID = comid, GNIS_ID = gnis_id, GNIS_NAME = gnis_name),
                           by = "COMID")
  
# Can add this back from it's not in the new nhgdplus attributes.  
  GNIS <- select(nhdplus_net, LevelPathI, Hydroseq, GNIS_ID, GNIS_NAME) %>%
    filter(GNIS_ID != " ") %>%
    group_by(LevelPathI) %>%
    filter(Hydroseq == min(Hydroseq)) %>%
    ungroup() %>%
    select(LevelPathI, outlet_GNIS_ID = GNIS_ID, outlet_GNIS_NAME = GNIS_NAME)

  ms <- data.frame("LevelPathI" = unique(nhdplus_net$LevelPathI))
  
  ms <- ms %>%
    left_join(GNIS, by = "LevelPathI") %>%
    left_join(select(get_lp_outlets(nhdplus_net),
                     LevelPathI, 
                     outlet_nhdpv2_COMID = COMID), by = "LevelPathI") %>%
    left_join(select(get_lp_heads(nhdplus_net),
                     LevelPathI, 
                     head_nhdpv2_COMID = COMID), by = "LevelPathI") %>%
    left_join(rename(get_hulp(olp, nhdplus_net, nhdp_wbd), outlet_nhdpv2HUC12 = outlet_HUC12, head_nhdpv2HUC12 = head_HUC12), by = "LevelPathI") %>%
    left_join(rename(get_hulp(nlp, nhdplus_net, wbd), outlet_latestHUC12 = outlet_HUC12, head_latestHUC12 = head_HUC12), by = "LevelPathI") %>%
    left_join(get_rf1lp(rf1, rfp, nhdplus_net), by = "LevelPathI")
  
  
}

get_lp_outlets <- function(nhdplus_net) {
  nhdplus_net %>%
    group_by(.data$LevelPathI) %>%
    filter(.data$Hydroseq == min(.data$Hydroseq)) %>%
    ungroup()
}

get_lp_heads <- function(nhdplus_net) {
  nhdplus_net %>%
    group_by(.data$LevelPathI) %>%
    filter(.data$Hydroseq == max(.data$Hydroseq)) %>%
    ungroup()
}

get_hulp <- function(hulp, nhdplus_net, wbd) {
  hulp <- select(st_set_geometry(hulp, NULL), .data$COMID, .data$HUC12, .data$LevelPathI) %>%
    filter(!is.na(.data$HUC12)) %>%
    left_join(select(nhdplus_net, .data$COMID, .data$Hydroseq), by = "COMID")
  
  heads <- get_lp_heads(hulp)
  
  outlets <- get_lp_outlets(hulp)
  
  hu_sort <- names(igraph::topo_sort(igraph::graph_from_data_frame(
    wbd, directed = TRUE), mode = "out"))
  
  hu_sort <- data.frame(HUC12 = hu_sort, sort = seq_along(hu_sort), 
                        stringsAsFactors = FALSE)
  
  hu <- hulp %>%
    left_join(select(outlets, outlet_HUC12 = .data$HUC12, LevelPathI), by = "LevelPathI") %>%
    left_join(select(heads, head_HUC12 = .data$HUC12, LevelPathI), by = "LevelPathI") %>%
    select(LevelPathI, outlet_HUC12, head_HUC12) %>%
    distinct() %>%
    left_join(hu_sort, by = c("head_HUC12" = "HUC12")) %>%
    group_by(LevelPathI) %>%
    filter(sort == max(sort)) %>%
    select(-sort) %>%
    ungroup() %>%
    left_join(hu_sort, by = c("outlet_HUC12" = "HUC12")) %>%
    group_by(LevelPathI) %>%
    filter(sort == min(sort)) %>%
    ungroup() %>%
    select(-sort)

}

get_rf1lp <- function(rf1, rfp, nhdplus_net) {
  rfp <- rfp %>%
    select(-.data$headwater_COMID) %>%
    distinct()
  
  rf1_sort <- data.frame(ID = as.integer(
    names(topo_sort(graph_from_data_frame(select(st_set_geometry(rf1, NULL), 
                                                 .data$ID, .data$toID)), 
                    mode = "out"))))
  
  rf1_sort["Hydroseq"] <- seq_len(nrow(rf1_sort))
  
  rfp <- left_join(rfp, rf1_sort, by = c("member_ID" = "ID"))
  
  heads <- get_lp_heads(rename(rfp, LevelPathI = .data$mr_LevelPathI))
  
  outlets <- get_lp_outlets(rename(rfp, LevelPathI = .data$mr_LevelPathI))
  
  rename(rfp, LevelPathI = .data$mr_LevelPathI) %>%
    left_join(select(outlets, outlet_rf1ID = .data$member_ID, LevelPathI), by = "LevelPathI") %>%
    left_join(select(heads, head_rf1ID = .data$member_ID, LevelPathI), by = "LevelPathI") %>%
    select(LevelPathI, outlet_rf1ID, head_rf1ID) %>%
    distinct()
}

add_v1 <- function(mainstems_table, nhdpv1_mapped) {
  mainstems_table <- mainstems_table %>%
    left_join(select(nhdpv1_mapped$mapped, LevelPathI, 
                     outlet_nhdpv1_COMID = outlet_V1_ComID, 
                     head_nhdpv1_COMID = head_V1_ComID), by = "LevelPathI")
}


make_ms_summary <- function(ms, nhdp_att) {
  nhdp_att <- align_nhdplus_names(nhdp_att)
  
  if(!"StreamLeve" %in% names(nhdp_att)) {
    
    nhdp_att$StreamLeve <-select(nhdp_att, 
                                 levelpathi = LevelPathI, 
                                 dnlevelpat = DnLevelPat) %>%
      nhdplusTools::get_streamlevel()
    
  }
  
  area <- select(nhdp_att, ID = COMID, toID = toCOMID, area = AreaSqKM) %>%
    calculate_total_drainage_area()
  
  ms_length <- nhdp_att %>%
    select(LevelPathI, LENGTHKM) %>%
    group_by(LevelPathI) %>%
    summarize(length = sum(LENGTHKM))

  summary <- data.frame(COMID = nhdp_att$COMID, totdasqkm = area) %>%
    left_join(select(nhdp_att, COMID, LevelPathI, level = StreamLeve), 
              by = "COMID") %>%
    left_join(ms_length, by = "LevelPathI") %>%
    select(-LevelPathI)
  
  ms$outlet_GNIS_ID <- as.integer(ms$outlet_GNIS_ID)
  ms$outlet_GNIS_NAME[ms$outlet_GNIS_NAME == " "] <- NA
  
  left_join(ms, summary, by = c("outlet_nhdpv2_COMID" = "COMID"))
}

get_hist_list <- function(ms) {

  make_hist_df <- function(h) {
    data.frame(breaks = h$breaks[2:length(h$breaks)], 
               counts = h$counts)
  }
  
  breaks <- seq(0,3000000, length.out = 11)
  
  # size distribution of GNIS_ID available basins
  gnis <- hist(filter(ms, !is.na(outlet_GNIS_ID))$totdasqkm, breaks = breaks)
  
  # size distribution of NA GNIS_ID basins
  na_gnis <- hist(filter(ms, is.na(outlet_GNIS_ID))$totdasqkm, breaks = breaks)
  
  huc12 <- hist(filter(ms, !is.na(outlet_nhdpv2HUC12))$totdasqkm, breaks = breaks)
  na_huc12 <- hist(filter(ms, is.na(outlet_nhdpv2HUC12))$totdasqkm, breaks = breaks)
  
  rf1 <- hist(filter(ms, !is.na(outlet_rf1ID))$totdasqkm, breaks = breaks)
  na_rf1 <- hist(filter(ms, is.na(outlet_rf1ID))$totdasqkm, breaks = breaks)
  
  list(gnis_hist = make_hist_df(gnis), 
       na_gnis = make_hist_df(na_gnis),
       huc12 = make_hist_df(huc12),
       na_huc12 = make_hist_df(na_huc12),
       rf1 = make_hist_df(rf1),
       na_rf1 = make_hist_df(na_rf1))
}

ds_till <- function(atts, h) {
  
}

find_v1_mainstems <- function(v1, ms, cores = NA) {
  ms <- filter(ms, !is.na(outlet_latestHUC12) & !is.na(head_latestHUC12))
  
  heads <- unique(ms$head_nhdpv1_COMID)
  
  heads <- heads[!is.na(heads)]
  
  # Get the v1 down hydroseq
  v1_dnhs <- left_join(distinct(select(v1, COMID, TONODE)), 
                       distinct(select(v1, FROMNODE, DnHydroseq = HYDROSEQ)), 
                       by = c("TONODE" = "FROMNODE")) %>%
    select(-TONODE) %>%
    distinct()
  
  # Deduplicate by removing the down minor hydro path
  v1_dnhs <- filter(v1_dnhs, !DnHydroseq %in% v1$DNMINHYDRO) %>%
    distinct()
  
  # Further deduplicate by removing the divergence = 2 toCOMIDs.
  div_heads <- select(v1, COMID, DIVERGENCE) %>%
    filter(DIVERGENCE == 2)
  
  v1_dnhs <- left_join(v1_dnhs, select(v1, toCOMID = COMID, HYDROSEQ), 
                       by = c("DnHydroseq" = "HYDROSEQ")) %>%
    filter(!toCOMID %in% div_heads$COMID)
  
  v1 <- left_join(v1, select(v1_dnhs, -toCOMID), by = "COMID")
  
  v1 <- select(v1, COMID, LENGTHKM, DnHydroseq, Hydroseq = HYDROSEQ, 
               LevelPathI = LEVELPATHI, DnLevelPat = DNLEVELPAT)
  
  hs <- select(v1, Hydroseq = Hydroseq) %>%
    distinct() %>%
    arrange(Hydroseq)
  
  hs$hs <- c(1:nrow(hs))
  
  v1 <- left_join(v1, hs, by = "Hydroseq") %>%
    select(-Hydroseq) %>%
    rename(Hydroseq = hs) %>%
    left_join(hs, by = c("DnHydroseq" = "Hydroseq")) %>%
    select(-DnHydroseq) %>%
    rename(DnHydroseq = hs) %>%
    left_join(hs, by = c("LevelPathI" = "Hydroseq")) %>%
    select(-LevelPathI) %>%
    rename(LevelPathI = hs) %>%
    left_join(hs, by = c("DnLevelPat" = "Hydroseq")) %>%
    select(-DnLevelPat) %>%
    rename(DnLevelPat = hs)
  
  cl <- NULL
  if(!is.na(cores)) {
    cl <- parallel::makeCluster(rep("localhost", cores), type = "SOCK")
  }
  
  nets <- pblapply(heads, function(x, net) nhdplusTools::get_DM(net, x), net = v1, cl = cl)
  
  parallel::stopCluster(cl)
  
  return(nets)
}

get_dlp <- function(x, dnlp) {
  out <- dnlp[x]
  
  if(out == 0) {
    return(out)
  }
  
  c(out, get_dlp(out, dnlp))
  
}

add_dm <- function(merit) {
  lp <- select(merit, levelpathi, dnlevelpat) %>%
    filter(levelpathi != dnlevelpat)
  
  dnlp <- data.frame(lp = seq(0, (max(lp$levelpathi)))) %>%
    left_join(lp, by = c("lp" = "levelpathi"))
  
  dnlp <- dnlp$dnlevelpat[2:nrow(dnlp)]
  
  all_lp <- pbapply::pblapply(unique(lp$levelpathi), get_dlp, dnlp = dnlp)
  
  all_lp <- data.frame(levelpathi = unique(lp$levelpathi),
                       dnlp = sapply(all_lp, paste, collapse = ","))
  
  merit %>%
    left_join(all_lp, by = "levelpathi") %>%
    rename(down_levelpaths = dnlp)
}

upload_sb <- function(mainstems_table_summary, geo_summary) {
  
  upload_list <- list(`nhdplusv2wbd.csv` = "out/nhdplus_oldwbd/map_joiner.csv",
                      `newwbd.csv` = "out/nhdplus_newwbd/map_joiner.csv",
                      `rf1.csv` = "out/rf1_out/map_joiner.csv", 
                      `mainstems_summary.gpkg` = "out/mainstems_summary.gpkg")
  
  
  upload_sb_fun(upload_list)
  
}

upload_sb_fun <- function(upload_list, sb_id = "60cb5edfd34e86b938a373f4") {
  sbtools::authenticate_sb()
  
  files <- sbtools::item_list_files(sb_id)
  
  upload_list <- upload_list[!names(upload_list) %in% files$fname]
  
  for(f in names(upload_list)) {
    fi <- upload_list[[f]]
    
    upload_file <- file.path(dirname(fi), f)
    
    if(file.rename(fi, upload_file)) {
      message(upload_file)
      
      try(sbtools::item_append_files(sb_id = sb_id, files = upload_file))
      file.rename(upload_file, fi)
    }
  }
}

make_geo_summary <- function(nhdplus_net, mainstems_table_summary, out_file) {
  
  geom <- select(nhdplus_net, LevelPathI, Hydroseq) %>%
    arrange(desc(Hydroseq)) %>%
    group_by(LevelPathI) %>%
    group_split(.keep = TRUE)

  cl <- parallel::makeCluster(8)
  
  geoms <- pbapply::pblapply(geom, get_single_line, cl = cl)
  
  geoms <- do.call(c, geoms)
  
  geoms <- st_sf(LevelPathI = sort(unique(nhdplus_net$LevelPathI)),
                 geom = geoms)
  
  mainstems_table_summary <- left_join(mainstems_table_summary, 
                                       geoms, by = "LevelPathI")
  
  mainstems_table_summary <- sf::st_sf(mainstems_table_summary)
  
  write_sf(mainstems_table_summary, out_file)
}

get_single_line <- function(x) {
  
  sf::st_sfc(
    sf::st_simplify(
    sf::st_linestring(
    sf::st_coordinates(sf::st_geometry(x))[, 1:2]),
    dTolerance = units::set_units(100, "m")), 
    crs = sf::st_crs(x))

}
