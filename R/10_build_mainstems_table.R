build_mainstem_table <- function(nhdplus_net, 
                                 olp, 
                                 nlp,
                                 rfp,
                                 rf1,
                                 nhdp_wbd,
                                 wbd) {
  
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

make_ms_summary <- function(ms, nhdp_att) {
  nhd_prep <- prepare_nhdplus(nhdp_att, 0, 0, 0, TRUE)
  
  area <- left_join(select(nhd_prep, ID = COMID, toID = toCOMID), select(nhdp_att, ID = COMID, area = AreaSqKM)) %>%
    calculate_total_drainage_area()
  
  ms_length <- nhdp_att %>%
    select(LevelPathI, LENGTHKM) %>%
    group_by(LevelPathI) %>%
    summarize(length = sum(LENGTHKM))
    
  summary <- data.frame(COMID = nhd_prep$COMID, totdasqkm = area) %>%
    left_join(select(nhdp_att, COMID, LevelPathI, 
                     level = StreamLeve, 
                     order = StreamOrde), by = "COMID") %>%
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
