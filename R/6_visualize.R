get_wbd_matched <- function(all, wbd) {

  no_intersect <- filter(all, trib_no_intersect) %>%
    group_by(corrected_LevelPathI) %>%
    mutate(group_size = n()) %>%
    select(corrected_LevelPathI, group_size) %>%
    distinct()
  
  if(is.numeric(all$HUC12)) all$HUC12 <- as.character(all$HUC12)
  
  return(left_join(wbd, select(all, -TOHUC), by = "HUC12"))
}

get_wbd_grouped <- function(wbd) { 
  wbd_grouped <- group_by(wbd, corrected_LevelPathI) %>%
    select(corrected_LevelPathI) %>%
    summarise(do_union = FALSE)
}

write_output_gpkg <- function(net, wbd, hu_joiner, points, prj, viz_simp, wbd_viz_gpkg) {
  
  if(names(net) == "NHDFlowline") net <- net$NHDFlowline
  
  write_sf(wbd, wbd_viz_gpkg, "wbd_viz")
  
  write_sf(st_simplify(st_transform(net, prj), 
                       dTolerance = viz_simp), 
           wbd_viz_gpkg, "net")
  
  wbd_matched <- get_wbd_matched(hu_joiner, wbd)
  
  write_sf(wbd_matched, wbd_viz_gpkg, "wbd_matched")
  
  # write_sf(get_wbd_grouped(wbd_matched), wbd_viz_gpkg, "wbd_grouped")
  
  write_sf(points, wbd_viz_gpkg, "linked_points")
}

write_rf1_output <- function(rf1, rf1_nhdplus, rf1_out) {
  ret <- select(rf1_nhdplus, ERF1_2. = member_ID, LevelPathI = mr_LevelPathI)
  readr::write_csv(ret,
                   file.path(rf1_out, "map_joiner.csv"))
  
  viz_gpkg <- file.path(rf1_out, "wbd_viz.gpkg")
  
  ret <- left_join(rf1, ret, by = c("ID" = "ERF1_2."))
  
  write_sf(ret, viz_gpkg)  
  
  ret
}

geom_plot_data <- function(hu12, net, lookup, filter) {
  
  lookup <- filter(lookup, grepl(filter, lookup$HUC12)) %>%
    select(-intersected_LevelPathI, -trib_intersect, -trib_no_intersect, -headwater_error, -TOHUC) %>%
    group_by(HUC12) %>%
    filter(corrected_LevelPathI == min(corrected_LevelPathI)) %>% # coastals with multiple get matched multiple times.
    ungroup()
  
  hu12 <- filter(hu12, grepl(filter, hu12$HUC12)) %>%
    left_join(lookup, by = "HUC12")
    
  hu_grouped <- get_wbd_grouped(hu12)
  hu_grouped <- filter(hu_grouped, corrected_LevelPathI %in% lookup$corrected_LevelPathI)
  hu12 <- filter(hu12, HUC12 %in% lookup$HUC12) %>%
    group_by(HUC12)
  
  hu12 <- dplyr::summarise(hu12, TOHUC = TOHUC[1])
  
  fromHUC <- sapply(hu12$HUC12, fromHUC_finder, 
                    hucs = hu12$HUC12, tohucs = hu12$TOHUC)
  aggrHUCs <- sapply(hu12$HUC12, HUC_aggregator, fromHUC = fromHUC)
  
  hu12_boundaries <- unionHUCSet(aggrHUCs, fromHUC, as_Spatial(hu12))
  
  net <- filter(net, LevelPathI %in% hu_grouped$corrected_LevelPathI)
  
  return(list(net = net, hu_grouped = hu_grouped, lookup = lookup, hu12_boundaries = hu12_boundaries))
}

create_png <- function(plot_data, hu_joiner, out_folder) {
  net <- plot_data$net
  hu_grouped <- plot_data$hu_grouped
  lookup <- plot_data$lookup
  
  hu12_boundaries <- st_as_sf(plot_data$hu12_boundaries)
  
  dir.create(out_folder, showWarnings = FALSE)
  unlink(paste0(out_folder, "*"))
  
  # Plot order:
  # 1) mainstem river and total upstream drainage boundary
  # 2) mainstem HU12s
  # 3) tributaries to mainstem with their drainage boundaries
  # 4) mainstem HU12s of tributaries added in 3
  # repeat till none left.
  
  lp_dnlp <- select(st_set_geometry(net, NULL), LevelPathI, DnLevelPat) %>%
    distinct() %>%
    filter(LevelPathI != DnLevelPat)
  
  # group_by(lookup, corrected_LevelPathI) %>% summarize(size = n())
  levelpath <- 290002124 # entire tombigbee
  # levelpath <- 290007998
  
  plotter <- function(i) {
    png(filename = paste0(out_folder, i, ".png"), width = 600, height = 768, units = "px")
    par(mar=c(0, 0, 0, 0))
    plot(st_geometry(viz_basinboundary), lwd = 0.1, border = "white", col = NA)
    plot(st_geometry(viz_mainstem), col = "blue", add = TRUE)
    plot(st_geometry(viz_hydrounits), lwd = 0.5, border = "grey", col = NA, add = TRUE)
    dev.off()
  }
  
  i <- 1
  viz_hydrounits <- filter(hu_grouped, corrected_LevelPathI %in% 0)
  mainstem_hus <- distinct(filter(lookup, corrected_LevelPathI %in% levelpath))
  viz_basinboundary <- filter(hu12_boundaries, HUC12 %in% unique(mainstem_hus$outlet_HUC12))
  
  for(j in c(1:4)) {
    viz_mainstem <- filter(net, LevelPathI %in% levelpath)
    mainstem_hus <- distinct(filter(lookup, corrected_LevelPathI %in% levelpath))
    plotter(i)
    i <- i + 1
    viz_hydrounits <- filter(hu_grouped, corrected_LevelPathI %in% levelpath)
    plotter(i)
    next_lp <- filter(lp_dnlp, DnLevelPat %in% levelpath)$LevelPathI
    levelpath <- unique(c(levelpath, next_lp))
    i <- i + 1
  }
  
  plotter <- function(i) {
    png(filename = paste0(out_folder, i, ".png"), width = 450, height = 768, units = "px")
    par(mar=c(0, 0, 0, 0))
    plot(st_geometry(viz_basinboundary), lwd = 0.5)
    plot(st_geometry(viz_mainstem), col = "blue", add = TRUE)
    plot(st_geometry(viz_hydrounits), lwd = 0.5, border = "grey", col = NA, add = TRUE)
    dev.off()
  }
  
  joiner <- hu_joiner %>%
    filter(grepl("^03.*", HUC12)) %>%
    select(HUC12, intersected_LevelPathI, corrected_LevelPathI)
  
  intersected <- left_join(select(hu12_boundaries, HUC12), joiner)
  
  levelpath <- 290003311
  viz_hydrounits <- filter(hu_grouped, corrected_LevelPathI %in% 0)
  viz_mainstem <- filter(net, LevelPathI %in% levelpath)
  mainstem_hus <- distinct(filter(lookup, corrected_LevelPathI %in% levelpath))
  viz_basinboundary <- filter(hu12_boundaries, HUC12 %in% unique(mainstem_hus$outlet_HUC12))
  
  plotter("fig1")
  intersected_hydrounits <- filter(intersected, intersected_LevelPathI %in% levelpath)
  corrected_hydrounits <- filter(intersected, corrected_LevelPathI %in% levelpath)
  
  png(filename = paste0(out_folder, "fig2.png"), width = 450, height = 768, units = "px")
  par(mar=c(0, 0, 0, 0))
  plot(st_geometry(viz_basinboundary), lwd = 0.5)
  plot(st_geometry(intersected_hydrounits), lwd = 0.1, border = "grey", col = NA, add = TRUE)
  plot(st_geometry(viz_mainstem), col = "blue", add = TRUE)
  plot(st_geometry(intersected_hydrounits), lwd = 0.5, border = "grey", col = NA, add = TRUE)
  dev.off()
  
  head_hu <- unique(filter(lookup, corrected_LevelPathI == levelpath)$head_HUC12)
  viz_head <- filter(hu12_boundaries, HUC12 == head_hu)
  
  png(filename = paste0(out_folder, "fig3.png"), width = 450, height = 768, units = "px")
  par(mar=c(0, 0, 0, 0))
  plot(st_geometry(viz_basinboundary), lwd = 0.5)
  plot(st_geometry(intersected_hydrounits), lwd = 0.1, border = "grey", col = NA, add = TRUE)
  plot(st_geometry(viz_mainstem), col = "blue", add = TRUE)
  plot(st_geometry(intersected_hydrounits), lwd = 0.5, border = "grey", col = NA, add = TRUE)
  plot(st_geometry(viz_head), lwd = 0.5, border = "red", col = "NA", add = TRUE)
  dev.off()
  
  png(filename = paste0(out_folder, "fig4.png"), width = 450, height = 768, units = "px")
  par(mar=c(0, 0, 0, 0))
  plot(st_geometry(viz_basinboundary), lwd = 0.5)
  plot(st_geometry(intersected_hydrounits), lwd = 0.5, border = "grey", col = NA, add = TRUE)
  plot(st_geometry(viz_mainstem), col = "blue", add = TRUE)
  plot(st_geometry(corrected_hydrounits), lwd = 0.5, border = "black", col = NA, add = TRUE)
  plot(st_geometry(viz_head), lwd = 0.5, border = "red", col = "NA", add = TRUE)
  dev.off()
  
  next_lp <- filter(lp_dnlp, DnLevelPat %in% levelpath)$LevelPathI
  levelpath <- unique(c(levelpath, next_lp))
  viz_mainstem <- filter(net, LevelPathI %in% levelpath)
  
  png(filename = paste0(out_folder, "fig5.png"), width = 450, height = 768, units = "px")
  par(mar=c(0, 0, 0, 0))
  plot(st_geometry(viz_basinboundary), lwd = 0.5)
  plot(st_geometry(intersected_hydrounits), lwd = 0.5, border = "grey", col = NA, add = TRUE)
  plot(st_geometry(corrected_hydrounits), lwd = 0.5, border = "black", col = NA, add = TRUE)
  plot(st_geometry(viz_mainstem), col = "blue", add = TRUE)
  dev.off()
  
  #######
  
  
  plotter <- function(i) {
    png(filename = paste0(out_folder, i, ".png"), width = 512, height = 768, units = "px")
    plot(st_geometry(viz_basinboundary), lwd = 0.1, border = "grey", col = NA)
    plot(st_geometry(viz_mainstem), col = "blue", add = TRUE)
    plot(st_geometry(viz_hydrounits), lwd = 0.5, border = "grey", col = NA, add = TRUE)
    dev.off()
  }
  
  levelpath <- 290001792 # ACF
  viz_hydrounits <- filter(hu_grouped, corrected_LevelPathI %in% 0)
  viz_mainstem <- filter(net, LevelPathI %in% levelpath)
  mainstem_hus <- distinct(filter(lookup, corrected_LevelPathI %in% levelpath))
  viz_basinboundary <- filter(hu12_boundaries, HUC12 %in% unique(mainstem_hus$outlet_HUC12))
  plotter("example1")
  next_lp <- filter(lp_dnlp, DnLevelPat %in% levelpath)$LevelPathI
  levelpath <- unique(c(levelpath, next_lp))
  viz_mainstem <- filter(net, LevelPathI %in% levelpath)
  mainstem_hus <- distinct(filter(lookup, corrected_LevelPathI %in% levelpath))
  viz_basinboundary <- filter(hu12_boundaries, HUC12 %in% unique(mainstem_hus$outlet_HUC12))
  plotter("example2")
  next_lp <- filter(lp_dnlp, DnLevelPat %in% levelpath)$LevelPathI
  levelpath <- unique(c(levelpath, next_lp))
  viz_mainstem <- filter(net, LevelPathI %in% levelpath)
  mainstem_hus <- distinct(filter(lookup, corrected_LevelPathI %in% levelpath))
  viz_basinboundary <- filter(hu12_boundaries, HUC12 %in% unique(mainstem_hus$outlet_HUC12))
  plotter("example3")
}

get_wbd_plot_data <- function(nhdplus_net, wbd_gdb_path, plumbing, viz_simp, prj, cores, out_gpkg) {
  if(file.exists(out_gpkg)) {
    message(out_gpkg, " already exists, skipping.")
    return(out_gpkg)
  }
  
  prep <- prepare_nhdplus(nhdplus_net, 0, 0, 0, FALSE)$COMID
  
  nhdplus_net <- select(nhdplus_net, COMID, LENGTHKM, LevelPathI, Pathlength, DnLevelPat, Hydroseq, DnHydroseq, UpHydroseq)
  
  net_atts <- st_set_geometry(nhdplus_net, NULL)
  
  hu_levels <- c(2, 4, 6, 8, 10)
  
  for(hu_level in hu_levels) {
    hu_att <- paste0("HUC", hu_level)
    gdb_layer <- paste0("WBDHU", hu_level)
    
    wbd <- read_sf(wbd_gdb_path, gdb_layer)
    wbd <- wbd[!grepl("^20.*|^19.*|^21.*|^22.*", wbd[[hu_att]]), ]
    wbd <- st_transform(wbd, prj)
    
    outlets <- plumbing[plumbing[[hu_att]] %in% wbd[[hu_att]], ]  
    
    outlet_ids <- unique(outlets$COMID)
    
    cl <- parallel::makeCluster(rep("localhost", cores), 
                                type = "SOCK", 
                                outfile = "viz.log")
    
    all <- c(do.call(c, parLapply(cl, outlet_ids, dm_fun, network = net_atts)), 
             do.call(c, parLapply(cl, outlet_ids, um_fun, network = net_atts)))
    
    stopCluster(cl)
    
    nt <- nhdplus_net[nhdplus_net[["COMID"]] %in% all &
                                 nhdplus_net[["COMID"]] %in% prep, ] %>%
      st_transform(prj)
    
    bb <- st_bbox(wbd)
    
    write_sf(wbd, out_gpkg, paste0("wbd_hu", hu_level))
    write_sf(nt, out_gpkg, paste0("nhd_hu", hu_level))
    write_sf(outlets, out_gpkg, paste0("out_hu", hu_level))
    
    # plot_fun(bb, st_geometry(wbd), nt, st_geometry(outlets), 
    #          paste0("png/wbd_level", hu_level, ".png"))
  }
  return(out_gpkg)
}


plot_wbd <- function(gpkg) {
  hu_levels <- c(2, 4, 6, 8, 10)
  
  for(hl in hu_levels) {
    out <- read_sf(gpkg, paste0("out_hu", hl))
    net <- read_sf(gpkg, paste0("nhd_hu", hl)) %>%
      rmapshaper::ms_simplify(keep = 0.05, sys = TRUE)
    wbd <- read_sf(gpkg, paste0("wbd_hu", hl)) %>%
      rmapshaper::ms_simplify(keep = 0.05, sys = TRUE)
    bb <- st_as_sfc(st_bbox(wbd), crs = st_crs(wbd))
    plot_fun(bb, st_geometry(wbd), st_geometry(net), st_geometry(out), paste0("out/png/wbd_viz_hu", hl))
  }
}

plot_fun <- function(bb, wb, nt, pt, fi) {
  png(fi, width = 4096, height = 3072)
  plot(bb, border = NA)
  plot(wb, lwd = .7, add = TRUE)
  plot(nt, col = "blue", add = TRUE)
  plot(pt, col = "black", pch = 16, add = TRUE)
  dev.off()
}

dm_fun <- function(id, network) {
  nhdplusTools::get_DM(network = network, comid = id)
}

um_fun <- function(id, network) {
  nhdplusTools::get_UM(network = network, comid = id)
}

get_lp_level <- function(order_f, ms, max_level) {
  level_f <- max_level - order_f + 1
  filter(ms, !is.na(outlet_nhdpv2HUC12) & 
           !is.na(outlet_GNIS_ID) & 
           order >= order_f &
           level <= level_f)
}

get_lp_data <- function(i, ms_lp, dat) {
  filter(dat, LevelPathI %in% ms_lp[[i]]$LevelPathI)
}

get_lp_plot_data <- function(nhdp_net, ms, super_simp) {
  ms <- filter(ms, !is.na(outlet_nhdpv2HUC12))

  nhdp_net <- select(nhdp_net, LevelPathI)
  nhdp_net <- nhdp_net[nhdp_net$LevelPathI %in% ms$LevelPathI, ]
  
  ms <- lapply(c(10:1), get_lp_level, ms = ms, max_level = 10)

  message("simplifying")
  nhdp_net <- st_simplify(nhdp_net, dTolerance = super_simp)
  
  list(bbox = st_as_sfc(st_bbox(nhdp_net), crs = st_crs(nhdp_net)), 
       lp_data = lapply(c(1:10), get_lp_data, 
                        dat = nhdp_net, ms_lp = ms),
       ms = ms)
}

get_lp_plot_data_wbd <- function(pd, wbd, wbd_lp, super_simp) {
  wbd <- left_join(wbd, select(st_set_geometry(wbd_lp, NULL), HUC12, LevelPathI), by = "HUC12")
  
  ms <- pd$ms
  
  wbd <- st_simplify(wbd, dTolerance = super_simp)
  
  c(pd, list(wbd = lapply(c(1:10), get_lp_data, 
                          dat = wbd, ms_lp = ms)))
}

get_lp_plot_data_rf1 <- function(pd, rf1, rf1_lp, super_simp) {
  rf1 <- left_join(rf1, select(rf1_lp, -headwater_COMID, LevelPathI = mr_LevelPathI), by = c("ID" = "member_ID"))
  
  ms <- pd$ms
  
  rf1 <- st_simplify(rf1, dTolerance = super_simp)
  
  c(pd, list(rf1 = lapply(c(1:10), get_lp_data, 
                          dat = rf1, ms_lp = ms)))
}

get_lp_plots <- function(lp_plot_data, delay, hu02, hu02_filter = NULL, bb = NULL) {
  
  if(!is.null(bb))  {
    class(bb) <- "bbox"
    attributes(bb)$crs <- st_crs(4326)
    lp_plot_data$bbox <- st_transform(st_as_sfc(bb, crs = st_crs(4326)), st_crs(hu02))
  } else if(!is.null(hu02_filter)) {
    lp_plot_data$bbox <- st_as_sfc(st_bbox(filter(hu02, HUC2 %in% hu02_filter)), crs = st_crs(hu02))
  }
  
  plotter <- function(i, bbox, lp_plot_data, hu02) {
    par(mar=c(0,0,0,0))
    plot(bbox, border = NA)
    plot(st_geometry(hu02), lwd = 1, col = NA, add = TRUE)
    plot(st_geometry(lp_plot_data$wbd[[i]]), lwd = .7, col = "lightgrey", add = TRUE)
    plot(st_geometry(lp_plot_data$lp_data[[i]]), lwd = 2.5, col = "blue", add = TRUE)
    plot(st_geometry(lp_plot_data$rf1[[i]]), lwd = 1.25, col = "red", add = TRUE)
  }
  
  gifski::save_gif({
    lapply(c(2:6), plotter, 
           bbox = lp_plot_data$bbox, 
           lp_plot_data = lp_plot_data, 
           hu02 = hu02)
  }, 
  gif_file = "out/gif/big_animation.gif", 
  width = 2048, height = 1536, delay = delay)
  
}

get_hw_fig <- function() {
  rf1 <- sf::read_sf(system.file("extdata/rf1_test.gpkg", package = "mainstems"), "rf1_test")
  nhdp <- sf::read_sf(system.file("extdata/rf1_test.gpkg", package = "mainstems"), "nhdp_flowline")
  hr <- nhdplusTools::download_nhdplushr("data/hr/", hu_list = "03")
  hr <- nhdplusTools::get_nhdplushr(hr_dir = hr, out_gpkg = "data/0314.gpkg", pattern = "0314_HU4_GDB.gdb")
  hr <- st_zm(read_sf(hr, "NHDFlowline"))
  hr_outlet <- hr$HydroSeq[which(hr$GNIS_ID == "00288720" & hr$TerminalFl == 1)]
  hr_sub <- hr[hr$TerminalPa == hr_outlet, ]
  hr_sub <- st_transform(hr_sub, st_crs(nhdp))
  
  hw_pairs <- structure(list(ID = c(7849L, 7851L, 7853L, 7856L, 7857L, 7858L, 
                                    7860L, 7862L, 7863L, 7864L), 
                             FEATUREID = c(2171225L, 2171903L, 2171013L, 
                                           2171859L, 2171883L, 2171917L, 
                                           2172119L, 2171999L, 2172059L, 
                                           2171563L)), 
                        row.names = c(NA, -10L), 
                        class = c("tbl_df", "tbl", "data.frame"))
  
  rf1_matches <- match_flowpaths(nhdp, rf1, hw_pairs)
  lps <- unique(rf1_matches$mr_LevelPathI)[1:9]
  
  bb <- c(xmin = 789495, ymin = 876900, 
          xmax = 794000, ymax = 898900)
  class(bb) <- "bbox"
  bb <- st_as_sfc(bb, crs = st_crs(nhdp))
  
  png("img/fig_hw.png", 500, 700)
  par(mar = c(0,0,0,0))
  plot(bb, border = NA)
  plot(st_geometry(hr_sub), lwd = 3, add = TRUE)
  plot(st_geometry(nhdp), lwd = 1.5, col = "blue", add = TRUE)
  plot(st_geometry(filter(nhdp, LevelPathI == 290016401)),
       col = "blue", lwd = 4, add = TRUE)
  plot(st_geometry(filter(rf1, ID %in% filter(rf1_matches, mr_LevelPathI %in% 290016401)$member_ID)), 
       col = "red", lwd = 3, add = TRUE)
  dev.off()
}

