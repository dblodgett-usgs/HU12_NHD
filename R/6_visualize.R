get_wbd_matched <- function(all, wbd) {

  no_intersect <- filter(all, trib_no_intersect) %>%
    group_by(corrected_LevelPathI) %>%
    mutate(group_size = n()) %>%
    select(corrected_LevelPathI, group_size) %>%
    distinct()
  
  return(left_join(wbd, select(all, -TOHUC), by = "HUC12"))
}

get_wbd_grouped <- function(wbd) { 
  wbd_grouped <- group_by(wbd, corrected_LevelPathI) %>%
    select(corrected_LevelPathI) %>%
    summarise(do_union = FALSE)
}

write_output_gpkg <- function(net, wbd, hu_joiner, points, prj, viz_simp, out_dir) {
  
  wbd_viz_gpkg <- file.path(out_dir, "wbd_viz.gpkg")
  
  write_sf(wbd, wbd_viz_gpkg, "wbd_viz")
  
  write_sf(st_simplify(st_transform(net, prj), 
                       dTolerance = viz_simp), 
           wbd_viz_gpkg, "net")
  
  wbd_matched <- get_wbd_matched(hu_joiner, wbd)
  
  write_sf(wbd_matched, wbd_viz_gpkg, "wbd_matched")
  
  # write_sf(get_wbd_grouped(wbd_matched), wbd_viz_gpkg, "wbd_grouped")
  
  write_sf(points, wbd_viz_gpkg, "linked_points")
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