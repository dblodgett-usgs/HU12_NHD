post_proc <- function(v2_net, hr_net, matched_lp) {
  
  matched_lp <- matched_lp %>%
    select(mr_LevelPathI, member_NHDPlusID, Hydroseq) %>%
    distinct()
  
  v2_net <- filter(v2_net, LevelPathI %in% matched_lp$mr_LevelPathI)
  
  hr_net <- filter(hr_net, COMID %in% matched_lp$member_NHDPlusID)
  
  matched_lp <- left_join(matched_lp, 
                          select(st_set_geometry(hr_net, NULL), 
                                 COMID, hydroseq = Hydroseq), 
                          by = c("member_NHDPlusID" = "COMID"))
  
  outlets <- group_by(matched_lp, mr_LevelPathI) %>%
    filter(Hydroseq == min(Hydroseq)) %>%
    ungroup() %>% select(outlet_NHDPlusID = member_NHDPlusID, mr_LevelPathI)
  
  matched_lp <- left_join(matched_lp, outlets, by = "mr_LevelPathI")
  
  hr_outlets <- hr_net %>%
    filter(COMID %in% outlets$outlet_NHDPlusID) %>%
    select(NHDPlusID = COMID) %>%
    get_mid_coord()
  
  st_crs(hr_outlets) <- st_crs(hr_net)
  
  gc()
  
  v2_outlets <- v2_net %>%
    select(COMID, LevelPathI, Hydroseq) %>%
    filter(LevelPathI %in% matched_lp$mr_LevelPathI) %>%
    group_by(LevelPathI) %>% 
    filter(Hydroseq == min(Hydroseq)) %>%
    ungroup() %>% select(LevelPathI)
  
  v2_outlets <- get_bottom_coord(v2_outlets)
  
  hr_outlets <- hr_outlets %>%
    left_join(distinct(select(matched_lp, NHDPlusID = outlet_NHDPlusID, mr_LevelPathI)), 
              by = "NHDPlusID")
  
  v2_outlets <- v2_outlets %>%
    filter(LevelPathI %in% hr_outlets$mr_LevelPathI)
  
  v2_outlets <- bind_cols(st_set_geometry(v2_outlets, NULL),
                          st_coordinates(v2_outlets) %>%
                            as.data.frame())
  
  hr_outlets <- st_transform(hr_outlets, st_crs(v2_net))
  
  hr_outlets <- bind_cols(st_set_geometry(hr_outlets, NULL),
                          st_coordinates(hr_outlets) %>%
                            as.data.frame())
  
  comp <- left_join(
    rename(v2_outlets, v2_X = X, v2_Y = Y),
    rename(hr_outlets, hr_X = X, hr_Y = Y),
    by = c("LevelPathI" = "mr_LevelPathI")) %>% 
    mutate(distance = sqrt((hr_X - v2_X)^2 + (hr_Y - v2_Y)^2))
  
  comp <- left_join(comp,
                    left_join(
                      select(st_set_geometry(v2_net, NULL), LevelPathI, LENGTHKM) %>%
                        group_by(LevelPathI) %>%
                        summarise(mr_lengthkm = sum(LENGTHKM)), 
                      select(st_set_geometry(hr_net, NULL), NHDPlusID = COMID, LengthKM = LENGTHKM) %>%
                        left_join(select(matched_lp, member_NHDPlusID, mr_LevelPathI), by = c("NHDPlusID" = "member_NHDPlusID")) %>%
                        group_by(mr_LevelPathI) %>%
                        summarise(hr_lengthkm = sum(LengthKM)), 
                      by = c("LevelPathI" = "mr_LevelPathI")),
                    by = "LevelPathI"
  )
  
  return(comp)
}

write_output <- function(compare, file) {
  readr::write_csv(compare, file)
  return(compare)
}

# for the HR flowpaths that are shorter than the MR flowpaths, 
# should find the network location where HR ends and should keep going.

# for the MR flowpaths that are shorter than the HR flowpaths, need to characterize still...

compare_table <- function(comp, hr_net, v2_net) {
  
  comp <- left_join(comp, 
                    select(st_set_geometry(hr_net, NULL), NHDPlusID = COMID, hr_LevelPathI = LevelPathI),
                    by = "NHDPlusID")
  
  comp <- select(comp, mr_LevelPathI = LevelPathI, 
                 hr_LevelPathI,
                 outlet_NHDPlusID = NHDPlusID, 
                 outlet_distance = distance, 
                 mr_lengthkm,
                 hr_lengthkm)
  
  DA <- st_set_geometry(v2_net, NULL) %>%
    filter(LevelPathI %in% comp$mr_LevelPathI) %>%
    select(LevelPathI, Hydroseq, TotDASqKM) %>%
    group_by(LevelPathI) %>%
    filter(Hydroseq == min(Hydroseq)) %>%
    ungroup() %>%
    select(-Hydroseq)
  
  comp <- left_join(comp, DA, by = c("mr_LevelPathI" = "LevelPathI")) %>%
    rename(mr_da = TotDASqKM)
  
  comp <- mutate(comp, mr_hr_lengthkm_diff = mr_lengthkm - hr_lengthkm, 
                 outlet_distance = outlet_distance / 1000) %>%
    rename(outlet_distance_km = outlet_distance)
  
  return(comp)
  
}

post_analysis <- function(comp, lower_thresh_km, upper_thresh_km) {
  
  comp_big <- filter(comp, mr_da > 1000)
  
  outlet_histo <- hist(comp_big$outlet_distance_km)
  length_hist <- hist(comp_big$mr_lengthkm - comp_big$hr_lengthkm)
  
  labs <- sapply(seq_len(length(outlet_histo$breaks) - 1), function(x) {
    paste(outlet_histo$breaks[x], "-", outlet_histo$breaks[x+1])
  })
  outlet_summary <- data.frame(distance_m = labs, count = outlet_histo$counts, stringsAsFactors = FALSE)
  
  labs <- sapply(seq_len(length(length_hist$breaks) - 1), function(x) {
    paste(length_hist$breaks[x], "-", length_hist$breaks[x+1])
  })
  length_summary <- data.frame(length_diff_m = labs, count = length_hist$counts)
  
  readr::write_csv(outlet_summary, "report/outlet_summary.csv")
  readr::write_csv(length_summary, "report/length_summary.csv")
  
  make_histo(comp, "out/report/all_outlet.png", "out/report/all_length.png", lower_thresh_km, upper_thresh_km)
  make_histo(comp_big, "out/report/big_outlet.png", "out/report/big_length.png", lower_thresh_km, upper_thresh_km)
}

plot_fun <- function(hr_net, v2_net, lp, comp, out_folder) {
  hr_net <- filter(hr_net, COMID %in% lp$member_NHDPlusID) %>%
    st_zm()
  v2_net <- filter(v2_net, LevelPathI %in% lp$mr_LevelPathI) %>%
    st_zm()
  
  gc()
  
  hr_net <- st_transform(hr_net, 3857) %>%
    st_simplify(dTolerance = 100)
  
  v2_net <- st_transform(v2_net, 3857) %>%
    st_simplify(dTolerance = 100)
  
  gc()
  
  bad <- comp[comp$outlet_distance_km > 20 & comp$outlet_distance_km < 800, ]
  
  if(nrow(bad) > 0) {
    lapply(1:nrow(bad), plot_comp, 
           bad = bad, lp = lp, hr_net = hr_net, 
           v2_net = v2_net, out_folder = out_folder)
  }
  
  return("done")
}

plot_comp <- function(row, bad, lp, hr_net, v2_net, states, out_folder) {
  mr_lp <- bad[row, ]$mr_LevelPathI
  
  hr <- filter(hr_net, COMID %in% lp$member_NHDPlusID[lp$mr_LevelPathI == mr_lp])
  mr <- filter(v2_net, LevelPathI == mr_lp)
  
  hr_length <- sum(hr$LengthKM)
  mr_length <- sum(mr$LENGTHKM)
  
  mr_bb <- st_bbox(mr)
  hr_bb <- st_bbox(hr)
  bbox <- c(min(mr_bb[1], hr_bb[1]), min(mr_bb[2], hr_bb[2]), 
            max(mr_bb[3], hr_bb[3]), max(mr_bb[4], hr_bb[4]))
  
  names(bbox) <- names(mr_bb)
  class(bbox) <- "bbox"
  length_scale <- bbox[4] - bbox[2]
  if(length_scale > 100000) {
    zoom = 9
  } else if(length_scale > 60000) {
    zoom = 10
  } else if(length_scale > 30000) {
    zoom = 11
  } else {
    zoom = 12
  }
  
  
  bbox <- st_as_sfc(bbox)
  st_crs(bbox) <- st_crs(mr)
  
  bgmap <- bbox %>%
    st_transform(4326) %>%
    st_bbox() %>%
    setNames(c("left", "bottom", "right", "top")) %>%
    ggmap::get_stamenmap(maptype = "terrain", zoom = zoom)
  
  png(file.path(out_folder, paste0(round(bad[row, ]$outlet_distance_km), "_", mr_lp, ".png")))
  plot(st_transform(bbox, 3857), 
       border = NA, 
       main = paste("Distance", bad[row, ]$outlet_distance_km, "km"), 
       sub = paste("MR Length:", mr_length, "km, HR Length:", hr_length, "km \n"),
       bgMap = bgmap)
  plot(st_geometry(mr), col = "blue", lwd = 4, add = TRUE)
  plot(st_geometry(hr), add = TRUE, col = "red", lwd = 2)
  usr <- par()$usr
  x_scl <- (usr[2] - usr[1])/100
  y_scl <- (usr[4] - usr[3])/100
  mr_l_y <- usr[3] + 5 * y_scl
  hr_l_y <- usr[3] + 5 * y_scl
  mr_l_x <- usr[1] + 20 * x_scl
  hr_l_x <- usr[1] + 60 * x_scl
  l_l <- 5 * x_scl
  segments(mr_l_x, mr_l_y, mr_l_x + l_l, mr_l_y, col = "blue", lwd = 2)
  segments(hr_l_x, hr_l_y, hr_l_x + l_l, hr_l_y, col = "red", lwd = 2)
  text(mr_l_x + 3 * x_scl, mr_l_y, "NHDPlusV2", pos = 4)
  text(hr_l_x + 3 * x_scl, hr_l_y, "NHDPlusHR", pos = 4)
  dev.off()
}

make_histo <- function(comp, out_file_outlet, out_file_length, lower_thresh_km, upper_thresh_km) {
  png(out_file_outlet, width = 512, height = 384)
  plot_data <- comp$outlet_distance_km[comp$outlet_distance_km > lower_thresh_km & 
                                         comp$outlet_distance_km < upper_thresh_km]
  lower <- sum(comp$outlet_distance_km < lower_thresh_km)
  upper <- sum(comp$outlet_distance_km > upper_thresh_km)
  outlet_histo <- hist(log10(plot_data), 
                       breaks = 20, 
                       main = paste("Distance between", length(plot_data), "MR and HR outlets more \n", 
                                    "than 1km and less than 100km apart \n",
                                    lower, "truncated below", 
                                    paste0(lower_thresh_km, "km."), upper, "truncated above", 
                                    paste0(upper_thresh_km, "km.")),
                       xaxt = "n", xlab = "Distance in km")
  axis(1, at = pretty(outlet_histo$breaks), labels = round(10^pretty(outlet_histo$breaks)))
  dev.off()
  
  png(out_file_length, width = 512, height = 384)
  l_diff <- comp$mr_lengthkm - comp$hr_lengthkm
  plot_data <- l_diff[l_diff > lower_thresh_km & l_diff < upper_thresh_km]
  lower <- sum(l_diff < lower_thresh_km)
  upper <- sum(l_diff > upper_thresh_km)
  length_hist <- hist(log10(plot_data),
                      breaks = 20, 
                      main = paste("Difference in length between", length(plot_data), "MR and HR \n", 
                                   "flowpaths more than 1km and less than 100km different. \n",
                                   lower, "truncated below", 
                                   paste0(lower_thresh_km, "km."), upper, "truncated above", 
                                   paste0(upper_thresh_km, "km.")),
                      xaxt = "n", xlab = "Distance in km")
  axis(1, at = pretty(length_hist$breaks), labels = round(10^pretty(length_hist$breaks)))
  dev.off()
}

get_mid_coord <- function(cats) {
  bind_cols(st_set_geometry(cats, NULL), 
            cats %>%
              st_coordinates() %>%
              as.data.frame() %>%
              group_by(L2) %>%
              filter(row_number() == round(n()/2)) %>%
              ungroup() %>%
              select(X, Y) %>%
              st_as_sf(coords = c("X", "Y"))) %>%
    st_sf()
}

get_bottom_coord <- function(cats) {
  bind_cols(st_set_geometry(cats, NULL), 
            cats %>%
              st_coordinates() %>%
              as.data.frame() %>%
              group_by(L2) %>%
              filter(row_number() == n()) %>%
              ungroup() %>%
              select(X, Y) %>%
              st_as_sf(coords = c("X", "Y"))) %>%
    st_sf()
}
