get_plot_data <- function(nhdp_net, hu02, order, proj, simp_meters) {
  pattern <- paste(paste0("^", hu02, ".*"), collapse = "|")
  filter(nhdp_net, 
         grepl(pattern, nhdp_net$REACHCODE) & 
           nhdp_net$StreamOrde >= order) %>%
    st_zm() %>%
    select(COMID) %>%
    st_transform(proj) %>%
    st_simplify(dTolerance = simp_meters)
}

get_v2_plot_data <- function(nhdp_net, lp, proj, simp_meters) {
  filter(nhdp_net, LevelPathI %in% lp$mr_LevelPathI) %>%
    st_zm() %>%
    select(COMID, LevelPathI) %>%
    st_transform(proj) %>%
    st_simplify(dTolerance = simp_meters)
}

get_hr_plot_data <- function(hr_net, lp, proj, simp_meters) {
  filter(hr_net, NHDPlusID %in% lp$member_NHDPlusID) %>%
    st_zm() %>%
    select(NHDPlusID) %>%
    st_transform(proj) %>%
    st_simplify(dTolerance = simp_meters)
}

create_plot_frames <- function(big_plot_data, v2_plot_data, 
                               hr_plot_data, matched_lp,
                               hr_vpu, wbd_hr, threads) {

  hu04 <- regexpr("[0-9][0-9][0-9][0-9]", hr_vpu)
  hu04 <- substr(hr_vpu, hu04, hu04 + attr(hu04, "match.length") - 1)
  hu_04 <- read_sf(wbd_hr) %>%
    filter(HUC4 == hu04) %>%
    st_transform(3857) %>%
    st_simplify(dTolerance = 200)
  
  hr_plot_data <- st_sf(hr_plot_data)
  big_plot_data <- st_sf(big_plot_data)
  
  big_plot_data <- st_intersection(big_plot_data, hu_04$geometry)
  v2_plot_data <- st_intersection(v2_plot_data, hu_04$geometry)
  
  lp_order <- matched_lp %>%
    group_by(mr_LevelPathI) %>%
    summarise(count = n()) %>%
    arrange(desc(count))
  
  lp <- lapply(seq_along(lp_order$mr_LevelPathI), function(p, lp) lp[1:p], 
               lp <- lp_order$mr_LevelPathI)
  
  title <- paste0("Hydrologic Unit ", hu04, ": ", hu_04$NAME, "\n NHDPlus V2 and High Res")
  subtitle <- ""
  
  gifski::save_gif(
    expr = lapply(lp[1:100], plot_it_up, hu_04 = hu_04, title = title, subtitle = subtitle, 
                  matched_lp = matched_lp, big_plot_data = big_plot_data, 
                  hr_plot_data = hr_plot_data, v2_plot_data = v2_plot_data),
    gif_file = paste0(hu04, ".gif"), delay = 0.3, loop = TRUE, progress = TRUE)
}

plot_it_up <- function(lp, hu_04, title = "Title", subtitle = "subtitle", 
                       matched_lp, big_plot_data, 
                       hr_plot_data, v2_plot_data) {
  hr_ids <- matched_lp$member_NHDPlusID[matched_lp$mr_LevelPathI %in% lp]
  
  hr <- hr_plot_data[hr_plot_data$NHDPlusID %in% hr_ids, ]
  mr <- v2_plot_data[v2_plot_data$LevelPathI %in% lp, ]
  
  par(mar = c(2, 0, 4, 0), xpd=TRUE)
  plot(st_geometry(hu_04),
       main = paste(title), 
       sub = paste(subtitle))
  
  plot(st_geometry(big_plot_data), col = "blue", lwd = 1, add = TRUE)
  
  plot(st_geometry(mr), col = "blue4", lwd = 3, add = TRUE)  
  plot(st_geometry(hr), col = "red", lwd = 1, add = TRUE)
  
  usr <- par()$usr
  x_scl <- (usr[2] - usr[1])/100
  y_scl <- (usr[4] - usr[3])/100
  mr_l_y <- usr[3] - 2 * y_scl
  hr_l_y <- usr[3] - 2 * y_scl
  mr_l_x <- usr[1] + 35 * x_scl
  hr_l_x <- usr[1] + 55 * x_scl
  l_l <- 3 * x_scl
  segments(mr_l_x, mr_l_y, mr_l_x + l_l, mr_l_y, col = "blue4", lwd = 3)
  segments(hr_l_x, hr_l_y, hr_l_x + l_l, hr_l_y, col = "red", lwd = 1)
  text(mr_l_x + 3 * x_scl, mr_l_y + 2, "NHDPlusV2", pos = 4)
  text(hr_l_x + 3 * x_scl, hr_l_y + 2, "NHDPlusHR", pos = 4)
  
}

get_bb <- function(mr = NA, hr = NA, geo = NA) {
  if(!is.na(geo)) {
    bbox <- st_bbox(geo)
    crs <- st_crs(geo)
  } else {
    mr_bb <- st_bbox(mr)
    hr_bb <- st_bbox(hr)
    bbox <- c(min(mr_bb[1], hr_bb[1]), min(mr_bb[2], hr_bb[2]), 
              max(mr_bb[3], hr_bb[3]), max(mr_bb[4], hr_bb[4]))
    
    names(bbox) <- names(mr_bb)
    class(bbox) <- "bbox"
    crs <- st_crs(mr)
  }
  x_size <- bbox$xmax - bbox$xmin
  y_size <- bbox$ymax - bbox$ymin
  if(x_size > y_size) {
    diff <- x_size - y_size
    bbox$xmin <- bbox$xmin - diff
    bbox$xmax <- bbox$xmax + diff
  } else {
    diff <- y_size - x_size
    bbox$xmin <- bbox$ymin - diff
    bbox$xmax <- bbox$ymax + diff
  }
  bbox <- st_as_sfc(bbox)
  st_crs(bbox) <- crs
  return(bbox)
}
