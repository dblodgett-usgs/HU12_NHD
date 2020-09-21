#' Get linked points
#' @description Given matched levelpaths, finds outlets of catchments along level paths.
#' @param hu_lp data.frame as returned by `match_levelpaths`
#' @param net sf data.frame NHDPlus Flowlines
#' @param wbd  sf data.frame HU12s
#' @param exclude character vector of HUs to exclude from consideration
#' @param cores numeric (optional) number of cores to use in paralel evaluation
#'
#' This function is preliminary and subject to revision. It has been tested thoroughly
#' but complete methods description have not yet been published.
#'
#' @importFrom dplyr bind_rows bind_cols
#' @importFrom sf st_coordinates st_as_sf st_crs<- st_transform st_point st_geometry<- st_intersection
#' @importFrom stats setNames
#' @export
#'
#' @examples
#' library(dplyr)
#' library(sf)
#' library(nhdplusTools)
#'
#' source(system.file("extdata/new_hope_data.R", package = "nhdplusTools"))
#' suppressWarnings(net_prep <- select(new_hope_flowline, COMID, DnLevelPat, 
#'                    AreaSqKM) %>%
#'                    right_join(prepare_nhdplus(new_hope_flowline, 20, 0, 10), 
#'                    by = "COMID") %>%
#'                    st_sf() %>%
#'                    group_by(LevelPathI) %>%
#'                    arrange(Hydroseq) %>%
#'                    mutate(DnLevelPat = DnLevelPat[1]) %>%
#'                    ungroup())
#'
#' net_prep["denTotalAreaSqKM"] <-
#'   nhdplusTools::calculate_total_drainage_area(select(st_set_geometry(net_prep, NULL),
#'                                        ID = COMID, toID = toCOMID,
#'                                        area = AreaSqKM))
#'
#' wbd <- select(new_hope_wbd, HUC12 = HUC_12, TOHUC = HU_12_DS) %>%
#'   st_transform(st_crs(net_prep))
#'
#' net_prep <- st_join(net_prep, wbd) %>%
#'   st_set_geometry(NULL)
#'
#' outlet_comid <- 8897784
#'
#' hu_lp <- match_levelpaths(net_prep, outlet_comid)
#'
#' linked_points <- get_linked_points(hu_lp, new_hope_flowline, wbd, exclude = c())
#'
#' lps <- filter(new_hope_flowline, LevelPathI %in% linked_points$LevelPathI)
#'
#' plot(st_geometry(new_hope_flowline), lwd = 0.4, col = "blue")
#' plot(st_geometry(lps), col = "blue", add = TRUE)
#' plot(linked_points$geometry, col = "red", pch = 19, add = TRUE)
#' plot(st_transform(new_hope_wbd$geom,
#'                   st_crs(linked_points)), lwd = 1.5, add = TRUE)

get_linked_points <- function(hu_lp, net, wbd, exclude, cores = NA) {
  
  if("NHDFlowline" %in% names(net)) net <- net$NHDFlowline
  
  lp_points <- get_lp_points(hu_lp, net, wbd, exclude)
    
  na_points <- lp_points$na
  lp_points <- lp_points$lp

  na_outlet_coords <- get_na_outlets_coords(na_points, net)
    
  lp_list <- unique(lp_points$lp)
  
  net <- select(net, .data$COMID, .data$LevelPathI, .data$REACHCODE, .data$ToMeas, 
                .data$FromMeas, .data$Hydroseq) %>%
    filter(.data$LevelPathI %in% lp_list)

  in_list <- lapply(lp_list, in_list_fun, net = net, lp_points = lp_points)
  
  return(get_linked_points_scalable(in_list, na_outlet_coords, cores, check_file = NULL))
}

#' Get LevelPath Points
#' @description finds 
#' @param hu_lp data.frame as returned by `match_levelpaths`
#' @param net sf data.frame ...
#' @param wbd  sf data.frame ...
#' @param exclude character vector of HUs to exclude from consideration
#' @export
get_lp_points <- function(hu_lp, net, wbd, exclude) {
  if("NHDFlowline" %in% names(net)) net <- net$NHDFlowline
  
  lp_points <- get_points_out(hu_lp, net, wbd, exclude) %>%
    get_lp_hu_points(st_crs(wbd)) %>%
    filter(!.data$hu12 %in% exclude) %>%
    mutate(lp = as.numeric(.data$lp)) %>%
    group_by(.data$hu12) %>%
    filter(.data$lp == min(.data$lp)) %>%
    ungroup()
  
  filter_na <- is.na(unname(st_coordinates(lp_points)[, 1]))
  na_points <- filter(lp_points, filter_na)
  lp_points <- filter(lp_points, !filter_na)
  
  na_points <- st_set_geometry(na_points, NULL)
  
  na_points <- distinct(na_points)
  
  both <- filter(na_points, .data$hu12 %in% lp_points$hu12) # Only broken border HUs included.
  
  na_points <- filter(na_points, !.data$hu12 %in% both)
  lp_points <- filter(lp_points, !.data$hu12 %in% both)
  
  return(list(lp = lp_points, na = na_points))
}

#' Get NA outlet coords
#' @description For HUs that have a levelpath but no intersection 
#' returns the outlet of the interior flowline.
#' @param na_points data.frame as returned in the "na" named return from 
#' get_lp_points.
#' @param net sf data.frame NHDPlus Flowlines
#' @return sf data.frame with outlet locations
#' @export
get_na_outlets_coords <- function(na_points, net) {
  if("NHDFlowline" %in% names(net)) net <- net$NHDFlowline
  
  if("na" %in% names(na_points)) na_points <- na_points$na
  
  if(nrow(na_points) == 0) return(data.frame())
  
  na_outlets <- net %>%
    filter(.data$LevelPathI %in% na_points$lp) %>%
    group_by(.data$LevelPathI) %>%
    filter(.data$Hydroseq == min(.data$Hydroseq, na.rm = TRUE)) %>%
    ungroup() %>%
    left_join(na_points, by = c("LevelPathI" = "lp"))
  
  problem_na <- filter(na_outlets, .data$FromMeas != 0)
  
  na_outlets <- filter(na_outlets, .data$FromMeas == 0)
  
  na_outlet_coords <- st_coordinates(na_outlets) %>%
    as.data.frame()
  
  if(!"L2" %in% names(na_outlet_coords)) {
    na_outlet_coords <- rename(na_outlet_coords, g = .data$L1)
  } else {
    na_outlet_coords <- rename(na_outlet_coords, g = .data$L2) %>%
      select(-.data$L1)
  }
  
  na_outlet_coords <- na_outlet_coords %>%
    group_by(.data$g) %>%
    filter(row_number() == n()) %>%
    ungroup() %>%
    select(-.data$g) %>%
    bind_cols(st_drop_geometry(na_outlets)) %>%
    st_as_sf(coords = c("X", "Y"), crs = st_crs(na_outlets)) %>%
    rename(geom = .data$geometry)
  
  na_outlet_coords$REACH_meas <- 0
  na_outlet_coords$offset <- 0
  
  select(na_outlet_coords, 
         .data$COMID, .data$REACHCODE, .data$REACH_meas, 
         .data$offset, HUC12 = .data$hu12, .data$LevelPathI)
}

#' @noRd
in_list_fun <- function(lp_search, net, lp_points) {
  list(lp_search = lp_search,
       lp_geom = filter(net, .data$LevelPathI == lp_search),
       hu_points = filter(lp_points, .data$lp == lp_search))
}

#' Scalable implementation of get linked points
#' @description Given matched levelpaths, finds outlets of catchments along level paths. 
#' Implements a parallel implementation of get_linked_points
#' @param in_list as returned by get_in_list
#' @param na_outlet_coords as returned in the "na" named output of get_lp_points
#' @param cores numeric (optional) number of cores to use in parallel evaluation
#' @param check_file geopackage file to check for "linked_points" layer -- will load layer and return if found.
#' @export
get_linked_points_scalable <- function(in_list, na_outlet_coords, cores = NA, check_file = NULL) {
  na_outlet_coords <- na_outlet_coords
  
  if(!is.null(check_file) && file.exists(check_file) && 
     "linked_points" %in% sf::st_layers(check_file)$name) {
    linked <- read_sf(check_file, "linked_points")
  } else {
    
    if(is.na(cores)) {
      linked <- lapply(in_list, par_linker)
    } else {
      cl <- parallel::makeCluster(rep("localhost", cores),
                                  type = "SOCK", outfile = "par.log")
      
      linked <- snow::parLapply(cl, in_list, par_linker)
      
      parallel::stopCluster(cl)
    }
    
    linked <- st_sf(do.call(rbind, linked), crs = st_crs(in_list[[1]]$lp_geom)) %>%
      select(.data$COMID, .data$REACHCODE, .data$REACH_meas, 
             .data$offset, HUC12 = .data$hu12, LevelPathI = .data$lp)
    
    if(nrow(na_outlet_coords) > 0) {
      names(linked)[names(linked) == attr(linked, "sf_column")] <-
        names(na_outlet_coords)[names(na_outlet_coords) == attr(na_outlet_coords, "sf_column")]
      
      attr(linked, "sf_column") <- attr(na_outlet_coords, "sf_column")
      
      linked <- rbind(linked, na_outlet_coords) %>%
        st_sf()
    }
  }
  return(linked)
}

#' get in_list for get_linked_points_scalable
#' @description Allows atomic access to in_list creation logic for calling functions used
#' in get_linked_points in a workflow.
#' @param lp_points as returned in the "lp" named output of get_lp_points
#' @param net sf data.frame NHDPlus Flowlines
#' @export
get_in_list <- function(lp_points, net) {
  
  if(names(net) == "NHDFlowline") net <- net$NHDFlowline
  
  lp_points <- lp_points$lp
  
  lp_list <- unique(lp_points$lp)
  
  net <- select(net, .data$COMID, .data$LevelPathI, .data$REACHCODE, 
                .data$ToMeas, .data$FromMeas, .data$Hydroseq) %>%
    filter(.data$LevelPathI %in% lp_list)
  
  pbapply::pblapply(lp_list, in_list_fun, net = net, lp_points = lp_points)
}

#' Parallel hu intersection point linking
#' @noRd
#' @importFrom dplyr group_size row_number
par_linker <- function(lp_list) {
  linked <- NULL
  tryCatch({
    linked <- nhdplusTools::get_flowline_index(lp_list$lp_geom, lp_list$hu_points, search_radius = 1000) %>%
      dplyr::bind_cols(lp_list$hu_points) %>%
      dplyr::left_join(select(sf::st_drop_geometry(lp_list$lp_geom), 
                              .data$COMID, .data$Hydroseq), by = "COMID") %>%
      dplyr::group_by(.data$hu12) %>%
      dplyr::filter(.data$Hydroseq == min(.data$Hydroseq, na.rm = TRUE))
    
    if(any(group_size(linked) > 1)) {
      linked <- linked %>%
        dplyr::group_by(.data$hu12, .data$REACHCODE) %>%
        dplyr::filter(.data$REACH_meas == min(.data$REACH_meas))
    }
    
    linked <- dplyr::ungroup(linked)
  },
  error = function(e) warning(paste(lp_list$lp_search, e)),
  warning = function(w) warning(paste(lp_list$lp_search, w)))
  return(linked)
}

#' HUC points function
#' @description Brute force create points for HU12 intersections 
#' @noRd
hu_points_fun <- function(hp) {
  if(length(unlist(hp)) > 0) {
    out <- do.call(rbind,
                   lapply(1:length(hp),
                          function(x) {
                            if(length(unlist(hp[[x]])) > 0) {
                              o <- lapply(1:length(hp[[x]]),
                                          function(y) st_cast(hp[[x]][y], "POINT"))
                              o <- data.frame(geometry = do.call(c, o))
                              o[["hu12"]] <- names(hp[x])
                            } else {
                              o <- data.frame(hu12 = names(hp[x]), stringsAsFactors = FALSE)
                              o[["geometry"]] <- list(st_point())
                            }
                            return(o)
                          }))
  } else {
    out <- data.frame(hu12 = names(hp), stringsAsFactors = FALSE)
    out[["geometry"]] <- list(st_point())
  }
  
  return(out)
}

#' Run LevelPath get_points_out 
#' @description Find intersection points between HU polygons and level paths
#' @noRd
run_lp <- function(lp_id, net, hu_lp, wbd) {
  out <- NULL
  tryCatch({
    lp <- net[net$LevelPathI == lp_id, ]
    
    outlet <- filter(lp, Hydroseq == LevelPathI) %>%
      nhdplusTools::get_node(position = "end")
    
    lp <- st_geometry(lp)
    
    hu_ids <- hu_lp[hu_lp$corrected_LevelPathI == lp_id, ]
    hus <- wbd[wbd$HUC12 %in% hu_ids$HUC12, ]
    
    outlet_hu <- hus$HUC12[st_contains(hus, outlet, sparse = FALSE)]
    
    st_geometry(hus) <- st_cast(st_geometry(hus), "MULTILINESTRING")
    
    out <- setNames(lapply(1:nrow(hus), 
                           function(i, lp, hus) {
                             suppressMessages(
                               st_intersection(lp, st_geometry(hus)[i])
                             )
                           }, 
                           lp = lp, hus = hus), 
                    hus$HUC12)
    
    if(length(outlet_hu) > 0) {
      out[outlet_hu][[length(out[outlet_hu])]] <- c(out[outlet_hu][[length(out[outlet_hu])]], 
                                                    st_geometry(outlet))
    }
  }, 
  error = function(e) {
    warning(paste(lp_id, e))
  }, 
  warning = function(w) {
    warning(paste(lp_id, w))
  })
  return(out)
}

#' get all points
#' @description deduplicates everything and executes run_lp for all levelpaths
#' @importFrom nhdplusTools get_node
#' @noRd
get_points_out <- function(hu_lp, net, wbd, exclude) {
  
  # Three cases where there are duplicates.
  # one or more terminates in HU one or more don't
  # all terminate in HU
  # None terminate in HU
  outlets <- net %>%
    filter(.data$Hydroseq == .data$LevelPathI) %>% # only level path outlet flowlines
    select(.data$LevelPathI)
  
  outlets <- st_sf(st_drop_geometry(outlets), 
                   geom = st_geometry(get_node(outlets, position = "end")))
  
  outlets <- st_join(outlets, select(wbd, .data$HUC12))
  
  outlets <- st_drop_geometry(rename(outlets, point_outlet = HUC12))
  
  hu_lp <- left_join(hu_lp, outlets, by = c("corrected_LevelPathI" = "LevelPathI"))
  
  terminates_in_hu <- filter(hu_lp, .data$head_HUC12 == .data$point_outlet)
  exits_head_hu <- filter(hu_lp, .data$head_HUC12 != .data$point_outlet)

  term_and_exits <- group_by(hu_lp, .data$HUC12) %>%
    filter(dplyr::n() > 1 &&
             any(.data$HUC12 %in% terminates_in_hu$HUC12) && 
             any(.data$HUC12 %in% exits_head_hu$HUC12))
  
  all_in_or_out <- group_by(hu_lp, .data$HUC12) %>%
    filter(!.data$HUC12 %in% term_and_exits$HUC12 &
             .data$corrected_LevelPathI == min(.data$corrected_LevelPathI)) %>%
    ungroup()
           
  if(nrow(term_and_exits) > 0) {
    term_and_exits <- term_and_exits %>%
      filter(.data$head_HUC12 != .data$point_outlet) %>%
      filter(.data$corrected_LevelPathI == min(.data$corrected_LevelPathI)) %>%
      ungroup()
  }
  
  hu_lp <- bind_rows(all_in_or_out, term_and_exits) %>%
    select(-.data$point_outlet)
  
  wbd <- filter(wbd, !.data$HUC12 %in% exclude)
  hu_lp <- filter(hu_lp, !.data$HUC12 %in% exclude)
  
  lp_ids <- unique(hu_lp$corrected_LevelPathI)
  
  if(st_crs(net) != st_crs(wbd)) {
    net <- st_transform(net, st_crs(wbd))
  }
  
  points <- setNames(pbapply::pblapply(X = lp_ids,
                            FUN = run_lp,
                            net = net, hu_lp = hu_lp, wbd = wbd),
                     lp_ids)
  
  return(points)
}

#' Get level path HUC outlet points
#' @noRd
get_lp_hu_points <- function(points, prj) {
  suppressWarnings({
    lp_points <- lapply(names(points),
                        function(lp, points) {
                          hu_points <- try(bind_rows(lapply(points[lp], hu_points_fun)))
                          
                          if(!is(hu_points, "data.frame")) browser()
                          
                          hu_points[["lp"]] <- lp
                          
                          return(hu_points)
                        }, points = points) 
    lp_points <- st_sf(do.call(rbind, lp_points))
  })
  
  st_crs(lp_points) <- prj
  
  return(lp_points)
}

