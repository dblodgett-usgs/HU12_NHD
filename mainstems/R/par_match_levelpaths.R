par_match_levelpaths_fun <- function(start_comid, net_atts, net_prep, wbd_atts, 
                                     temp_dir, cluster = NULL) {
  
  out_file <- file.path(temp_dir, paste0(start_comid, ".rds"))
  
  if(!file.exists(out_file)) {
    
    all_comid <- nhdplusTools::get_UT(net_atts, start_comid)
    
    sub_net <- dplyr::filter(net_prep, .data$COMID %in% all_comid)
    
    out <- list(NULL)
    
    if(nrow(sub_net) > 0) {
      sub_net <- dplyr::bind_rows(sub_net,
                                  dplyr::select(filter(wbd_atts,
                                                       !.data$HUC12 %in% sub_net$HUC12),
                                                .data$HUC12, .data$TOHUC))
      
      out <- list(mainstems::match_levelpaths(sub_net, start_comid, 
                                              add_checks = TRUE, cl = cluster))
      
    }
    
    saveRDS(out, out_file)
    
  }
  return(out_file)
}

#' Parallel match levelpaths
#' @details Calls match_levelpaths in a parallel mode for use in large process workflows
#' @param net NHDPlus network
#' @param wbd WBD HU12 polygons
#' @param simp simplification in units of WBD and NHDPlus to limit memory usage
#' @param cores Number of parallel processes to run. (1 for now process splitting)
#' @param temp_dir directory to write temprary files to
#' @param out_file file to write results to.
#' @param net_int intermediate process artifact to be used in drake plans output 
#' of get_process_data
#' @param purge_temp should temp directory be purged?
#' @export
par_match_levelpaths <- function(net, wbd, simp, cores, temp_dir = "temp/", 
                                 out_file = "temp.csv", net_int = NULL, 
                                 purge_temp = TRUE) {
  
  if(length(names(net)) == 1 && names(net) == "NHDFlowline") net <- net$NHDFlowline
  
  if(file.exists(out_file)) {
    all <- readr::read_csv(out_file)
  } else {
    
    if(purge_temp) unlink(file.path(temp_dir, "*"), recursive = TRUE)
    
    if(is.null(net_int)) {
      net_int <- get_process_data(net, wbd, simp)
    }
    
    net <- st_set_geometry(net, NULL)
    wbd <- st_set_geometry(wbd, NULL)
    
    terminals <- net %>%
      select(.data$TerminalPa) %>%
      distinct() %>%
      left_join(select(net, .data$COMID, .data$LevelPathI, .data$Hydroseq), 
                by = c("TerminalPa" = "LevelPathI")) %>%
      filter(.data$COMID %in% net_int$COMID) %>%
      group_by(.data$TerminalPa) %>%
      filter(.data$Hydroseq == min(.data$Hydroseq, na.rm= TRUE)) %>%
      ungroup()
    
    big_terminals <- terminals %>%
      left_join(select(net, .data$COMID, .data$TotDASqKM), by = "COMID") %>%
      filter(.data$TotDASqKM > 100000) %>%
      select(-.data$TotDASqKM)
    
    terminals <- filter(terminals, !.data$COMID %in% big_terminals$COMID)
    
    if(cores > 1) {
      cl <- parallel::makeCluster(cores, outfile = "hu_joiner.log")
    } else {
      cl <- NULL
    }
    
    to_run <- terminals$COMID
    to_run_big <- big_terminals$COMID
    
    dir.create(temp_dir, showWarnings = FALSE)
    already_run <- list.files(temp_dir, pattern = "*.rds")
    already_run <- as.numeric(gsub(".rds", "", already_run))
    
    to_run <- to_run[!to_run %in% already_run]
    to_run_big <- to_run_big[!to_run_big %in% already_run]
    
    small_outlets <- pblapply(to_run, par_match_levelpaths_fun,
                              net_atts = net,
                              net_prep = net_int,
                              wbd_atts = wbd,
                              temp_dir = temp_dir, cl = cl)
    
    big_outlets <- pblapply(to_run_big, par_match_levelpaths_fun,
                            net_atts = net,
                            net_prep = net_int,
                            wbd_atts = wbd,
                            temp_dir = temp_dir,
                            cluster = cl)
    
    if(!is.null(cl)) {
      parallel::stopCluster(cl)
    }
    
    out_files <- list.files(temp_dir, full.names = TRUE)
    
    all <- sapply(out_files, readRDS, USE.NAMES = TRUE)
    
    names(all) <- basename(names(all))
    names(all) <- gsub(".rds", "", names(all))
    
    all <- do.call(rbind, all)
    
    # Once everything is done there are some orphans where only one
    # matching network element was found. We can just add those.
    missed <- net_int %>%
      filter(!.data$HUC12 %in% all$HUC12) %>%
      group_by(.data$HUC12)
    
    add_match <- missed %>%
      filter(n() == 1) %>%
      ungroup()
    
    missed <- missed %>%
      filter(!.data$HUC12 %in% add_match$HUC12)
    
    add_match <- select(add_match, HUC12 = .data$HUC12, TOHUC = .data$TOHUC, 
                        intersected_LevelPathI = .data$LevelPathI)
    
    if(nrow(add_match) > 0) {
      add_match$corrected_LevelPathI <- add_match$intersected_LevelPathI
      add_match$head_HUC12 <- add_match$HUC12
      add_match$outlet_HUC12 <- add_match$HUC12
      add_match$trib_intersect <- FALSE
      add_match$trib_no_intersect <- FALSE
      add_match$headwater_error <- FALSE
    }
    
    all <- rbind(all, add_match)
    
    readr::write_csv(all, out_file)
    
  }
  
  return(all)
}

get_process_data <- function(net, wbd, simp) {
  
  if(length(names(net)) == 1 && names(net) == "NHDFlowline") net <- net$NHDFlowline
  
  net_prep <- prep_net(net, simp)
  
  wbd <- select(st_simplify(wbd, dTolerance = simp), .data$HUC12, .data$TOHUC)
  
  net_prep <- st_join(net_prep, wbd) %>%
    st_set_geometry(NULL)
  
  return(net_prep)
}

#' @import nhdplusTools sf dplyr
prep_net <- function(net, simp) {
  
  if(!"StreamOrde" %in% names(net)) {
    net$StreamOrde <- 1
    net$StreamCalc <- 1
  }
  
  net_prep <- rename(net, toCOMID = tocomid) %>%
    st_sf() %>%
    group_by(.data$LevelPathI) %>%
    arrange(.data$Hydroseq) %>%
    mutate(DnLevelPat = .data$DnLevelPat[1]) %>%
    ungroup()
  
  net_prep["denTotalAreaSqKM"] <-
    nhdplusTools::calculate_total_drainage_area(select(st_set_geometry(net_prep, NULL),
                                                       ID = .data$COMID, toID = .data$toCOMID,
                                                       area = .data$AreaSqKM))
  
  net_prep <- st_simplify(net_prep, dTolerance = simp)
  
  return(net_prep)
}