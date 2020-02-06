#' Match Level Paths
#' @description Attempts to match dendritic hydrologic unit networks to level paths.
#' @param fline_hu sf data.frame flowlines intersected with hydrologic units containing
#' COMID, Hydroseq, LevelPathI, DnLevelPat, denTotalAreaSqKM, HUC12, TOHUC attributes.
#' @param start_comid integer COMID to start search from.
#' @param add_checks boolean if TRUE, checks for toHUC errors are added.
#' @details Match level paths compares the set of hydrologic units found through
#' spatial intersection to the path identified by navigating the TOHUC codes from
#' headwater to outlet and attempts of reconcile anomolies found.
#'
#' This function is preliminary and subject to revision. It has been tested thoroughly
#' but complete methods description have not yet been published.
#'
#' @export
#' @importFrom sf st_cast st_union st_geometry st_sfc st_sf st_crs st_set_geometry st_line_merge st_geometry_type
#' @importFrom dplyr filter mutate left_join select distinct case_when rename arrange ungroup
#' @importFrom tidyr unnest
#'
#' @examples
#' net_prep <- readRDS(system.file("extdata/match_levelpaths_2279159.rds", package = "hyRefactor"))
#' match_levelpaths(net_prep, 2279159, add_checks = TRUE)
#'
match_levelpaths <- function(fline_hu, start_comid, add_checks = FALSE) {

  #############################################################################
  # get_lp_hu gets levelpath / hydrologic unit pairs from the intersection set.
  #############################################################################
  hu <- get_lp_hu(fline_hu, start_comid)

  if(nrow(hu) == 0) {
    warning(paste("No matches found for", start_comid, "\n"))
    return(NULL)
  }

  #############################################################################
  # create hu with head_hu for each levelpath. label levelpath as "intersected"
  #############################################################################
  hu <- hu %>%
    right_join(distinct(select(fline_hu, HUC12, TOHUC)), by = "HUC12") %>%
    left_join(get_head_hu(hu, fline_hu), by = "LevelPathI") %>%
    mutate(LevelPathI = as.numeric(LevelPathI)) %>%
    rename(intersected_LevelPathI = LevelPathI)

  # In case the outlet flowline spans multiple HUs near the outlet.
  outlet_hu <- sort(fline_hu[which(fline_hu$COMID == start_comid),]$TOHUC,
                     decreasing = FALSE)[1]
  
  if(is.na(outlet_hu)) { # This is sketchy but only for the old WBD.
    outlet_hu <- filter(fline_hu, !is.na(HUC12) & HUC12 != "UNKNOWN" & TOHUC != "UNKNOWN") %>%
      filter(Hydroseq == min(Hydroseq, na.rm = TRUE))
    outlet_hu <- sort(outlet_hu$TOHUC, 
                      decreasing = TRUE)[1]
  }


  #################################################################
  # trace from head_hu to outlet, match paths and reconcile issues.
  #################################################################
  hu <- trace_hu_network(hu, outlet_hu, fline_hu)
  funky_headwaters <- hu$funky_headwaters
  hu <- hu$hu

  #############################################################################
  # Need to fix head_HUC12 where mainstem outlets were changed in trace_hu_network.
  # These are the ones that are part of a tributary so main_LeveLPath got updated above.
  # main_LevelPath == 0 is a different case.
  #############################################################################
  correct_heads <- hu %>%
    filter(intersected_LevelPathI == main_LevelPathI & main_LevelPathI != 0) %>%
    select(main_LevelPathI, correct_head_HUC12 = head_HUC12) %>%
    distinct()

  hu <- left_join(hu, correct_heads, by = "main_LevelPathI") %>%
    mutate(head_HUC12 = ifelse(intersected_LevelPathI != main_LevelPathI & main_LevelPathI != 0,
                               correct_head_HUC12,
                               head_HUC12)) %>%
    select(-correct_head_HUC12) %>%
    distinct()

  #############################################################################
  # Odd situation where a single headwater catchment intersects two HUs causes duplicate HUC12s
  #############################################################################
  hu <- group_by(hu, HUC12) %>%
    arrange(head_HUC12) %>%
    filter(dplyr::row_number() == 1) %>%
    ungroup()

  #############################################################################
  # main_levelpath was found by tracing toHUC downstream.
  # intersected_LevelPath was passed in from spatial intersection
  # Not all "main_levelpaths" got filled in.
  # The goal of the following code is to get a "corrected_levelpath"
  #############################################################################

  #############################################################################
  # For outlets, just set main_levelpath to intersected.
  #############################################################################
  hu <- mutate(hu, main_LevelPathI = ifelse(HUC12 == outlet_hu, intersected_LevelPathI, main_LevelPathI))

  #############################################################################
  # Run corrections on results
  #############################################################################
  return(correct_hu(hu, fline_hu, funky_headwaters, add_checks))
}

gather_ds <- function(hu_data, huc12) {
  next_huc12 <- unique(hu_data[["TOHUC"]][hu_data[["HUC12"]] %in% huc12])
  if(length(next_huc12) == 0) {
    return()
  }
  return(c(huc12, gather_ds(hu_data, next_huc12)))
}

get_head_hu <- function(lp_hu, fline_hu) {
  lp_hu %>%
    left_join(select(fline_hu, Hydroseq,
                     nhd_LevelPath = LevelPathI, HUC12), by = "HUC12") %>%
    filter(LevelPathI == nhd_LevelPath) %>%
    group_by(LevelPathI) %>%
    filter(Hydroseq == max(Hydroseq)) %>%
    ungroup() %>%
    select(-Hydroseq, -nhd_LevelPath, head_HUC12 = HUC12)
}

get_lp_hu <- function(fline_hu, start_comid) {

  nlp <- unique(filter(fline_hu, COMID == start_comid)[["LevelPathI"]])

  check <- TRUE
  lp_hu <- list() # List to store levelpath/HU pairs
  next_lp <- c() # Vector for sets of levelpaths that are at the next level from current.
  nlp_tracker <- c() # Tracker for levelpaths that need to be descended into later.
  count <- 0 # Stop checker for while loop.
  none_count <- 0 # performance improvement to not check too much stuff.
  keep_going <- FALSE # Solves an edge case with small paths downstream of large paths.

  # Nothing to do.
  if(all(is.na(unique(fline_hu$HUC12[fline_hu$LevelPathI == nlp])))) check <- FALSE

  # There's a chance that this search could be done with an artfully crafted
  # grouped filter but I've not been able to wrap my head around getting it
  # right in all cases.
  while(check == TRUE & count < 100000) {
    # get the HUC12s that intersect the nlp we are looking for.
    lp_hu_temp <- unique(fline_hu$HUC12[which(fline_hu$LevelPathI == nlp)])

    if(length(lp_hu_temp) == 1) if(is.na(lp_hu_temp)) lp_hu_temp <- character(0)

    if(length(lp_hu_temp) > 0) { # if hu12s are found
      lp_hu[as.character(nlp)] <- list(lp_hu_temp) # save that list.

      # filter the found hu12s out of the set we search next time.
      # Remove all that match the HU12s and current nlp.
      remove <- fline_hu$LevelPathI == nlp & fline_hu$HUC12 %in% lp_hu[[as.character(nlp)]]
      
      fline_hu <- filter(fline_hu, !remove)
      
      nlp_tracker <- c(nlp_tracker, nlp) # record this one to zoom in on later.

    } else {
      if(any(fline_hu$DnLevelPat == nlp, na.rm = TRUE)) {
        nlp_tracker <- c(nlp, nlp_tracker)
        keep_going <- TRUE
      }
    }

    if(length(next_lp) == 0 | keep_going) { # If on the last nlp next_lp will be empty.
      
      keep_going <- FALSE
      
      i <- 0
      while(length(next_lp) == 0 & length(nlp_tracker) > 0) {
        
        next_lp <- get_next_lp(fline_hu, nlp_tracker)

        next_lp <- next_lp[["LevelPathI"]]
        
        i <- i + 1
        if(i > 1000) stop("runaway loop?")

      if(length(nlp_tracker) > 1) { # maintain backlog that needs to be worked through.
        nlp_tracker <- nlp_tracker[2:length(nlp_tracker)]
      } else {
        nlp_tracker <- c()
      }
      }
    }

    nlp <- next_lp[1]

    if(length(next_lp) > 1) {
      next_lp <- next_lp[2:length(next_lp)]
    } else {
      next_lp <- c()
    }

    if(length(next_lp) == 0 & length(nlp_tracker) == 0 & (is.null(nlp) || is.na(nlp))) check <- FALSE

    count <- count + 1
  }

  return(data.frame(LevelPathI = names(lp_hu),
                    HUC12 = I(lp_hu),
                    stringsAsFactors = FALSE) %>%
           tidyr::unnest(cols = c(HUC12)))
}

get_next_lp <- function(fline_hu, nlp_tracker) {
  # Grab all the levelpaths that intersect the one we are on.
  filter(fline_hu,
                    fline_hu$DnLevelPat == nlp_tracker[1] &
                      !LevelPathI == nlp_tracker[1]) %>%
    group_by(LevelPathI) %>%
    # Pick only the outlet catchment/flowline
    filter(Hydroseq == min(Hydroseq)) %>%
    # Sort from biggest drainage area to smallest.
    arrange(desc(denTotalAreaSqKM))
}

trace_hu_network <- function(hu, outlet_hu, fline_hu) {
  # Make sure we don't have the outlet in scope so the network breaks where we want.
  hu_destructive <- filter(hu, !HUC12 %in% outlet_hu)

  hu[["main_LevelPathI"]] <- 0
  hu[["outlet_HUC12"]] <- ""

  lps <- sort(unique(hu$intersected_LevelPathI))
  funky_headwaters <- c()
  broken_path <- c()

  count <- 0
  for(lp in lps) {
    lp <- filter(hu, intersected_LevelPathI == lp)
    # could gather_ds levelpath instead of hu but some toHUCs go outside the levelpath intersection!
    main_stem <- tryCatch(gather_ds(hu_destructive, lp$head_HUC12[1]),
                          error = function(e) {
                            warning(paste("something real bad happened with levelpath",
                                          lp$intersected_LevelPathI[1]))
                            funky_headwaters <- c(funky_headwaters, lp$head_HUC12[1])})

    # If any of the HUs found are still available to be allocated...
    if(any(main_stem %in% hu_destructive$HUC12)) {
      # If the main stem found doesn't follow the expected path at all.
      # Then something is wrong with the head_HUC12 and we need to try a different one.
      candidates <- main_stem[!main_stem == lp$head_HUC12[1]]
      if(!any(candidates %in% lp$HUC12) &
         length(main_stem) > 1) {

        checker <- filter(fline_hu, HUC12 %in% candidates & LevelPathI == lp$intersected_LevelPathI[1])
        if(!any(candidates %in% checker$HUC12)) {
          # Need to find the actual top HU12 first grab all that intersect this LP.
          top_cat <- filter(fline_hu, LevelPathI == lp$intersected_LevelPathI[1]) %>%
            # In case there is a catchment completely outside the HU
            # We have to find the one that overlaps the boundary!!
            group_by(COMID) %>% filter(n() > 1) %>% ungroup() %>%
            # Now grab the top most row that isn't the one we found above.
            filter(Hydroseq == max(Hydroseq) & HUC12 != lp$head_HUC12[1])

          # No test for this, but it shouldn't happen so throw error!!
          if(nrow(top_cat) > 1) {
            warning(paste("something very wrong with headwaters around HUC",
                          lp$head_HUC12[1], "and levelpath", lp$intersected_LevelPathI[1]), "\n")
            funky_headwaters <- c(funky_headwaters, top_cat$HUC12)
            main_stem <- c()
          } else if(nrow(top_cat) == 0) {
            warning(paste("broken path along levelpath", lp$intersected_LevelPathI, "passing by."))
            broken_path <- c(broken_path, lp$head_HUC12)
          } else {
            # set and rerun.
            lp$head_HUC12 <- top_cat$HUC12
            funky_headwaters <- c(funky_headwaters, top_cat$HUC12)
            main_stem <- gather_ds(hu_destructive, lp$head_HUC12[1])
          }
        }
      }

      hu <- mutate(hu, main_LevelPathI = ifelse(HUC12 %in% main_stem &
                                                  main_LevelPathI == 0,
                                                lp$intersected_LevelPathI[1],
                                                main_LevelPathI),
                   outlet_HUC12 = ifelse(HUC12 %in% main_stem,
                                         main_stem[length(main_stem)],
                                         outlet_HUC12))
      hu_destructive <- filter(hu_destructive, !HUC12 %in% main_stem)
    }
    count <- count + 1
    if (count %% 100 == 0) {
      message(paste(count, "of", length(lps), "levelpaths"))
    }
  }
  return(list(hu = hu, funky_headwaters = funky_headwaters, broken_path = broken_path))
}

correct_hu <- function(hu, fline_hu, funky_headwaters, add_checks) {
  ################################################################################
  # HUs found through intersection but not through main path trace belong with trib
  ################################################################################
  hu_trib <- hu %>%
    left_join(distinct(select(fline_hu, HUC12, check_LevelPathI = LevelPathI)), by = "HUC12") %>%
    # only modify ones not found to be on the main path in the loop above
    # and where the LevelPathI assigned is not equal to the original one assigned
    filter(main_LevelPathI == 0 & intersected_LevelPathI != check_LevelPathI) %>%
    group_by(HUC12) %>%
    # This grabs the biggest trib in the HU after filtering out the originally assigned one.
    filter(check_LevelPathI == min(check_LevelPathI))

  hu <- hu %>%
    left_join(select(hu_trib, HUC12, check_LevelPathI), by = "HUC12") %>%
    # Kill head_HUC12 which is now wrong.
    mutate(head_HUC12 = ifelse(main_LevelPathI == 0, "", head_HUC12),
           # The change to main_LevelPath makes this one an outlet.
           outlet_HUC12 = ifelse(main_LevelPathI == 0 & !is.na(check_LevelPathI), HUC12, outlet_HUC12),
           # Update the main_levelpath where needed per hu_trib from above.
           main_LevelPathI = as.numeric(ifelse(main_LevelPathI == 0, check_LevelPathI, main_LevelPathI)),
           main_LevelPathI = ifelse(is.na(main_LevelPathI), intersected_LevelPathI, main_LevelPathI)) %>%
    select(-check_LevelPathI)

  ################################################################################
  # HUs found along trace but not along intersection need to be flagged and checked
  ################################################################################
  hu_trib2 <- left_join(hu, select(fline_hu, HUC12, check_LevelPathI = LevelPathI), by = "HUC12") %>%
    group_by(HUC12) %>%
    # When nothing in the group was found intersecting the levelpath
    # it should actually me on another level path.
    filter(!any(main_LevelPathI == check_LevelPathI) & check_LevelPathI == min(check_LevelPathI)) %>%
    distinct() %>%
    ungroup()

  hu <- hu %>%
    left_join(select(hu_trib2, HUC12, check_LevelPathI), by = "HUC12") %>%
    mutate(main_LevelPathI = as.numeric(ifelse(!is.na(check_LevelPathI), check_LevelPathI, main_LevelPathI))) %>%
    select(HUC12, TOHUC, intersected_LevelPathI, corrected_LevelPathI = main_LevelPathI, head_HUC12, outlet_HUC12)

  hu <- filter(hu, !is.na(intersected_LevelPathI))

  ################################################################################
  # Update head_hu where they were broken before
  ################################################################################
  head_hu <- hu %>%
    filter(head_HUC12 == "") %>%
    select(LevelPathI = corrected_LevelPathI, HUC12) %>%
    get_head_hu(fline_hu) %>%
    select(LevelPathI, update_head_HUC12 = head_HUC12)

  hu <- hu %>%
    left_join(head_hu, by = c("corrected_LevelPathI" = "LevelPathI")) %>%
    mutate(head_HUC12 = ifelse(head_HUC12 == "", update_head_HUC12, head_HUC12)) %>%
    select(-update_head_HUC12)

  # need to deduplicate resulting head_HUC12 in some edge cases
  lp_head <- select(hu, corrected_LevelPathI, head_HUC12) %>%
    distinct() %>%
    group_by(corrected_LevelPathI) %>%
    arrange(head_HUC12) %>%
    filter(dplyr::row_number() == 1) %>%
    ungroup()

  hu <- hu %>%
    left_join(select(lp_head, corrected_LevelPathI, updated_head_HUC12 = head_HUC12),
              by = "corrected_LevelPathI") %>%
    select(-head_HUC12) %>% rename(head_HUC12 = updated_head_HUC12)

  ################################################################################
  # When HUC12 doesn't have anything going to it -- force it to be a head_HUC12
  ################################################################################
  broken_head <- !hu$HUC12 %in% hu$TOHUC & !hu$head_HUC12 == hu$HUC12
  hu <- mutate(hu, head_HUC12 = ifelse(broken_head, HUC12, head_HUC12),
               corrected_LevelPathI = ifelse(broken_head, -1, corrected_LevelPathI)) # sketchy but should work?

  ################################################################################
  # Fix up outlets that got broken in edge cases above.
  ################################################################################
  if(any(hu$outlet_HUC12 == "")) {
    lp_outlet <- select(hu, corrected_LevelPathI, outlet_HUC12) %>%
      filter(outlet_HUC12 != "") %>%
      distinct()

    hu <- hu %>%
      left_join(select(lp_outlet, corrected_LevelPathI, updated_outlet_HUC12 = outlet_HUC12),
                by = "corrected_LevelPathI") %>%
      select(-outlet_HUC12) %>% rename(outlet_HUC12 = updated_outlet_HUC12) %>%
      filter(outlet_HUC12 != "" | corrected_LevelPathI == -1) %>%
      mutate(outlet_HUC12 = ifelse(corrected_LevelPathI == -1, HUC12, outlet_HUC12))
  }

  lp_outlet <- select(hu, corrected_LevelPathI, outlet_HUC12) %>%
    distinct() %>%
    group_by(corrected_LevelPathI) %>%
    arrange(dplyr::desc(outlet_HUC12)) %>%
    filter(dplyr::row_number() == 1) %>%
    ungroup()

  hu <- hu %>%
    left_join(select(lp_outlet, corrected_LevelPathI, updated_outlet_HUC12 = outlet_HUC12),
              by = "corrected_LevelPathI") %>%
    select(-outlet_HUC12) %>% 
    rename(outlet_HUC12 = updated_outlet_HUC12) %>%
    filter(corrected_LevelPathI != -1) %>%
    distinct()

  ################################################################################
  # Add checks
  ################################################################################
  if(add_checks) {
    hu <- mutate(hu, trib_intersect = HUC12 %in% hu_trib$HUC12,
                 trib_no_intersect = HUC12 %in% hu_trib2$HUC12,
                 headwater_error = HUC12 %in% funky_headwaters)
  }
  return(hu)
}

par_match_levelpaths_fun <- function(start_comid, net_atts, net_prep, wbd_atts, temp_dir) {
  library(nhdplusTools)
  library(sf)
  library(dplyr)
  
  out_file <- file.path(temp_dir, paste0(start_comid, ".rds"))
  
  if(!file.exists(out_file)) {

    all_comid <- get_UT(net_atts, start_comid)
    
    sub_net <- filter(net_prep, COMID %in% all_comid)
    
    out <- list(NULL)
    
    if(nrow(sub_net) > 0) {
      sub_net <- bind_rows(sub_net,
                           select(filter(wbd_atts,
                                         !HUC12 %in% sub_net$HUC12),
                                  HUC12, TOHUC))
      
      out <- list(match_levelpaths(sub_net, start_comid, add_checks = TRUE))
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
#' @param temp_dir directory to write temprary files to
#' @param out_dir directory to check for cached output.
#' @export
par_match_levelpaths <- function(net, wbd, simp, cores, temp_dir = "temp/", 
                                 out_file = "temp.csv", net_int = NULL, purge_temp = TRUE) {
  
  if(names(net) == "NHDFlowline") net <- net$NHDFlowline
  
  if(file.exists(out_file)) {
    all <- readr::read_csv(out_file)
  } else {
    
    if(purge_temp) unlink(temp_dir, recursive = TRUE)
    
    if(is.null(net_int)) {
      net_int <- get_process_data(net, wbd, simp)
    }
    
    net <- st_set_geometry(net, NULL)
    wbd <- st_set_geometry(wbd, NULL)
    
    terminals <- net %>%
      select(TerminalPa) %>%
      distinct() %>%
      left_join(select(net, COMID, LevelPathI, Hydroseq), 
                by = c("TerminalPa" = "LevelPathI")) %>%
      filter(COMID %in% net_int$COMID) %>%
      group_by(TerminalPa) %>%
      filter(Hydroseq == min(Hydroseq))
    
    if(cores > 1) {
      cl <- parallel::makeCluster(rep("localhost", cores), 
                                  type = "SOCK", 
                                  outfile = "hu_joiner.log")
    } else {
      cl <- NULL
    }
    
    to_run <- terminals$COMID
    
    dir.create(temp_dir, showWarnings = FALSE)
    already_run <- list.files(temp_dir, pattern = "*.rds")
    already_run <- as.numeric(gsub(".rds", "", already_run))
    
    to_run <- to_run[!to_run %in% already_run]
    
    all_outlets <- pblapply(to_run, par_match_levelpaths_fun,
                             net_atts = net,
                             net_prep = net_int,
                             wbd_atts = wbd,
                             temp_dir = temp_dir, cl = cl)
    if(!is.null(cl)) {
      parallel::stopCluster(cl)
    }
    
    out_files <- list.files(temp_dir, full.names = TRUE)
    
    all <- sapply(out_files, readRDS, USE.NAMES = TRUE)
    
    names(all) <- gsub(temp_dir, "", names(all))
    names(all) <- gsub(".rds", "", names(all))
    
    all <- bind_rows(all)
    
    # Once everything is done there are some orphans where only one
    # matching network element was found. We can just add those.
    missed <- net_int %>%
      filter(!HUC12 %in% all$HUC12) %>%
      group_by(HUC12)
    
    add_match <- missed %>%
      filter(n() == 1) %>%
      ungroup()
    
    missed <- missed %>%
      filter(!HUC12 %in% add_match$HUC12)
    
    add_match <- select(add_match, HUC12 = HUC12, TOHUC = TOHUC, 
                        intersected_LevelPathI = LevelPathI)
    
    if(nrow(add_match) > 0) {
      add_match$corrected_LevelPathI <- add_match$intersected_LevelPathI
      add_match$head_HUC12 <- add_match$HUC12
      add_match$outlet_HUC12 <- add_match$HUC12
      add_match$trib_intersect <- FALSE
      add_match$trib_no_intersect <- FALSE
      add_match$headwater_error <- FALSE
    }
    
    all <- bind_rows(all, add_match)
    
    readr::write_csv(all, out_file)
    
  }
  
  return(all)
}

get_process_data <- function(net, wbd, simp) {
  
  if(names(net) == "NHDFlowline") net <- net$NHDFlowline
  
  net_prep <- prep_net(net, simp)
  
  wbd <- select(st_simplify(wbd, dTolerance = simp), HUC12, TOHUC)
  
  net_prep <- st_join(net_prep, wbd) %>%
    st_set_geometry(NULL)
  
  return(net_prep)
}

#' @import nhdplusTools sf dplyr
prep_net <- function(net, simp) {
  
  net_prep <- prepare_nhdplus(net, 
                              min_network_size = 0, # sqkm
                              min_path_length = 0, # sqkm
                              min_path_size = 0, # sqkm
                              purge_non_dendritic = TRUE,
                              warn =  TRUE, error = FALSE) %>%
    left_join(select(net, COMID, DnLevelPat, AreaSqKM), by = "COMID") %>%
    st_sf() %>%
    group_by(LevelPathI) %>%
    arrange(Hydroseq) %>%
    mutate(DnLevelPat = DnLevelPat[1]) %>%
    ungroup()
  
  net_prep["denTotalAreaSqKM"] <-
    nhdplusTools::calculate_total_drainage_area(select(st_set_geometry(net_prep, NULL),
                                         ID = COMID, toID = toCOMID,
                                         area = AreaSqKM))
  
  net_prep <- st_simplify(net_prep, dTolerance = simp)
  
  return(net_prep)
}
