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
    right_join(distinct(select(fline_hu, .data$HUC12, .data$TOHUC)), by = "HUC12") %>%
    left_join(get_head_hu(hu, fline_hu), by = "LevelPathI") %>%
    mutate(LevelPathI = as.numeric(.data$LevelPathI)) %>%
    rename(intersected_LevelPathI = .data$LevelPathI)

  # In case the outlet flowline spans multiple HUs near the outlet.
  outlet_hu <- sort(fline_hu[which(fline_hu$COMID == start_comid),]$TOHUC,
                     decreasing = FALSE)[1]
  
  if(is.na(outlet_hu)) { # This is sketchy but only for the old WBD.
    outlet_hu <- filter(fline_hu, !is.na(.data$HUC12) & .data$HUC12 != "UNKNOWN" & 
                          .data$TOHUC != "UNKNOWN") %>%
      filter(.data$Hydroseq == min(.data$Hydroseq, na.rm = TRUE))
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
    filter(.data$intersected_LevelPathI == .data$main_LevelPathI & .data$main_LevelPathI != 0) %>%
    select(.data$main_LevelPathI, correct_head_HUC12 = .data$head_HUC12) %>%
    distinct()

  hu <- left_join(hu, correct_heads, by = "main_LevelPathI") %>%
    mutate(head_HUC12 = ifelse(.data$intersected_LevelPathI != .data$main_LevelPathI & 
                                 .data$main_LevelPathI != 0,
                               .data$correct_head_HUC12,
                               .data$head_HUC12)) %>%
    select(-.data$correct_head_HUC12) %>%
    distinct()

  #############################################################################
  # Odd situation where a single headwater catchment intersects two HUs causes duplicate HUC12s
  #############################################################################
  hu <- group_by(hu, .data$HUC12) %>%
    arrange(.data$head_HUC12) %>%
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
  hu <- mutate(hu, main_LevelPathI = ifelse(.data$HUC12 == outlet_hu, .data$intersected_LevelPathI, .data$main_LevelPathI))

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
  out <- lp_hu %>%
    left_join(select(fline_hu, .data$Hydroseq,
                     nhd_LevelPath = .data$LevelPathI, .data$HUC12), by = "HUC12") %>%
    filter(.data$LevelPathI == .data$nhd_LevelPath & !is.na(.data$HUC12)) %>%
    group_by(.data$LevelPathI)
  
  if(nrow(out) > 0) {
    filter(out, .data$Hydroseq == max(.data$Hydroseq)) %>%
      ungroup() %>%
      select(-.data$Hydroseq, -.data$nhd_LevelPath, head_HUC12 = .data$HUC12)
  } else {
    select(ungroup(out), -.data$Hydroseq, -.data$nhd_LevelPath, head_HUC12 = .data$HUC12)
  }
}

get_next_lp <- function(fline_hu, nlp_tracker) {
  # Grab all the levelpaths that intersect the one we are on.
  out <- filter(fline_hu,
                .data$DnLevelPat == nlp_tracker[1] &
                  !.data$LevelPathI == nlp_tracker[1]) %>%
    group_by(.data$LevelPathI)
    
  if(nrow(out) > 0) {
    # Pick only the outlet catchment/flowline
    filter(out, .data$Hydroseq == min(.data$Hydroseq, na.rm = TRUE)) %>%
    # Sort from biggest drainage area to smallest.
    arrange(desc(.data$denTotalAreaSqKM))
  } else {
    out
  }
}

trace_hu_network <- function(hu, outlet_hu, fline_hu) {
  # Make sure we don't have the outlet in scope so the network breaks where we want.
  hu_destructive <- filter(hu, !.data$HUC12 %in% outlet_hu)

  hu[["main_LevelPathI"]] <- 0
  hu[["outlet_HUC12"]] <- ""

  lps <- sort(unique(hu$intersected_LevelPathI))
  funky_headwaters <- c()
  broken_path <- c()

  count <- 0
  for(lp in lps) {
    
    lp <- filter(hu, .data$intersected_LevelPathI == lp)
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

        checker <- filter(fline_hu, .data$HUC12 %in% candidates & .data$LevelPathI == lp$intersected_LevelPathI[1])
        if(!any(candidates %in% checker$HUC12)) {
          # Need to find the actual top HU12 first grab all that intersect this LP.
          top_cat <- filter(fline_hu, .data$LevelPathI == lp$intersected_LevelPathI[1]) %>%
            # In case there is a catchment completely outside the HU
            # We have to find the one that overlaps the boundary!!
            group_by(.data$COMID) %>% filter(n() > 1) %>% ungroup() %>%
            # Now grab the top most row that isn't the one we found above.
            filter(.data$Hydroseq == max(.data$Hydroseq) & .data$HUC12 != lp$head_HUC12[1])

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

      hu <- mutate(hu, main_LevelPathI = ifelse(.data$HUC12 %in% main_stem &
                                                  .data$main_LevelPathI == 0,
                                                lp$intersected_LevelPathI[1],
                                                .data$main_LevelPathI),
                   outlet_HUC12 = ifelse(.data$HUC12 %in% main_stem,
                                         main_stem[length(main_stem)],
                                         .data$outlet_HUC12))
      hu_destructive <- filter(hu_destructive, !.data$HUC12 %in% main_stem)
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
    left_join(distinct(select(fline_hu, .data$HUC12, check_LevelPathI = .data$LevelPathI)), 
              by = "HUC12") %>%
    # only modify ones not found to be on the main path in the loop above
    # and where the LevelPathI assigned is not equal to the original one assigned
    filter(.data$main_LevelPathI == 0 & .data$intersected_LevelPathI != .data$check_LevelPathI) %>%
    group_by(.data$HUC12)
  
  if(nrow(hu_trib) > 0) {
    # This grabs the biggest trib in the HU after filtering out the originally assigned one.
    hu_trib <- filter(hu_trib, .data$check_LevelPathI == min(.data$check_LevelPathI))
  }

  hu <- hu %>%
    left_join(select(hu_trib, .data$HUC12, .data$check_LevelPathI), by = "HUC12") %>%
    # Kill head_HUC12 which is now wrong.
    mutate(head_HUC12 = ifelse(.data$main_LevelPathI == 0, "", .data$head_HUC12),
           # The change to main_LevelPath makes this one an outlet.
           outlet_HUC12 = ifelse(.data$main_LevelPathI == 0 & !is.na(.data$check_LevelPathI), 
                                 .data$HUC12, .data$outlet_HUC12),
           # Update the main_levelpath where needed per hu_trib from above.
           main_LevelPathI = as.numeric(ifelse(.data$main_LevelPathI == 0, .data$check_LevelPathI, 
                                               .data$main_LevelPathI)),
           main_LevelPathI = ifelse(is.na(.data$main_LevelPathI), 
                                    .data$intersected_LevelPathI, .data$main_LevelPathI)) %>%
    select(-.data$check_LevelPathI)

  ################################################################################
  # HUs found along trace but not along intersection need to be flagged and checked
  ################################################################################
  hu_trib2 <- left_join(hu, select(fline_hu, .data$HUC12, check_LevelPathI = .data$LevelPathI), 
                        by = "HUC12") %>%
    group_by(.data$HUC12) %>%
    # When nothing in the group was found intersecting the levelpath
    # it should actually me on another level path.
    filter(!any(.data$main_LevelPathI == .data$check_LevelPathI) & 
             .data$check_LevelPathI == min(.data$check_LevelPathI)) %>%
    distinct() %>%
    ungroup()

  hu <- hu %>%
    left_join(select(hu_trib2, .data$HUC12, .data$check_LevelPathI), by = "HUC12") %>%
    mutate(main_LevelPathI = as.numeric(ifelse(!is.na(.data$check_LevelPathI), 
                                               .data$check_LevelPathI, .data$main_LevelPathI))) %>%
    select(.data$HUC12, .data$TOHUC, .data$intersected_LevelPathI, 
           corrected_LevelPathI = .data$main_LevelPathI, .data$head_HUC12, .data$outlet_HUC12)

  hu <- filter(hu, !is.na(.data$intersected_LevelPathI))

  ################################################################################
  # Update head_hu where they were broken before
  ################################################################################
  head_hu <- hu %>%
    filter(.data$head_HUC12 == "") %>%
    select(LevelPathI = .data$corrected_LevelPathI, .data$HUC12) %>%
    get_head_hu(fline_hu) %>%
    select(.data$LevelPathI, update_head_HUC12 = .data$head_HUC12)

  hu <- hu %>%
    left_join(head_hu, by = c("corrected_LevelPathI" = "LevelPathI")) %>%
    mutate(head_HUC12 = ifelse(.data$head_HUC12 == "", 
                               .data$update_head_HUC12, .data$head_HUC12)) %>%
    select(-.data$update_head_HUC12)

  # need to deduplicate resulting head_HUC12 in some edge cases
  lp_head <- select(hu, .data$corrected_LevelPathI, .data$head_HUC12) %>%
    distinct() %>%
    group_by(.data$corrected_LevelPathI) %>%
    arrange(.data$head_HUC12) %>%
    filter(dplyr::row_number() == 1) %>%
    ungroup()

  hu <- hu %>%
    left_join(select(lp_head, .data$corrected_LevelPathI, updated_head_HUC12 = .data$head_HUC12),
              by = "corrected_LevelPathI") %>%
    select(-.data$head_HUC12) %>% rename(head_HUC12 = .data$updated_head_HUC12)

  ################################################################################
  # When HUC12 doesn't have anything going to it -- force it to be a head_HUC12
  ################################################################################
  broken_head <- !hu$HUC12 %in% hu$TOHUC & !hu$head_HUC12 == hu$HUC12
  hu <- mutate(hu, head_HUC12 = ifelse(broken_head, .data$HUC12, .data$head_HUC12),
               corrected_LevelPathI = ifelse(broken_head, -1, .data$corrected_LevelPathI)) # sketchy but should work?

  ################################################################################
  # Fix up outlets that got broken in edge cases above.
  ################################################################################
  if(any(hu$outlet_HUC12 == "")) {
    lp_outlet <- select(hu, .data$corrected_LevelPathI, .data$outlet_HUC12) %>%
      filter(.data$outlet_HUC12 != "") %>%
      distinct()

    hu <- hu %>%
      left_join(select(lp_outlet, .data$corrected_LevelPathI, updated_outlet_HUC12 = .data$outlet_HUC12),
                by = "corrected_LevelPathI") %>%
      select(-.data$outlet_HUC12) %>% rename(outlet_HUC12 = .data$updated_outlet_HUC12) %>%
      filter(.data$outlet_HUC12 != "" | .data$corrected_LevelPathI == -1) %>%
      mutate(outlet_HUC12 = ifelse(.data$corrected_LevelPathI == -1, .data$HUC12, .data$outlet_HUC12))
  }

  lp_outlet <- select(hu, .data$corrected_LevelPathI, .data$outlet_HUC12) %>%
    distinct() %>%
    group_by(.data$corrected_LevelPathI) %>%
    arrange(dplyr::desc(.data$outlet_HUC12)) %>%
    filter(dplyr::row_number() == 1) %>%
    ungroup()

  hu <- hu %>%
    left_join(select(lp_outlet, .data$corrected_LevelPathI, updated_outlet_HUC12 = .data$outlet_HUC12),
              by = "corrected_LevelPathI") %>%
    select(-.data$outlet_HUC12) %>% 
    rename(outlet_HUC12 = .data$updated_outlet_HUC12) %>%
    filter(.data$corrected_LevelPathI != -1) %>%
    distinct()

  ################################################################################
  # Add checks
  ################################################################################
  if(add_checks) {
    hu <- mutate(hu, trib_intersect = .data$HUC12 %in% hu_trib$HUC12,
                 trib_no_intersect = .data$HUC12 %in% hu_trib2$HUC12,
                 headwater_error = .data$HUC12 %in% funky_headwaters)
  }
  return(hu)
}

#' get_length_per_hu
#' @param net nhdplus
#' @param wbd wbd
#' @param simp simplification
#' @export
get_length_per_hu <- function(net, wbd, simp) {
  
  if(file.exists("temp_intersect.rds")) {
    intersect_net <- readRDS("temp_intersect.rds")
    
    warning("using cache")
    
  } else {
    
    wbd <- st_simplify(wbd, dTolerance = simp)
    
    intersect_net <- st_intersection(select(net, COMID), select(wbd, HUC12))
    
  }
  
  intersect_net$length <- st_length(intersect_net)
  
  return(intersect_net)
}
