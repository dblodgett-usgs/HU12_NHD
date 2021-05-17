#' Get headwater points
#' @description Returns headwater locations for line geometries 
#' @param fline flowline data.frame from NHDPlus or with an ID and toID column.
#' @importFrom sf st_sf st_coordinates st_as_sf st_crs st_drop_geometry
#' @importFrom dplyr left_join select filter group_by ungroup bind_cols
#' @importFrom nhdplusTools prepare_nhdplus
#' @export
get_hw_points <- function(fline) {
  if("COMID" %in% names(fline)) {
    
    if(!"StreamOrde" %in% names(fline)) {
      fline$StreamOrde <- 1
      fline$StreamCalc <- 1
    }
    
    fline <- select(fline, .data$COMID) %>%
      right_join(prepare_nhdplus(fline, 0, 0, 0, warn = FALSE), 
                by = "COMID") %>%
      st_sf() %>%
      filter(!.data$COMID %in% .data$toCOMID)
  } else {
    fline <- filter(fline, !.data$ID %in% .data$toID)
  }

  outlets <-  fline %>%
    st_coordinates() %>%
    as.data.frame()
  
  if("L2" %in% names(outlets)) {
    outlets <- group_by(outlets, .data$L2)
  } else {
    outlets <- group_by(outlets, .data$L1)
  }
  
  outlets <- outlets %>%
    filter(row_number() == round(n()/2)) %>%
    ungroup() %>%
    select(.data$X, .data$Y) %>%
    st_as_sf(coords = c("X", "Y"))

  if("COMID" %in% names(fline)) {
    bind_cols(outlets, select(st_drop_geometry(fline), .data$COMID)) %>%
      st_sf(crs = st_crs(fline))
  } else {
    bind_cols(outlets, select(st_drop_geometry(fline), .data$ID)) %>%
      st_sf(crs = st_crs(fline))
  }
}

clean_geom <- function(x) {
  if("sf" %in% class(x)) {
    st_set_geometry(x, NULL)
  } else {
    x
  }
}

#' Match Flowpaths
#' @description Implements a flowpath-matching algorithm that traces downstream along
#' the target flowline network and determines which levelpath from the source flowlines
#' best matches the resulting downstream traces. The algorithm starts from the outlet
#' location of the upstream most catchment in the source flowlines to stay away from
#' complexity that occurs near drainage divides.
#'
#' @param source_flowline sf data.frame with COMID and LevelPathI attributes.
#' @param target_flowline sf data.frame either NHDPlusHR with NHDPlusID, LENGTHKM, DnHydroseq, 
#' Hydroseq, and LevelPathI or an ad-hoc network with ID and toID attributes.
#' @param hw_pair data.frame with headwater pairs
#' @param cores integer number of cores to use in NHDPlus downstream network navigation.
#' @return data.frame with the headwater id from the source_flowline input and 
#' levelpath from the source_flowline input for each identifier matched from the 
#' target dataset.
#' @export
#' @importFrom sf st_join st_set_geometry st_within
#' @importFrom tidyr unnest
#' @importFrom dplyr select distinct  left_join bind_rows
#' @importFrom nhdplusTools get_DM align_nhdplus_names
#' @importFrom pbapply pblapply pboptions
#' @examples
#' source(system.file("extdata/nhdplushr_data.R", package = "mainstems"))
#' source(system.file("extdata/new_hope_data.R", package = "nhdplusTools"))
#'
#' hr_catchment <- nhdplusTools::align_nhdplus_names(hr_catchment)
#' hw_pair <- sf::st_set_geometry(sf::st_join(get_hw_points(new_hope_flowline),
#'                               dplyr::select(hr_catchment, FEATUREID),
#'                               join = sf::st_within), NULL)
#'
#' lp_df_df <- match_flowpaths(new_hope_flowline, hr_flowline, hw_pair)
#' matched <- dplyr::left_join(dplyr::select(hr_flowline, COMID),
#'                             dplyr::select(lp_df_df, member_NHDPlusID,
#'                                           MR_LevelPathI = mr_LevelPathI), 
#'                                           by = c("COMID" = "member_NHDPlusID"))
#'
#' lp <- min(matched$MR_LevelPathI, na.rm = TRUE)
#' mr_lp <- dplyr::filter(new_hope_flowline, LevelPathI <= lp)
#' hr_lp <- dplyr::filter(matched, MR_LevelPathI <= lp)
#' plot(sf::st_geometry(matched), col = "blue", lwd = 0.5)
#' plot(sf::st_geometry(mr_lp), col = "red", lwd = 3, add = TRUE)
#' plot(sf::st_geometry(hr_lp), col = "black", add = TRUE)
#'
match_flowpaths <- function(source_flowline, target_flowline, hw_pair, cores = NA) {

  source_flowline <- align_nhdplus_names(source_flowline)
  source_flowline <- clean_geom(source_flowline)
  target_flowline <- clean_geom(target_flowline)
  
  if("ID" %in% names(target_flowline)) {
  
    hw_pair <- filter(hw_pair, !is.na(.data$FEATUREID) & !is.na(.data$ID) & 
                        .data$FEATUREID %in% source_flowline$COMID)
  
    source_flowline <- source_flowline[source_flowline$TerminalPa %in% 
                                         source_flowline[source_flowline$COMID %in% 
                                                           hw_pair$FEATUREID, ]$TerminalPa, ]
    
    source_flowline <- select(source_flowline, c("COMID", "LevelPathI", "DnLevelPat", "DnHydroseq", "Hydroseq"))
      
    lpt <- data.frame(headwater_COMID = hw_pair$FEATUREID)
    
    cl <- NULL
    if(!is.na(cores)) {
      cl <- parallel::makeCluster(rep("localhost", cores), type = "SOCK")
    }
    
    lpt["member_ID"] <- list(pblapply(hw_pair$ID,
                                      function(x, fa) get_dwn(x, fa),
                                      fa = target_flowline, cl = cl))
    
    lpt <- unnest(lpt, cols = c("member_ID"))
    
    lps <- data.frame(headwater_COMID = hw_pair$FEATUREID)
    
    if(!is.na(cores)) {
      parallel::stopCluster(cl)
      cl <- parallel::makeCluster(rep("localhost", cores), type = "SOCK")
    }
    
    lps["member_COMID"] <- list(pblapply(hw_pair$FEATUREID,
                                         function(x, fa) get_DM(fa, x),
                                         fa = source_flowline, cl = cl))
    
    lps <- unnest(lps, cols = c("member_COMID"))
    
    # With all navigations downstream main we need to find the 
    # dominant downstream levelpath for each headwater location.
    # First filter all results so we have only the top of each 
    # levelpath found.
    lps <- lps %>%
      left_join(select(source_flowline, .data$LevelPathI, .data$Hydroseq, .data$COMID), 
                by = c("member_COMID" = "COMID")) %>%
      group_by(.data$LevelPathI) %>%
      filter(.data$Hydroseq == max(.data$Hydroseq)) %>% 
      ungroup()
    
    # Also grab the most upstream by headwater_Hydroseq to deal with assumption noted below.
    lps <- left_join(lps, select(source_flowline, .data$COMID, headwater_Hydroseq = .data$Hydroseq),
                     by = c("headwater_COMID" = "COMID")) %>%
      group_by(.data$LevelPathI) %>%
      filter(.data$headwater_Hydroseq == max(.data$headwater_Hydroseq)) %>%
      ungroup()
    
    # Now for each headwater, choose the minimum (largest) levelpath.
    # This assumes that the largest levelpaths extend upstream of the
    # rest of the network.
    lps <- group_by(lps, .data$headwater_COMID) %>%
      filter(.data$LevelPathI == min(.data$LevelPathI)) %>%
      select(-.data$Hydroseq, mr_LevelPathI = .data$LevelPathI, 
             -.data$member_COMID, -.data$headwater_Hydroseq) %>%
      ungroup()
    
    target_flowline <- select(clean_geom(target_flowline), .data$ID)
    source_flowline <- distinct(select(source_flowline, .data$COMID, .data$LevelPathI))
    
    gc()
    
    hw_pair <- hw_pair %>%
      rename(COMID = .data$FEATUREID) %>%
      left_join(lps, by = c("COMID" = "headwater_COMID"))
    
    # Join so we have RF1 IDs found downstream of a given MR LevelPath headwater.
    lpt <- left_join(lpt, select(hw_pair, -.data$ID),
                     by = c("headwater_COMID" = "COMID"))
    
    group_by(lpt, .data$member_ID) %>%
      filter(!is.na(.data$mr_LevelPathI)) %>%
      filter(.data$mr_LevelPathI == min(.data$mr_LevelPathI)) %>%
      ungroup()
      
  } else {

  target_flowline <- select(clean_geom(target_flowline),
                            NHDPlusID = .data$COMID, 
                            HydroSeq = .data$Hydroseq, 
                            DnHydroSeq = .data$DnHydroseq, 
                            LevelPathI = .data$LevelPathI, 
                            DnLevelPat = .data$DnLevelPat)
  
  gc()

  hw_pair <- rename(hw_pair, NHDPlusID = .data$FEATUREID) %>%
    filter(.data$NHDPlusID %in% target_flowline$NHDPlusID)
  
  dm_NHDPlusID <- lapply(hw_pair$NHDPlusID,
                         function(x, fa) get_DM(fa, x),
                         fa = target_flowline)
  
  target_flowline <- select(target_flowline, .data$NHDPlusID)

  # Expand into data.frame
  lp_df <- data.frame(headwater_COMID = hw_pair$COMID)
  lp_df["member_NHDPlusID"] <- list(dm_NHDPlusID)

  lp_df <- unnest(lp_df, cols = c("member_NHDPlusID"))

  # Get MR levelpaths for headwater COMIDs
  hw_pair <- left_join(hw_pair, 
                       rename(source_flowline, mr_LevelPathI = .data$LevelPathI), 
                       by = "COMID")

  # Join so we have HR NHDPlusIDs found downstream of a given MR LevelPath headwater.
  lp_df <- left_join(lp_df, select(hw_pair, -.data$NHDPlusID), 
                     by = c("headwater_COMID" = "COMID"))

  group_by(lp_df, .data$member_NHDPlusID) %>% # Group by level paths present in HR.
    filter(.data$mr_LevelPathI == min(.data$mr_LevelPathI)) %>% # Filter so only one (largest) MR levelpath is linked to each HR path.
    ungroup() %>%
    left_join(target_flowline,
              by = c("member_NHDPlusID" = "NHDPlusID")) %>%
    select(-.data$headwater_COMID)
  }
}

get_dwn_rec <- function(ID, target_flowline) {
  next_dn <- target_flowline$toID[target_flowline$ID == ID]
  if(is.na(next_dn)) return()
  c(ID, get_dwn(next_dn, target_flowline))
}

get_dwn <- function(ID, target_flowline) {
  next_dn <- target_flowline$toID[target_flowline$ID == ID]
  if(is.na(next_dn)) {
    return(ID)
  } else {
    return(get_dwn_rec(ID, target_flowline))
  }
}
