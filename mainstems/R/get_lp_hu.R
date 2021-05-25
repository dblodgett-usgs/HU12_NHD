get_lp_hu <- function(fline_hu, start_comid, one_level = FALSE) {
  
  nlp <- unique(filter(fline_hu, .data$COMID == start_comid)[["LevelPathI"]])
  
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
      
      if(one_level) {
        keep_going <- FALSE
        check <- FALSE
      } else {
        
        keep_going <- FALSE
        
        i <- 0
        while(length(next_lp) == 0 & length(nlp_tracker) > 0) {
          
          next_lp <- get_next_lp(fline_hu, nlp_tracker)
          
          next_lp <- unique(next_lp[["LevelPathI"]])
          
          i <- i + 1
          if(i > 1000) stop("runaway loop?")
          
          if(length(nlp_tracker) > 1) { # maintain backlog that needs to be worked through.
            nlp_tracker <- nlp_tracker[2:length(nlp_tracker)]
          } else {
            nlp_tracker <- c()
          }
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
  
  if(one_level) {
    return(list(fline_hu, 
                data.frame(LevelPathI = names(lp_hu),
                           HUC12 = I(lp_hu),
                           stringsAsFactors = FALSE) %>%
                  tidyr::unnest(cols = c("HUC12"))))
  }
  
  return(data.frame(LevelPathI = names(lp_hu),
                    HUC12 = I(lp_hu),
                    stringsAsFactors = FALSE) %>%
           tidyr::unnest(cols = c("HUC12")))
}