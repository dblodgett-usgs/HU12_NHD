get_exclusions <- function(wbd_gdb) {
  
  if("HUC12" %in% st_layers(wbd_gdb)$name) {
    wbd <- read_sf(wbd_gdb, "HUC12")
    wbd <- rename(wbd, HUC12 = HUC_12, TOHUC = HU_12_DS)
  } else {
    wbd <- read_sf(wbd_gdb, "WBDHU12")
  }
  
  names(wbd)[names(wbd) != attr(wbd, "sf_column")] <- toupper(names(wbd)[names(wbd) != attr(wbd, "sf_column")])

  wbd_type <- st_set_geometry(wbd, NULL)
    
  if("HUTYPE" %in% names(wbd_type)) {
    wbd_type <- distinct(select(wbd_type, HUC12, HUTYPE))
  } else {
    wbd_type <- distinct(select(wbd_type, HUC12, HUTYPE = HU_12_TYPE))
  }
  
  wbd <- sf::st_drop_geometry(wbd) %>%
    group_by(HUC12) %>%
    summarise(TOHUC = TOHUC[1])
  
  # Exclusions where river-flow does not apply:
  exclude_type <- wbd_type$HUC12[wbd_type$HUTYPE %in% c("F", "I", "C", "U")] # frontal closed or island
  exclude_first_order_toHUC <- wbd$HUC12[wbd$TOHUC %in% c("OCEAN", "CANADA", "GEATLAKES", "UNKNOWN") & 
                                           !wbd$HUC12 %in% wbd$TOHUC] # Unless it has something flowing to it.
  
  exclude <- unique(c(exclude_type, exclude_first_order_toHUC))
  
  return(exclude)
}
