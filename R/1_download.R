download_v1 <- function(out_dir) {
  dir.create(out_dir, showWarnings = FALSE)
  system(paste('wget -P', out_dir, '-r -np -nc -A "*atshape.zip" -A "*NHD.zip" ftp://ftp.horizon-systems.com/NHDplus/NHDPlusV1/'))
  out_dir
}

compile_v1_fline <- function(out_dir, out_gpkg) {
  layer <- "flowline"
  if(!file.exists(out_gpkg) || !layer %in% st_layers(out_gpkg)) {
    fl <- get_fl(out_dir)
    
    fline_files <- fl[grepl("nhdflowline.shp", fl, ignore.case = TRUE)]
    vaa_files <- fl[grepl(".*VAA.dbf$", fl, ignore.case = TRUE)]
    
    flines <- lapply(fline_files, function(x) st_zm(read_sf(x)))
    
    flines[[21]] <- dplyr::select(flines[[21]], -Shape_Le_1)
    
    flines <- st_sf(bind_rows(flines))
    
    vaa <- lapply(vaa_files, read.dbf)
    
    vaa <- bind_rows(vaa)
    
    flines <- left_join(flines, vaa, by = "COMID")
    
    write_sf(flines, out_gpkg, layer)
  } else {
    flines <- read_sf(out_gpkg, layer)
  }
  return(flines)
}

get_fl <- function(out_dir) {
  fl_zip <- list.files(out_dir, recursive = TRUE, full.names = TRUE, pattern = ".zip")
  
  dev_null <- lapply(fl_zip, unzip, overwrite = FALSE, exdir = out_dir)
  
  fl <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
  
  for(i in seq_along(fl)) {
    new_fi <- NULL
    
    if(grepl("nhdf", basename(fl[i]), ignore.case = FALSE)) new_fi <- gsub("nhdf", "NHDF", fl[i])
    if(grepl("nhda", basename(fl[i]), ignore.case = FALSE)) new_fi <- gsub("nhda", "NHDA", fl[i])
    if(grepl("nhdw", basename(fl[i]), ignore.case = FALSE)) new_fi <- gsub("nhdw", "NHDW", fl[i])
    if(grepl("nhdp", basename(fl[i]), ignore.case = FALSE)) new_fi <- gsub("nhdp", "NHDP", fl[i])
    if(grepl("nhdl", basename(fl[i]), ignore.case = FALSE)) new_fi <- gsub("nhdl", "NHDL", fl[i])
    
    if(!is.null(new_fi)) {
      file.rename(fl[i], new_fi)
      
      fl[i] <- new_fi
    }
  }
  
  fl[!fl %in% fl_zip]
}

compile_v1_cats <- function(out_dir, out_gpkg) {
  layer <- "catchment"
  if(!file.exists(out_gpkg) || !layer %in% st_layers(out_gpkg)) {

    fl <- get_fl(out_dir)
    
    cat_files <- fl[grepl("catchment.shp", fl, ignore.case = TRUE)]
    
    cats <- lapply(cat_files, read_sf)
    
    cats <- st_sf(bind_rows(cats))
    
    write_sf(cats, out_gpkg, layer)
  } else {
    cats <- read_sf(out_gpkg, layer)
  }
  return(cats)
}

get_nhdp_crosswalk <- function(nhdplus_dir) {
  foreign::read.dbf(file.path(nhdplus_dir, "NHDPlusV1Network_V2Network_Crosswalk.dbf"))
}

get_merit_basins <- function(merit_dir) {
  dir.create(merit_dir, recursive = TRUE, showWarnings = FALSE)
  url <- "http://hydrology.princeton.edu/data/mpan/MERIT_Basins/MERIT_Hydro_v07_Basins_v01/zip/pfaf_level_01/"
  f1 <- "pfaf_"
  f2 <- "_MERIT_Hydro_v07_Basins_v01.zip"
  out <- c()
  for(i in 1:9) {
    f <- file.path(merit_dir, paste0(f1, i, f2))
    if(!file.exists(f)) {
      download.file(paste0(url, f1, i, f2), 
                    destfile = f)
      unzip(f, exdir = merit_dir)
    }
    out <- c(out, f)
  }
}

