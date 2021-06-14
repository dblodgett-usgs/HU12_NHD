source("R/0_setup.R")
source("R/1_download.R")
source("R/3_setup.R")

plan <- drake_plan(cores = 12, 
                   naturalearth_dir = "data/naturalearth/", 
                   naturalearth_shp = download_naturalearth(naturalearth_dir),
                   naturalearth_data = sf::read_sf(naturalearth_shp),
                   merit_dir = "data/merit/",
                   merit_out = "out/merit/merit_plus.gpkg",
                   merit_raw = get_merit(merit_dir),
                   merit_cats = get_merit_cats(merit_dir),
                   merit_atts = get_merit_atts(merit_raw, cores, merit_cache = merit_out),
                   merit = write_merit(merit_raw, merit_atts, merit_out))

make(plan)
