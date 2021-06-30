source("R/0_setup.R")
source("R/1_download.R")
source("R/3_setup.R")
source("R/4_find_match.R")
source("R/10_build_mainstems_table.R")

plan <- drake_plan(cores = 12, 
                   naturalearth_dir = "data/naturalearth/", 
                   naturalearth_shp = download_naturalearth(naturalearth_dir),
                   naturalearth_data = sf::read_sf(naturalearth_shp),
                   naturalearth_data_trim = get_trimmed(naturalearth_data),
                   naturalearth_data_ends = get_fp_ends(naturalearth_data_trim),
                   merit_dir = "data/merit/",
                   merit_out = "out/merit/merit_plus.gpkg",
                   merit_raw = get_merit(merit_dir),
                   merit_cats = get_merit_cats(merit_dir),
                   merit_atts = get_merit_atts(merit_raw, cores, merit_cache = merit_out),
                   merit = write_merit(merit_raw, merit_atts, merit_out),
                   merit_natearth_ends = join_naturalearth_ends(naturalearth_data_ends,
                                                                merit),
                   merit_natearth_heads = get_naturalearth_heads(merit_natearth_ends, 
                                                                 merit),
                   merit_w_names = add_names(merit_natearth_heads, 
                                             select(st_drop_geometry(naturalearth_data),
                                                    name, rivernum),
                                             merit_atts),
                   merit_w_dm = add_dm(merit_w_names),
                   write_merit_share_ = write_merit_share(merit, merit_w_dm, "out/merit/merit_plus_simplify.gpkg"),
                   upload = upload_sb_fun(list(`merit_plus_simplify.gpkg` = "out/merit/merit_plus_simplify.gpkg")))

make(plan, memory_strategy = "autoclean", garbage_collection = TRUE)
