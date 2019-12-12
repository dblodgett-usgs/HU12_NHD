source("sourcer.R")
source("R/3_setup.R")
source("R/4_find_match.R")
source("R/5_find_outlets.R")
source("R/6_visualize.R")
source("R/10_build_mainstems_table.R")
hr_hu02 <- c("01", "02", "03", "07", "08", "05", "06", "10", 
        "11", "17", "12", "13", "14", "15", "16", "18")
# hr_hu02 <- c("06")
hr_dir <- "data/hr/"
out <- "nhdplushr_newwbd"

fixes <- read_csv("fixes/hu_fixes.csv") %>%
  bind_rows(list(HUC12 = "180102040904", TOHUC = "180102041003", comment = "misdirected"))

plan <- drake_plan(
  ##### Constants
  cores = NA,
  prj = 5070,
  viz_simp = 30,
  proc_simp = 2,
  national_viz_simp = 500,
  temp_dir = "temp/",
  #### Constants for newest WBD.
  wbd_dir = "data/wbd",
  wbd_zip_file = "WBD_National_GDB.zip",
  wbd_gdb_file = "WBD_National_GDB.gdb",
  wbd_url = "https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/WBD/National/GDB/NATIONAL_WBD_GDB.zip",
  ##### Static dependencies for newest WBD
  wbd_fixes = fixes,
  wbd_gdb_path = download_wbd(wbd_dir, url = wbd_url),
  wbd_exclusions = get_exclusions(wbd_gdb_path),
  wbd = get_wbd(wbd_gdb_path, wbd_fixes, prj),
  hu02 = st_simplify(st_transform(read_sf(wbd_gdb_path, "WBDHU2"), prj), dTolerance = national_viz_simp),
  hr_path = target(download_nhdplushr(hr_dir, hr_huset),
                      transform = map(hr_huset = !!hr_hu02)),
  hr_net = target(get_nhdplushr(hr_path, 
                               !!file.path(hr_dir, sprintf("%s.gpkg", .id_chr)), layers = "NHDFlowline",
                               patter = ".*[0-9][0-9][0-9][0-9].*.gdb$", min_size_sqkm = 6, simp = 2, 
                               proj = prj, check_terminals = TRUE),
                     transform = map(hr_path)),
  hu_joiner = target(par_match_levelpaths(hr_net, wbd, proc_simp, 1, temp_dir, 
                                          !!file.path(out, sprintf("joiner_%s.csv", .id_chr))),
                     transform = map(hr_net)),
  lp_points = target(get_lp_points(hu_joiner, hr_net, wbd, wbd_exclusions),
                     transform = map(hu_joiner, hr_net)),
  na_ol = target(get_na_outlets_coords(lp_points, hr_net),
                            transform = cross(lp_points)),
  in_list = target(get_in_list(lp_points, hr_net),
                   transform = cross(lp_points)),
  lp = target(get_linked_points_scalable(in_list, na_ol, cores, !!file.path(out, sprintf("%s.gpkg", hr_hu02))),
              transform = map(in_list, na_ol, hr_hu02 = !!hr_hu02)),
  write = target(write_output_gpkg(hr_net, wbd, hu_joiner, lp, prj, viz_simp,
                                   file_out(!!file.path(out, sprintf("%s.gpkg", hr_hu02)))),
                 transform = map(hr_net, hu_joiner, lp, hr_hu02 = !!hr_hu02))
)


config <- drake_config(plan = plan,
                       memory_strategy = "autoclean",
                       garbage_collection = TRUE,
                       parallelism = "future", 
                       jobs = 8)

make(config = config)

all_outlets <- do.call(rbind, lapply(hr_hu02, function(x) {
  read_sf(file.path(out, paste0(x, ".gpkg")), "linked_points")
}))

all_outlets <- st_transform(all_outlets, "+proj=longlat +datum=WGS84 +no_defs")
names(all_outlets)[1] <- "NHDPlusID"
write_sf(all_outlets, file.path(out, "all_outlets.gpkg"))

all_joiner <- do.call(rbind, lapply(hr_hu02, function(x) {
  read_csv(file.path(out, paste0("joiner_hu_joiner_hr_net_hr_path_", x, ".csv")))
}))

all_joiner <- write_csv(all_joiner, file.path(out, "all_joiner.csv"))

sb <- sbtools::authenticate_sb()

sbtools::item_append_files("5dd73836e4b06957976525c7", 
                           c(file.path(out, "all_joiner.csv"), 
                             file.path(out, "all_outlets.gpkg")), 
                           session = sb)


sbtools::item_append_files("5dd73836e4b06957976525c7", file.path(out, "all_outlets.gpkg"), session = sb)

