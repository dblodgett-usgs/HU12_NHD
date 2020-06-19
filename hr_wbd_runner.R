source("R/0_setup.R")
source("R/2_fixes.R")
source("R/3_setup.R")
source("R/4_find_match.R")
source("R/5_find_outlets.R")
source("R/6_visualize.R")
source("R/10_build_mainstems_table.R")

hr_hu02 <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", 
             "10", "11", "17", "12", "13", "14", "15", "16", "18")

pattern <- ".*[0-9][0-9][0-9][0-9].*.gdb$"

hr_dir <- "data/hr/"
out <- "out/nhdplushr_newwbd"

wbd_fixes <- get_fixes("latest") %>%
  bind_rows(list(HUC12 = "180102040904", TOHUC = "180102041003", comment = "misdirected"))

cores <- NA
prj <- 5070
min_size_sqkm <- 6
proc_simp <- 2
temp_dir <- paste0("temp", hr_hu02, "/")
wbd_dir <- "data/wbd"
wbd_url <- "https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/WBD/National/GDB/WBD_National_GDB.zip"

plan <- drake_plan(
  wbd_gdb_path = download_wbd(wbd_dir, url = wbd_url),
  wbd_exclusions = get_exclusions(wbd_gdb_path[1]),
  wbd = get_wbd(wbd_gdb_path[1], wbd_fixes, prj),
  hr_path = target(download_nhdplushr(hr_dir, hr_huset),
                      transform = map(hr_huset = !!hr_hu02)),
  hr_net = target(get_nhdplushr(hr_path, 
                               !!file.path(hr_dir, sprintf("%s.gpkg", .id_chr)), 
                               layers = "NHDFlowline",
                               pattern = pattern, 
                               min_size_sqkm = min_size_sqkm, simp = proc_simp, 
                               proj = prj, check_terminals = TRUE),
                     transform = map(hr_path, .id = FALSE), hpc = FALSE),
  hr_net_fix = target(fix_hr(hr_net), transform = map(hr_net, .id = FALSE), hpc = FALSE),
  hu_joiner = target(par_match_levelpaths(hr_net_fix, wbd, proc_simp, 1, temp_dir, 
                                          !!file.path(out, sprintf("joiner_%s.csv", .id_chr))),
                     transform = map(hr_net_fix, prep_data, temp_dir = !!temp_dir, .id = FALSE), hpc = TRUE),
  lp_points = target(get_lp_points(hu_joiner, hr_net_fix, wbd, wbd_exclusions),
                     transform = map(hu_joiner, hr_net_fix, .id = FALSE), hpc = TRUE),
  na_ol = target(get_na_outlets_coords(lp_points, hr_net_fix),
                            transform = map(lp_points, .id = FALSE), hpc = TRUE),
  in_list = target(get_in_list(lp_points, hr_net_fix),
                   transform = map(lp_points, .id = FALSE), hpc = TRUE),
  lp = target(get_linked_points_scalable(in_list, na_ol, cores, !!file.path(out, sprintf("%s.gpkg", hr_hu02))),
              transform = map(in_list, na_ol, hr_hu02 = !!hr_hu02, .id = FALSE), hpc = TRUE),
  write = target(write_output_gpkg(hr_net_fix, wbd, hu_joiner, lp, prj, 30,
                                   file_out(!!file.path(out, sprintf("%s.gpkg", hr_hu02)))),
                 transform = map(hr_net_fix, hu_joiner, lp, hr_hu02 = !!hr_hu02,.id = FALSE), hpc = TRUE)
)

config <- drake_config(plan = plan,
                       memory_strategy = "autoclean",
                       garbage_collection = TRUE)

future::plan(future::multiprocess)

make(config = config)

all_outlets <- do.call(rbind, lapply(hr_hu02, function(x) {
  read_sf(file.path(out, paste0(x, ".gpkg")), "linked_points")
}))

all_outlets <- st_transform(all_outlets, "+proj=longlat +datum=WGS84 +no_defs")
names(all_outlets)[1] <- "NHDPlusID"
write_sf(all_outlets, file.path(out, "all_outlets.gpkg"))

all_joiner <- do.call(rbind, lapply(c("", paste0("_", c(2:length(hr_hu02)))), function(x) {
  read_csv(file.path(out, paste0("joiner_hu_joiner", x, ".csv")))
}))

all_joiner <- write_csv(all_joiner, file.path(out, "all_joiner.csv"))

sb <- sbtools::authenticate_sb()

sbtools::item_append_files("5dd73836e4b06957976525c7", 
                           c(file.path(out, "all_joiner.csv"), 
                             file.path(out, "all_outlets.gpkg")), 
                           session = sb)


sbtools::item_append_files("5dd73836e4b06957976525c7", file.path(out, "all_outlets.gpkg"), session = sb)

