source("sourcer.R")
EC <- c("01", "02", "03")
MS <- c("07", "08")
OH <- c("05", "06")
MO <- c("10", "11")
MW <- c("17")
SO <- c("12", "13")
CO <- c("14", "15")
CA <- c("16", "18")

hr_hu02 <- list(EC = EC, MS = MS, OH = OH, MO = MO, MW = MW, SO = SO, CO = CO, CA = CA)

hr_hu02 <- list(OH = "06")
hr_hu02_n <- names(hr_hu02)
hr_dir <- "data/hr"
out <- "nhdplushr_newwbd"

plan <- drake_plan(
  ##### Constants
  cores = 1,
  prj = 5070,
  viz_simp = 30,
  proc_simp = 10,
  national_viz_simp = 500,
  temp_dir = "temp/",
  #### Constants for newest WBD.
  wbd_dir = "data/wbd",
  wbd_zip_file = "WBD_National_GDB.zip",
  wbd_gdb_file = "WBD_National_GDB.gdb",
  wbd_url = "https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/WBD/National/GDB/WBD_National_GDB.zip",
  ##### Static dependencies for newest WBD
  wbd_fixes = get_fixes("latest"),
  wbd_gdb_path = download_wbd(wbd_dir),
  wbd_exclusions = get_exclusions(wbd_gdb_path),
  wbd = get_wbd(wbd_gdb_path, wbd_fixes, prj),
  hu02 = st_simplify(st_transform(read_sf(wbd_gdb_path, "WBDHU2"), prj), dTolerance = national_viz_simp),
  hr_path = target(download_nhdplushr(hr_dir, hr_huset),
                      transform = map(hr_huset = !!hr_hu02, .id = FALSE)),
  hr_net = target(nhdhr_mod(hr_path, 
                               !!file.path(hr_dir, sprintf("%s.gpkg", .id_chr)), 
                               min_size = 6, simp = 10, proj = prj, force_terminal = TRUE, fix_terminals = TRUE),
                     transform = map(hr_path)),
  hu_joiner = target(par_match_levelpaths(hr_net, wbd, proc_simp, 1, temp_dir, out),
                     transform = cross(hr_net)),
  lp_points = target(get_lp_points(hu_joiner, hr_net, wbd, wbd_exclusions),
                     transform = map(hu_joiner, hr_net)),
  na_ol = target(get_na_outlets_coords(lp_points, hr_net),
                            transform = cross(lp_points)),
  in_list = target(get_in_list(lp_points, hr_net),
                   transform = cross(lp_points)),
  lp = target(get_linked_points_scalable(in_list, na_ol, cores, !!file.path(out, sprintf("%s.gpkg", hr_hu02_n))),
              transform = map(in_list, na_ol, hr_hu02_n = !!hr_hu02_n)),
  write = target(write_output_gpkg(hr_net, wbd, hu_joiner, lp, prj, viz_simp,
                                   file_out(!!file.path(out, sprintf("%s.gpkg", hr_hu02_n)))),
                 transform = map(hr_net, hu_joiner, lp, hr_hu02_n = !!hr_hu02_n))
)

config <- drake_config(plan = plan,
                       memory_strategy = "autoclean",
                       garbage_collection = TRUE)

make(config = config)
