library(nhdplusTools)
library(hyRefactor)
library(sf)
library(dplyr)
library(HUCAgg)
library(drake)
library(snow)
library(xml2)
library(readr)
library(rmapshaper)
library(R.utils)
library(tidyr)
library(igraph)
library(readr)

source("R/1_download_data.R")
source("R/2_fixes.R")
source("R/3_setup.R")
source("R/4_find_match.R")
source("R/5_find_outlets.R")
source("R/6_visualize.R")
source("R/10_build_mainstems_table.R")

plan <- drake_plan(
  ##### Constants
  cores = 3,
  prj = 5070,
  viz_simp = 30,
  proc_simp = 10,
  national_viz_simp = 500,
  temp_dir = "temp/",
  nhdplus_dir = "data/nhdp/",
  nhdplus_gdb = "NHDPlusNationalData/NHDPlusV21_National_Seamless.gdb",
  nhdplus_url = "https://s3.amazonaws.com/nhdplus/NHDPlusV21/Data/NationalData/NHDPlusV21_NationalData_CONUS_Seamless_Geodatabase_05.7z",
  rf1_url = "https://water.usgs.gov/GIS/dsdl/erf1_2.e00.gz",
  rf1_dir = "data/RF1/",
  rf1_file = download_rf1(rf1_url, rf1_dir),
  rf1 = clean_rf1(read_sf(rf1_file)),
  nhdplus_cats = st_transform(read_sf(nhdplus_gdb_path, "CatchmentSP"), prj),
  ##### Load static dependencies
  nhdplus_wbd_fixes = get_fixes("nhdplusv2"),
  nhdplus_gdb_path = download_nhd_gdb(nhdplus_dir, nhdplus_gdb, nhdplus_url),
  nhdplus_wbd_exclusions = get_exclusions(nhdplus_gdb_path),
  ##### Load Data
  nhdplus_wbd = get_wbd(nhdplus_gdb_path, nhdplus_wbd_fixes, prj),
  nhdplus_net = get_net(read_sf(nhdplus_gdb_path, "NHDFlowline_Network"), prj),
  nhdplus_net_atts = st_set_geometry(nhdplus_net, NULL),
  ##### Match NHDPlusV2 with stable (old) WBD
  nhdplus_oldwbd_out = "nhdplus_oldwbd_out/",
  nhdplus_oldwbd_hu_joiner = get_hu_joiner(nhdplus_net, nhdplus_wbd, proc_simp, cores, temp_dir, nhdplus_oldwbd_out),
  nhdplus_oldwbd_lp_points = get_lp_points(nhdplus_oldwbd_hu_joiner, nhdplus_net, wbd, wbd_exclusions),
  nhdplus_oldwbd_na_outlet_coords = get_na_outlets_coords(nhdplus_oldwbd_lp_points, nhdplus_net),
  nhdplus_oldwbd_in_list = get_in_list(nhdplus_oldwbd_lp_points, nhdplus_net),
  nhdplus_oldwbd_linked_points = get_linked_points(nhdplus_oldwbd_in_list, nhdplus_oldwbd_na_outlet_coords, cores,
                                                   file.path(nhdplus_oldwbd_out, "wbd_viz.gpkg")),
  nhdplus_oldwbd_write = write_output_gpkg(nhdplus_net, nhdplus_wbd, nhdplus_oldwbd_hu_joiner,
                                           nhdplus_oldwbd_linked_points, prj, viz_simp, nhdplus_oldwbd_out),
  wbd_plumbing = get_hu_outlets(nhdplus_wbd, nhdplus_oldwbd_linked_points, file.path(nhdplus_oldwbd_out, "wbd_viz.gpkg")),
  #### Constants for newest WBD.
  wbd_dir = "data/wbd",
  wbd_zip_file = "WBD_National_GDB.zip",
  wbd_gdb_file = "WBD_National_GDB.gdb",
  wbd_url = "https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/WBD/National/GDB/WBD_National_GDB.zip",
  ##### Static dependencies for newest WBD
  wbd_fixes = get_fixes("latest"),
  wbd_gdb_path = download_wbd(wbd_dir, wbd_zip_file, wbd_gdb_file, wbd_url),
  wbd_exclusions = get_exclusions(wbd_gdb_path),
  wbd = get_wbd(wbd_gdb_path, wbd_fixes, prj),
  ##### Match newest WBD to NHDPlusV2.
  nhdplus_newwbd_out = "nhdplus_newwbd_out/",
  nhdplus_newwbd_hu_joiner = get_hu_joiner(nhdplus_net, wbd, proc_simp, cores, temp_dir, nhdplus_newwbd_out),
  nhdplus_newwbd_lp_points = get_lp_points(nhdplus_newwbd_hu_joiner, nhdplus_net, wbd, wbd_exclusions),
  nhdplus_newwbd_na_outlet_coords = get_na_outlets_coords(nhdplus_newwbd_lp_points, nhdplus_net),
  nhdplus_newwbd_in_list = get_in_list(nhdplus_newwbd_lp_points, nhdplus_net),
  nhdplus_newwbd_linked_points = get_linked_points(nhdplus_newwbd_in_list, nhdplus_newwbd_na_outlet_coords, cores,
                                                   file.path(nhdplus_newwbd_out, "wbd_viz.gpkg")),
  nhdplus_newwbd_write = write_output_gpkg(nhdplus_net, wbd, nhdplus_newwbd_hu_joiner,
                                           nhdplus_newwbd_linked_points, prj, viz_simp, nhdplus_newwbd_out),
  nhdplus_newwbd_plumbing = get_hu_outlets(wbd, nhdplus_newwbd_linked_points, file.path(nhdplus_newwbd_out, "wbd_viz.gpkg")),
  # Create plots for newest WBD matches.
  plot_data = geom_plot_data(wbd, nhdplus_net, nhdplus_newwbd_hu_joiner, "^03.*"),
  out_png = create_png(plot_data, nhdplus_newwbd_hu_joiner, "png/"),
  out_wbd_plot_data = get_wbd_plot_data(nhdplus_net, wbd_gdb_path, nhdplus_newwbd_plumbing, viz_simp, prj, cores, file.path(nhdplus_newwbd_out, "wbd_plots.gpkg")),
  out_wbd_plots = plot_wbd(out_wbd_plot_data),
  # Match RF1 to NHDPlusV2
  rf1_out = "rf1_out",
  rf1_hw = get_hw_points(rf1),
  rf1_nhdplus_hw_pairs = get_hw_pairs(rf1_hw, nhdplus_cats),
  rf1_nhdplus = match_flowpaths(left_join(select(prepare_nhdplus(nhdplus_net_atts, 100, 20, 40, FALSE), COMID),
                                          nhdplus_net_atts, by = "COMID"), #[, 1:40]
                                st_set_geometry(rf1, NULL),
                                rf1_nhdplus_hw_pairs, 4),
  rf1_output = write_rf1_output(rf1, rf1_nhdplus, rf1_out),
  mainstems_table = build_mainstem_table(nhdplus_net_atts, 
                                         nhdplus_oldwbd_linked_points, 
                                         nhdplus_newwbd_linked_points,
                                         rf1_nhdplus,
                                         rf1,
                                         st_set_geometry(nhdplus_wbd, NULL),
                                         st_set_geometry(wbd, NULL)),
  mainstems_table_summary = make_ms_summary(mainstems_table, nhdplus_net_atts),
  hist_list = get_hist_list(mainstems_table_summary),
  plot_lps_data = get_lp_plot_data(nhdplus_net, mainstems_table_summary, national_viz_simp),
  plot_lps_data_wbd = get_lp_plot_data_wbd(plot_lps_data, nhdplus_wbd, nhdplus_oldwbd_linked_points, national_viz_simp),
  plot_lps_data_all = get_lp_plot_data_rf1(plot_lps_data_wbd, rf1, rf1_nhdplus, national_viz_simp),
  plot_lps = get_lp_plots(plot_lps_data_all, 2)
  ##### NHDPlsuHR Stuff
  # nhdhr_hu02 = c("01", "02"),
  # nhdhr_dir = "data/hr",
  # nhdhr_path = download_nhdhr(nhdhr_dir, nhdhr_hu02),
  # nhdhr = nhdhr_mod(nhdhr_path, file.path(nhdhr_dir, "nhdplushr.gpkg"), force_terminal = TRUE),
  # nhdhr_net = get_net(nhdhr, prj),
  # nhdplushr_newwbd_out = "nhdplushr_newwbd",
  # nhdplushr_newwbd_hu_joiner = get_hu_joiner(nhdhr_net, wbd, proc_simp, cores, temp_dir, nhdplushr_newwbd_out),
  # nhdplushr_newwbd_linked_points = get_linked_points(nhdplushr_newwbd_hu_joiner, nhdhr_net, wbd, wbd_exclusions, cores,
  #                                                    file.path(nhdplushr_newwbd_out, "wbd_viz.gpkg")),
  # nhdplushr_newwbd_write = write_output_gpkg(nhdhr_net, wbd, nhdplushr_newwbd_hu_joiner, 
  #                                           nhdplus_oldwbd_linked_points, prj, viz_simp, nhdplushr_newwbd_out),
)

config <- drake_config(plan = plan,
                       memory_strategy = "autoclean", 
                       garbage_collection = TRUE)

make(config = config)
