source("R/0_setup.R")
source("R/1_download.R")
source("R/2_fixes.R")
source("R/3_setup.R")
source("R/4_find_match.R")
source("R/5_find_outlets.R")
source("R/6_visualize.R")
source("R/10_build_mainstems_table.R")
plan <- drake_plan(
  ##### Constants
  cores = 12,
  prj = 5070,
  viz_simp = 30,
  proc_simp = 10,
  national_viz_simp = 500,
  temp_dir = "temp/",
  nhdplus_dir = "data/nhdp/NHDPlusNationalData",
  nhdplus_update = file_in("data/nhdp/enhd_nhdplusatts.csv"),
  nhdplus_gdb = "NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07.7z",
  nhdp_v1_dir = "data/nhdpv1/",
  rf1_url = "https://water.usgs.gov/GIS/dsdl/erf1_2.e00.gz",
  rf1_dir = "data/RF1/",
  ##### Load static dependencies
  rf1_file = download_rf1(rf1_dir, rf1_url),
  rf1 = clean_rf1(read_sf(rf1_file)),
  nhdplus_cats = st_transform(read_sf(nhdplus_gdb_path, "CatchmentSP"), prj),
  nhdplus_wbd_fixes = get_fixes("nhdplusv2"),
  nhdplus_gdb_path = download_nhdplusv2(nhdplus_dir),
  nhdplus_wbd_exclusions = get_exclusions(nhdplus_gdb_path),
  nhdpv1_path = download_v1(nhdp_v1_dir),
  nhdpv1_cat = compile_v1_cats(nhdpv1_path, "data/nhdpv1/nhdpv1.gpkg"),
  nhdpv1_fline = compile_v1_fline(nhdpv1_path, "data/nhdpv1/nhdpv1.gpkg"),
  nhdpv1_atts = st_drop_geometry(nhdpv1_fline),
  nhdpv1_2_xwalk = get_nhdp_crosswalk(nhdplus_dir),
  nhdpv1_mapped = map_nhdpv1(mainstems_table, nhdpv1_atts, nhdplus_net_atts, nhdpv1_2_xwalk),
  ##### Load Data
  nhdplus_wbd = get_wbd(nhdplus_gdb_path, nhdplus_wbd_fixes, prj),
  nhdplus_net = get_net(read_sf(nhdplus_gdb_path, "NHDFlowline_Network"), prj, nhdplus_update),
  nhdplus_net_atts = st_set_geometry(nhdplus_net, NULL),
  ##### Match NHDPlusV2 with stable (old) WBD
  nhdplus_oldwbd_out = "out/nhdplus_oldwbd/",
  nhdplus_oldwbd_net_int = mainstems:::get_process_data(nhdplus_net, nhdplus_wbd, proc_simp),
  nhdplus_oldwbd_net_int_fix = nhdplus_oldwbd_net_int[!is.na(nhdplus_oldwbd_net_int$HUC12), ],
  nhdplus_oldwbd_hu_joiner = par_match_levelpaths(nhdplus_net, nhdplus_wbd, proc_simp, 
                                                  cores, "oldwbdtemp/", 
                                                  paste0(nhdplus_oldwbd_out, 
                                                         "map_joiner.csv"),
                                                  nhdplus_oldwbd_net_int_fix, FALSE,
                                                  split_terminals),
  nhdplus_oldwbd_hu_joiner_lengths = get_length_per_hu(nhdplus_net, nhdplus_wbd, proc_simp),
  nhdplus_oldwbd_hu_joiner_override = filter_dominant_length(nhdplus_oldwbd_hu_joiner, 
                                                             st_drop_geometry(nhdplus_oldwbd_hu_joiner_lengths), 
                                                             st_drop_geometry(readd(nhdplus_net)), 
                                                             factor_corrected = 2, 
                                                             factor_intersected = 10, 
                                                             remove_length_m = 200),
  nhdplus_oldwbd_lp_points = get_lp_points(nhdplus_oldwbd_hu_joiner_override, nhdplus_net, nhdplus_wbd, nhdplus_wbd_exclusions),
  nhdplus_oldwbd_na_outlet_coords = get_na_outlets_coords(nhdplus_oldwbd_lp_points$na, nhdplus_net),
  nhdplus_oldwbd_in_list = get_in_list(nhdplus_oldwbd_lp_points, nhdplus_net),
  nhdplus_oldwbd_linked_points = get_linked_points_scalable(nhdplus_oldwbd_in_list, nhdplus_oldwbd_na_outlet_coords, cores,
                                                   file.path(nhdplus_oldwbd_out, "wbd_viz.gpkg")),
  nhdplus_oldwbd_write = write_output_gpkg(nhdplus_net, nhdplus_wbd, nhdplus_oldwbd_hu_joiner_override,
                                           nhdplus_oldwbd_linked_points, prj, viz_simp, file.path(nhdplus_oldwbd_out, "wbd_viz.gpkg")),
  wbd_plumbing = get_hu_outlets(nhdplus_wbd, nhdplus_oldwbd_linked_points, file.path(nhdplus_oldwbd_out, "hu_outlets.gpkg")),
  #### Constants for newest WBD.
  wbd_dir = "data/wbd",
  wbd_zip_file = "WBD_National_GDB.zip",
  wbd_gdb_file = "WBD_National_GDB.gdb",
  wbd_url = "https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/WBD/National/GDB/WBD_National_GDB.zip",
  ##### Static dependencies for newest WBD
  wbd_fixes = get_fixes("latest"),
  wbd_gdb_path = download_wbd(wbd_dir),
  wbd_exclusions = get_exclusions(wbd_gdb_path[1]),
  wbd = get_wbd(wbd_gdb_path[1], wbd_fixes, prj),
  hu02 = get_hu02(wbd_gdb_path[1], prj, national_viz_simp),
  ##### Match newest WBD to NHDPlusV2.
  nhdplus_newwbd_out = "out/nhdplus_newwbd/",
  nhdplus_newwbd_net_int = mainstems:::get_process_data(nhdplus_net, wbd, proc_simp),
  nhdplus_newwbd_hu_joiner = par_match_levelpaths(nhdplus_net, wbd, proc_simp, cores, 
                                                  "newwbd_temp/", file.path(nhdplus_newwbd_out,
                                                                      "map_joiner.csv"), 
                                                  nhdplus_newwbd_net_int, FALSE,
                                                  split_terminals),
  nhdplus_newwbd_hu_joiner_lengths = get_length_per_hu(nhdplus_net, wbd, proc_simp),
  nhdplus_newwbd_hu_joiner_override = filter_dominant_length(nhdplus_newwbd_hu_joiner, 
                                                             st_drop_geometry(nhdplus_newwbd_hu_joiner_lengths), 
                                                             st_drop_geometry(readd(nhdplus_net)), 
                                                             factor_corrected = 2, 
                                                             factor_intersected = 10, 
                                                             remove_length_m = 200),
  nhdplus_newwbd_lp_points = get_lp_points(nhdplus_newwbd_hu_joiner_override, nhdplus_net, wbd, wbd_exclusions),
  nhdplus_newwbd_na_outlet_coords = get_na_outlets_coords(nhdplus_newwbd_lp_points$na, nhdplus_net),
  nhdplus_newwbd_in_list = get_in_list(nhdplus_newwbd_lp_points, nhdplus_net),
  nhdplus_newwbd_linked_points = get_linked_points_scalable(nhdplus_newwbd_in_list, nhdplus_newwbd_na_outlet_coords, cores,
                                                   file.path(nhdplus_newwbd_out, "wbd_viz.gpkg")),
  nhdplus_newwbd_write = write_output_gpkg(nhdplus_net, wbd, nhdplus_newwbd_hu_joiner_override,
                                           nhdplus_newwbd_linked_points, prj, viz_simp, file.path(nhdplus_newwbd_out, "wbd_viz.gpkg")),
  nhdplus_newwbd_plumbing = get_hu_outlets(wbd, nhdplus_newwbd_linked_points, file.path(nhdplus_newwbd_out, "wbd_viz.gpkg")),
  # Create plots for newest WBD matches.
  # plot_data = geom_plot_data(wbd, nhdplus_net, nhdplus_newwbd_hu_joiner, "^03.*"),
  # out_png = create_png(plot_data, nhdplus_newwbd_hu_joiner_override, "png/"),
  # out_wbd_plot_data = get_wbd_plot_data(nhdplus_net, wbd_gdb_path, nhdplus_newwbd_plumbing, viz_simp, prj, cores, file.path(nhdplus_newwbd_out, "wbd_plots.gpkg")),
  # out_wbd_plots = plot_wbd(out_wbd_plot_data),
  # Match RF1 to NHDPlusV2
  rf1_out = "out/rf1_out",
  rf1_hw = get_hw_points(rf1),
  rf1_nhdplus_hw_pairs = get_hw_pairs(rf1_hw, nhdplus_cats),
  rf1_nhdplus = match_flowpaths(nhdplus_net_atts,
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
  mainstems_with_v1 = add_v1(mainstems_table, nhdpv1_mapped),
  v1_mainstems = find_v1_mainstems(nhdpv1_atts, mainstems_with_v1, 6),
  mainstems_table_summary = make_ms_summary(mainstems_with_v1, nhdplus_net_atts),
  hist_list = get_hist_list(mainstems_table_summary),
  plot_lps_data = get_lp_plot_data(nhdplus_net, mainstems_table_summary, national_viz_simp),
  plot_lps_data_wbd = get_lp_plot_data_wbd(plot_lps_data, nhdplus_wbd, nhdplus_oldwbd_linked_points, national_viz_simp),
  plot_lps_data_all = get_lp_plot_data_rf1(plot_lps_data_wbd, rf1, rf1_nhdplus, national_viz_simp),
  plot_lps = get_lp_plots(plot_lps_data_all, 3, hu02, hu02_filter = "10",
                          bb = c(xmin = -103.5, ymin = 44.5, xmax = -101.5, ymax = 47)),
  plot_hw = get_hw_fig(),
  geo_summary = make_geo_summary(nhdplus_net, mainstems_table_summary, "out/mainstems_summary.gpkg"),
  upload = upload_sb(mainstems_table_summary, geo_summary)
)

make(plan = plan, memory_strategy = "autoclean", garbage_collection = TRUE)
