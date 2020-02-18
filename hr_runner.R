library(nhdplusTools)
library(sf)
library(dplyr)
library(HUCAgg)
library(drake)
library(snow)
library(xml2)
library(purrr)
library(gifski)
library(mainstems)

source("R/3_setup.R")
source("R/4_find_match.R")
source("R/8_hr_proc.R")
source("R/9_viz_plot.R")

nhdhr_dir <- "data/hr/"
nhdhr_hu02 <- c("01", "02", "03", "04", "05", "06", "07", "08", "09",
               "10", "11", "12", "13", "14", "15", "16", "17", "18")
# nhdhr_hu02 <- c("10")
prj = 5070
nhdp_dir <- "data/nhdp"
nhdp_gdb <- "NHDPlusNationalData/NHDPlusV21_National_Seamless.gdb"
nhdp_url <- "https://s3.amazonaws.com/nhdplus/NHDPlusV21/Data/NationalData/NHDPlusV21_NationalData_CONUS_Seamless_Geodatabase_05.7z"
nhdp_gdb_path <- download_nhdplusv2(nhdp_dir)
# hr_path <- download_nhdplushr(nhdhr_dir, nhdhr_hu02)
hr_path <- "data/hr"
hr_vpus <- list.files(hr_path, pattern = ".*10[0-9][0-9].*.gdb$",
                         full.names = TRUE, recursive = TRUE, 
                         include.dirs = TRUE)

hr_vpus = rename_hr_fl(hr_vpus)

plan <- drake_plan(
  nhdp_net = target(get_net(read_sf(nhdp_gdb_path, "NHDFlowline_Network"), prj), hpc = FALSE),
  
  hr_vpu = target(nhdplusTools:::get_hr_data(vpu_f, "NHDFlowline"), 
                     transform = cross(vpu_f = !!hr_vpus), hpc = TRUE),
  
  vpu_prep = target(prep_nhdplushr(hr_vpu), 
                    transform = cross(hr_vpu), hpc = TRUE),
  
  nhdp_hw_outlets = target(get_hw_points(nhdp_net), hpc = FALSE),

  hr_pairs = target(get_hr_pairs(hr_path, nhdp_hw_outlets, 5070, 3), hpc = FALSE),

  matched_lp = target(match_flowpaths(source_flowline = nhdp_net,
                                      target_flowline = hr_vpu,
                                      hw_pair = hr_pairs),
                       transform = cross(hr_vpu), hpc = TRUE),

  compare = target(post_proc(v2_net = nhdp_net,
                             hr_net = hr_vpu,
                             matched_lp = matched_lp),
                   transform = map(hr_vpu, matched_lp), hpc = FALSE),

  combine = target(write_output(bind_rows(compare_2), "out/report/compare.csv"),
                   transform = combine(compare_2), hpc = FALSE),
  combine_matched_lp = target(write_output(bind_rows(matched_lp), "out/report/matched.csv"),
                              transform = combine(matched_lp), hpc = FALSE),

  compare_2 = target(compare_table(comp = compare,
                                 hr_net = hr_vpu,
                                 v2_net = nhdp_net),
                   transform = map(compare, hr_vpu), hpc = FALSE),
  compare_3 = target(post_analysis(comp = bind_rows(compare_2), lower_thresh_km = 1, upper_thresh_km = 100),
                     transform = combine(compare_2), hpc = FALSE),
  plot_bad = target(plot_fun(hr_net = hr_vpu, v2_net = nhdp_net,
                             lp = matched_lp, comp = compare_2, out_folder = "out/report"),
                    transform = map(matched_lp, compare_2, hr_vpu), hpc = FALSE),
  big_plot_data = target(get_plot_data(nhdp_net = nhdp_net, hu02 = nhdhr_hu02,
                                       order = 4, proj = 3857, simp_meters = 400),
                         hpc = FALSE),
  v2_plot_data = target(get_v2_plot_data(nhdp_net = nhdp_net, combine_matched_lp,
                                         proj = 3857, simp_meters = 400),
                        hpc = FALSE),
  hr_plot_data = target(get_hr_plot_data(hr_net = hr_vpu, lp = matched_lp,
                                         proj = 3857, simp_meters = 400),
                        transform = map(matched_lp, hr_vpu),
                        hpc = FALSE),
  plot_frames = target(create_plot_frames(big_plot_data, v2_plot_data,
                                          hr_plot_data, matched_lp, vpu_f,
                                          wbd_hr = "data/WBD_hr/WBDHU4_DA_06072019.shp",
                                          threads = 4),
                       transform = map(matched_lp, hr_plot_data, hr_vpu),
                       hpc = TRUE)
)

config <- drake_config(plan = plan,
                       memory_strategy = "autoclean", 
                       garbage_collection = TRUE, 
                       parallelism = "future", 
                       jobs = 4)

future::plan(future::multiprocess)

make(config = config)

# make(plan = plan, "compare_hr_vpu_data.hr.10.1025.gdb_matched_lp_hr_vpu_data.hr.10.1025.gdb")

# v2 LP: 590023769
