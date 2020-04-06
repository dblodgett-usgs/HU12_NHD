# HU12_NHD R source

This directory contains R source files for the HU12_NHD workflow. The project uses a drake workflow which uses numerous R packages and the functions contained in this folder.

- `0_setup.R`: Ensure all required packages are installed and loaded.
- `1_download.R`: Download functions not included in nhdplusTools.
- `2_fixes.R`: Modifications to HUC12 "TOHUC" codes required to make network function.
- `3_setup.R`: Data prep functions to rename and clean up inputs.
- `4_find_match.R`: Pre-processing functions for headwater and outlet matching
- `5_find_outlets.R`: Mostly moved to mainstems package. WBD Exclusions are still in this file.
- `6_visualize.R`: Output gpkg writers and functions to create png/gif plots.
- `7_compare.R`: (Not called by Drake -- could remove) comparison between work here and another method. 
- `8_hr_proc.R`: NHDPlusHR-Specific workflow functions for output post processing and visualiztion.
- `9_viz_plot.R`: Plotting functions for hr workflow.
- `10_build_mainstems_table.R`: Functions to build final summary output.
