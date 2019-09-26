# nolint start
proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

hr_path <- system.file("extdata/03_sub.gpkg", package = "mainstems")

hr_catchment <- sf::read_sf(hr_path, "NHDPlusCatchment")
hr_catchment <- sf::st_transform(hr_catchment, proj)

hr_flowline <- nhdplusTools:::get_hr_data(hr_path, "NHDFlowline")
hr_flowline <- sf::st_transform(hr_flowline, proj)

hr_wbd <- sf::read_sf(system.file("extdata/new_hope_wbd.gpkg", package = "nhdplusTools"),
                      "HUC12")
# nolint end
