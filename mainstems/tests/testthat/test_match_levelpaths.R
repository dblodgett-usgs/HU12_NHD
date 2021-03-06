context("match levelpaths")
clean_huc12 <- function(huc12) {
  bad_lps <- select(huc12, corrected_LevelPathI, trib_no_intersect, headwater_error) %>%
    filter(trib_no_intersect | headwater_error)

  huc12 %>%
    filter(!corrected_LevelPathI %in% bad_lps$corrected_LevelPathI & !is.na(outlet_HUC12)) %>%
    distinct()
}

test_that("match levelpaths and get_linked_points", {
  source(system.file("extdata/new_hope_data.R", package = "nhdplusTools"))

  suppressWarnings(net_prep <- select(new_hope_flowline, COMID, DnLevelPat, AreaSqKM) %>%
    left_join(prepare_nhdplus(new_hope_flowline,
                              min_network_size = 20, # sqkm
                              min_path_length = 0, # sqkm
                              min_path_size = 10, # sqkm
                              purge_non_dendritic = TRUE,
                              warn =  TRUE), by = "COMID") %>%
    st_sf() %>%
    group_by(LevelPathI) %>%
    arrange(Hydroseq) %>%
    mutate(DnLevelPat = DnLevelPat[1]) %>%
    ungroup())

  net_prep["denTotalAreaSqKM"] <-
    nhdplusTools::calculate_total_drainage_area(select(st_set_geometry(net_prep, NULL),
                                         ID = COMID, toID = toCOMID,
                                         area = AreaSqKM))

  wbd <- select(new_hope_wbd, HUC12 = HUC_12, TOHUC = HU_12_DS) %>%
    st_transform(st_crs(net_prep))

  net_prep <- st_join(net_prep, wbd) %>%
    st_set_geometry(NULL)

  outlet_comid <- 8897784

  hu_lp <- match_levelpaths(net_prep, outlet_comid)

  linked_points <- get_linked_points(hu_lp, new_hope_flowline, wbd, exclude = c())

  expect_equal(nrow(linked_points), 8)
  
  new_hope_flowline$TerminalPa <- min(new_hope_flowline$LevelPathI)
  
  temp <- file.path(tempdir(check = TRUE), "mainstems_temp")
  unlink(temp, recursive = TRUE)
  dir.create(temp, recursive = TRUE)
  
  suppressWarnings(hu_lp_2 <- par_match_levelpaths(new_hope_flowline, wbd, 10, 1, 
                                                   temp_dir = temp, 
                                                   out_file = file.path(temp, "temp.csv")))

  expect_true(all(names(hu_lp) %in% names(hu_lp_2)))
  
  expect_true(nrow(hu_lp) == nrow(hu_lp_2))
  
  split_terminals <- data.frame(COMID = c(8894324, 8897784),
                                stop = c(0, 8894324))
  
  unlink(temp, recursive = TRUE)
  dir.create(temp, recursive = TRUE)
  
  suppressWarnings(hu_lp_2 <- par_match_levelpaths(new_hope_flowline, wbd, 10, 1, 
                                                   temp_dir = temp, 
                                                   out_file = file.path(temp, "temp.csv"),
                                                   split_terminals = split_terminals))
  
  expect_true(all(names(hu_lp) %in% names(hu_lp_2)))
  
  expect_true(nrow(hu_lp) == nrow(hu_lp_2))
  
})

test_that("match levelpaths runs 2279159", {
  start_comid <- 2279159
  net_prep <- readRDS(system.file("extdata/match_levelpaths_2279159.rds", package = "hyRefactor"))
  matched <- match_levelpaths(net_prep, start_comid, add_checks = TRUE)
  expect_true(matched$intersected_LevelPathI[which(matched$HUC12 == "102702020404")] == 550002171)
  expect_true(matched$intersected_LevelPathI[which(matched$HUC12 == "102702050201")] == 550002171)
  expect_true(matched$intersected_LevelPathI[which(matched$HUC12 == "102702020205")] == 550002171)

  # This tests for bad toHUC corrections where the intersection doesn't match the
  # TOHUC coding.
  expect_true(matched$intersected_LevelPathI[which(matched$HUC12 == "102702020102")] == 550031020)
  expect_true(matched$trib_no_intersect[which(matched$HUC12 == "102702020102")])

  # This could be fixed by fixing toHUC codes, but will test it as is for now to verify behavior.
  expect_true(matched$outlet_HUC12[which(matched$HUC12 == "102702020102")] == "102702050705")

  expect_equal(sum(matched$trib_no_intersect), 2)

  # Corrected levelpath and head huc when first order tributary gets corrected.
  expect_true(matched$corrected_LevelPathI[matched$HUC12 == "102702050701"] == 550049901)
  expect_true(matched$head_HUC12[matched$HUC12 == "102702050701"] == "102702050701")

  # Corrected head huc when part of a multi-HU tributary.
  expect_true(matched$head_HUC12[matched$HUC12 == "102702060404"] == "102702060403")

  # If we remove the trib_no_intersect errors this should be unique.
  matched <- filter(matched, !trib_no_intersect & matched$corrected_LevelPathI != -1) %>%
    select(corrected_LevelPathI, head_HUC12, outlet_HUC12) %>%
    distinct()
  expect_true(length(which(duplicated(matched$corrected_LevelPathI))) == 0)
})

test_that("match levelpaths multi-overlap-outlet 931010009", {
  matched <- match_levelpaths(readRDS("data/match_levelpaths_931010009.rds"), 931010009)
  expect_true(nrow(matched) == 4)

  huc12 <- dplyr::select(matched, levelpath = corrected_LevelPathI, head_huc12 = head_HUC12, outlet_huc12 = outlet_HUC12) %>%
    dplyr::filter(!is.na(outlet_huc12)) %>%
    dplyr::distinct()

  expect_true(length(unique(huc12$levelpath)) == nrow(huc12))
})

test_that("match levelpaths funky heatwater 4292649", {
  net <- read_sf("data/test_4292649.gpkg", "NHDFlowline_Network") %>%
    st_zm()
  hu <- read_sf("data/test_4292649.gpkg", "HUC12_new")
    
  matched <- par_match_levelpaths(net, hu, 0, 1, 
                                  tempdir(check = TRUE), 
                                  "temp.csv")
  
  unlink("temp.csv")
  
  #verified manually.
  expect_true(matched$corrected_LevelPathI[matched$HUC12 == "010100040901"] == 150020702)

  # not much to do with this one. 010100040905 has multiple overlaps from the wrong watershed.
  expect_true(!150067066 %in% matched$corrected_LevelPathI)
  expect_equal(nrow(matched), 155)
  expect_equal(sum(matched$headwater_error), 2)

  huc12 <- dplyr::select(clean_huc12(matched), corrected_LevelPathI, head_HUC12, outlet_HUC12) %>%
    dplyr::filter(!is.na(outlet_HUC12)) %>%
    dplyr::filter(!is.na(corrected_LevelPathI)) %>%
    dplyr::distinct()

  expect_true(length(unique(huc12$corrected_LevelPathI)) == nrow(huc12))

  # strange headwater behavior.
  expect_true(huc12$outlet_HUC12[huc12$corrected_LevelPathI == 150014576] == "010100030308")
  expect_true(huc12$head_HUC12[huc12$corrected_LevelPathI == 150014576] == "010100030308")
})

test_that("match levelpaths 20204804", {
  matched <- match_levelpaths(readRDS("data/match_levelpaths_20204804.rds"), 20204804, add_checks = TRUE)
  expect_equal(sum(matched$trib_no_intersect), 5)

  # goofy closed basin breaks things.
  expect_true(matched$intersected_LevelPathI[matched$HUC12 == "180901030502"] == 10004787)
  expect_true(matched$corrected_LevelPathI[matched$HUC12 == "180901030502"] == 10030198)
})

test_that("match levelpaths 10055266", {
  # headwaters of levelpath 200011667
  matched <- match_levelpaths(readRDS("data/match_levelpaths_10055266.rds"), 
                              10055266, add_checks = TRUE)
  expect_true(all(matched$head_HUC12[matched$corrected_LevelPathI == 200011667] == "020802010601"))

  expect_equal(nrow(matched), 297)
  
  huc12 <- dplyr::select(matched, levelpath = corrected_LevelPathI, head_huc12 = head_HUC12, outlet_huc12 = outlet_HUC12) %>%
    dplyr::filter(!is.na(outlet_huc12)) %>%
    dplyr::distinct()

  expect_true(length(unique(huc12$levelpath)) == nrow(huc12))

})

test_that("match levelpaths 10390202", {
  # uncorrected outlet_HUC12 on levelpath: 800015797
  matched <- match_levelpaths(readRDS("data/match_levelpaths_10390202.rds"), 10390202, add_checks = TRUE)

  huc12 <- dplyr::select(matched,
                         levelpath = corrected_LevelPathI,
                         head_huc12 = head_HUC12,
                         outlet_huc12 = outlet_HUC12) %>%
    dplyr::filter(!is.na(outlet_huc12)) %>%
    dplyr::distinct()

  expect_true(length(unique(huc12$levelpath)) == nrow(huc12))

})

test_that("match levelpaths runs 12228521", {
  start_comid <- 12228521
  net_prep <- readRDS("data/match_levelpaths_12228521.rds")
  matched <- match_levelpaths(net_prep, start_comid, add_checks = TRUE)

  expect_true(length(matched$head_HUC12[matched$HUC12 == "040601040203"]) == 0)
})

test_that("match levelpaths doesn't miss levelpath 250031924", {
  start_comid <- 9074086
  net_prep <- readRDS("data/match_levelpaths_9074086.rds")
  matched <- match_levelpaths(net_prep, start_comid, add_checks = TRUE)
  
  expect_true(matched$corrected_LevelPathI[matched$HUC12 == "030402060502"] == "250031924")
})

test_that("high_res problem oputlet", {
  start_comid <- 5000200014717
  net_prep <- readRDS("data/hr_lp_poroblem.rds")

  matched <- match_levelpaths(net_prep, start_comid, add_checks = TRUE)

  expect_true(matched$head_HUC12[matched$HUC12 == "010200051009"] == "010200020201")
})

test_that("alaska scale issues", {
  start_comid <- 81020021
  net_prep <- readRDS(list.files(pattern = "*81020021.rds",
                                 full.names = TRUE, recursive = TRUE))

  net_prep <- dplyr::mutate(net_prep, TOHUC = ifelse(HUC12 == "190801081608", 
                                                     "190801081700", TOHUC))
  
  warnings <- capture_warnings(matched <- match_levelpaths(net_prep, start_comid, add_checks = TRUE))
  
  expect_equal(length(warnings), 9)
  
  expect_equal(nrow(matched), 4246)
})

test_that("duplicate in r01", {
  start_comid <- 1737140
  
  net_prep <- readRDS(list.files(pattern = "*1737140.rds",
                                 full.names = TRUE, 
                                 recursive = TRUE))
  
  net_prep <- dplyr::mutate(net_prep, TOHUC = ifelse(HUC12 == "010200040401", "010200040402", TOHUC))
  net_prep <- dplyr::mutate(net_prep, TOHUC = ifelse(HUC12 == "010200040603", "010200040501", TOHUC))
  net_prep <- dplyr::mutate(net_prep, TOHUC = ifelse(HUC12 == "010200040501", "010200040503", TOHUC))
  
  
  matched <- match_levelpaths(net_prep, start_comid, add_checks = TRUE)
})
