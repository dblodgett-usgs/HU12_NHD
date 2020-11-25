test_that("headwater terminals discriminate correctly", {
  load("data/disciminate_outlets.rda")
  
  hr_net_sub <- filter(hr_net_sub, FCode != 33600)
  hu_joiner_sub <- filter(hu_joiner_sub, corrected_LevelPathI %in% hr_net_sub$LevelPathI)
  
  lp <- get_lp_points(hu_joiner_sub, hr_net_sub, wbd_sub, wbd_exclusions_sub)
  
  expect_equal(lp$lp[lp$lp$hu12 == "170501030205", ]$lp, 55000700010254)
  expect_equal(lp$lp[lp$lp$hu12 == "170501030202", ]$lp, 55000700000702)

})

test_that("outlet linke points", {
  load("data/lp_outlet_point.rda")

  point <- get_lp_points(joiner, net, hu, "")  
  
  expect_equal(nrow(point$lp), 4)
})