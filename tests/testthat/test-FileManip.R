
# comb
test_that("comboRDS: load File adds experimentName", {
  ds = loadFile("./data/im/input/comb.rds")
  expect_equal(getExperimentFromSeurat(ds), "im")
})

test_that("comboRDS: load File adds datasetName", {
  ds = loadFile("./data/im/input/comb.rds")
  expect_equal(getNameFromSeurat(ds), "comb")
})

test_that("comboRDS: getFname works", {
  expect_equal(getFname("./data/im/input/comb.rds"), "comb")
})