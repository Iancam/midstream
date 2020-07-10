test_that("Pipe Runs Without Errors", {
  tenX2Seurat("data/im", percent.mt = 3, minNFeature = 2000)
})

test_that("QCFiles exist", {
  expect_equal(file.exists("data/im/seuratOutput/d1/afterThresh.vln.pdf"), T)
})


test_that("PCAFiles exist", {
  expect_equal(file.exists("data/im/seuratOutput/d1/dim.loadings.pdf"), T)
})

# file.remove("./data/im/input/comb.rds")
# unlink("data/im/seuratOutput")