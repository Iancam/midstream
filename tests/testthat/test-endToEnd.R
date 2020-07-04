test_that("Pipe Runs Without Errors", {
  tenX2Seurat("data/im", percent.mt = 1, minNFeature = 2000)
})

test_that("QCFiles exist", {
  file.exists("data/im/seuratOutput/d1/d1.afterThresh.vln.pdf")
})


test_that("PCAFiles exist", {
  file.exists("data/im/seuratOutput/d1/d1.dim.loadings.pdf")
})

file.remove("./data/im/input/comb.rds")
unlink("data/im/seuratOutput")