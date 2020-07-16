test_that("combo creates files", {
  seurat2ClusterMap("data/im/")
  expect_equal(dir.exists("data/im/clusterMapOutput/figures/"), T)
})
