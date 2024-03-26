test_that("run_svd() correctly performs SVD on input data", {
  sif_path <- testthat::test_path("dummydata", "raw_inputs", "sample_information.txt")
  normal_path <- testthat::test_path("dummydata", "raw_inputs", "normal_log2RCN.txt")

  run_svd_out <- run_svd(sif_path, normal_path)

  run_svd_expect_out <- readRDS(testthat::test_path("dummydata", "n.autox.svd.rds"))
  testthat::expect_equal(run_svd_out, run_svd_expect_out)
})
