test_that("run_svd() correctly performs SVD on input data", {
  sif_path <- testthat::test_path("dummydata", "raw_inputs", "sample_information.txt")
  transformed_normal_path <- testthat::test_path("dummydata", "normal_log2RCN_linearTrans.txt")

  run_svd_out <- run_svd(sif_path, transformed_normal_path)

  run_svd_expect_out <- readRDS(testthat::test_path("dummydata", "n.autox.svd.rds"))
  testthat::expect_equal(run_svd_out, run_svd_expect_out)
})
