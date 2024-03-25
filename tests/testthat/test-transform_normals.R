test_that("transform_normals() correctly performs linear transformation on input normals", {
  sif_path <- testthat::test_path("dummydata", "raw_inputs", "sample_information.txt")
  normal_path <- testthat::test_path("dummydata", "raw_inputs", "normal_log2RCN.txt")

  tn_test <- transform_normals(sif_path, normal_path)

  tn_expect_out <- readRDS(testthat::test_path("dummydata", "normal_log2RCN_linearTrans.rds"))
  testthat::expect_equal(tn_test, tn_expect_out)
})
