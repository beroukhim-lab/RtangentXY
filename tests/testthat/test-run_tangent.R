test_that("run_tangent() correctly runs with 1 latent factor", {
  sif_path <- testthat::test_path("dummydata", "raw_inputs", "sample_information.txt")
  normal_path <- testthat::test_path("dummydata", "raw_inputs", "normal_log2RCN.txt")
  tumor_path <- testthat::test_path("dummydata", "raw_inputs", "tumor_log2RCN.txt")

  run_tangent_out <- run_tangent(sif_path, normal_path, tumor_path, 1)

  run_tangent_expect_out <- readRDS(testthat::test_path("dummydata", "TangentXYnormalized_tumor_log2RCN_1latentFactors.rds"))
  testthat::expect_equal(run_tangent_out, run_tangent_expect_out)
})

test_that("run_tangent() correctly runs with 5 latent factors", {
  sif_path <- testthat::test_path("dummydata", "raw_inputs", "sample_information.txt")
  normal_path <- testthat::test_path("dummydata", "raw_inputs", "normal_log2RCN.txt")
  tumor_path <- testthat::test_path("dummydata", "raw_inputs", "tumor_log2RCN.txt")

  run_tangent_out <- run_tangent(sif_path, normal_path, tumor_path, 5)

  run_tangent_expect_out <- readRDS(testthat::test_path("dummydata", "TangentXYnormalized_tumor_log2RCN_5latentFactors.rds"))
  testthat::expect_equal(run_tangent_out, run_tangent_expect_out)
})

test_that("run_tangent() correctly runs with 10 latent factors", {
  sif_path <- testthat::test_path("dummydata", "raw_inputs", "sample_information.txt")
  normal_path <- testthat::test_path("dummydata", "raw_inputs", "normal_log2RCN.txt")
  tumor_path <- testthat::test_path("dummydata", "raw_inputs", "tumor_log2RCN.txt")

  run_tangent_out <- run_tangent(sif_path, normal_path, tumor_path, 10)

  run_tangent_expect_out <- readRDS(testthat::test_path("dummydata", "TangentXYnormalized_tumor_log2RCN_10latentFactors.rds"))
  testthat::expect_equal(run_tangent_out, run_tangent_expect_out)
})
