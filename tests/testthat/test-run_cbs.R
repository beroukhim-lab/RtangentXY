test_that("run_cbs() correctly runs with 1 latent factor Tangent", {
  tangent_out <- readRDS(testthat::test_path("dummydata", "TangentXYnormalized_tumor_log2RCN_1latentFactors.rds"))
  cbs_out <- run_cbs(tangent_out)

  cbs_expect_out <- readRDS(testthat::test_path("dummydata", "tumor_log2RCN_1latentFactors_TANGENT_CBS.rds"))
  testthat::expect_equal(cbs_out, cbs_expect_out)
})

test_that("run_cbs() correctly runs with 5 latent factor Tangent", {
  tangent_out <- readRDS(testthat::test_path("dummydata", "TangentXYnormalized_tumor_log2RCN_5latentFactors.rds"))
  cbs_out <- run_cbs(tangent_out)

  cbs_expect_out <- readRDS(testthat::test_path("dummydata", "tumor_log2RCN_5latentFactors_TANGENT_CBS.rds"))
  testthat::expect_equal(cbs_out, cbs_expect_out)
})

test_that("run_cbs() correctly runs with 10 latent factor Tangent", {
  tangent_out <- readRDS(testthat::test_path("dummydata", "TangentXYnormalized_tumor_log2RCN_10latentFactors.rds"))
  cbs_out <- run_cbs(tangent_out)

  cbs_expect_out <- readRDS(testthat::test_path("dummydata", "tumor_log2RCN_10latentFactors_TANGENT_CBS.rds"))
  testthat::expect_equal(cbs_out, cbs_expect_out)
})
