test_that("convert_cbs() correctly runs with 1 latent factor Tangent", {
  cbs_out <- readRDS(testthat::test_path("test_data", "tumor_log2RCN_1latentFactors_TANGENT_CBS.rds"))
  sif_path <- testthat::test_path("test_data", "raw_inputs", "set_1", "sample_information.txt")
  tumor_path <- testthat::test_path("test_data", "raw_inputs", "set_1", "tumor_log2RCN.txt")

  cbs_reformatted <- convert_cbs(cbs_out, sif_path, tumor_path)
  cbs_reformatted_expected <- readRDS(testthat::test_path("test_data", "tumor_log2RCN_1latentFactors_TANGENT_CBS_REF.rds"))

  testthat::expect_equal(cbs_reformatted, cbs_reformatted_expected)
})

test_that("convert_cbs() correctly runs with 5 latent factor Tangent", {
  cbs_out <- readRDS(testthat::test_path("test_data", "tumor_log2RCN_5latentFactors_TANGENT_CBS.rds"))
  sif_path <- testthat::test_path("test_data", "raw_inputs", "set_1", "sample_information.txt")
  tumor_path <- testthat::test_path("test_data", "raw_inputs", "set_1", "tumor_log2RCN.txt")

  cbs_reformatted <- convert_cbs(cbs_out, sif_path, tumor_path)
  cbs_reformatted_expected <- readRDS(testthat::test_path("test_data", "tumor_log2RCN_5latentFactors_TANGENT_CBS_REF.rds"))

  testthat::expect_equal(cbs_reformatted, cbs_reformatted_expected)
})

test_that("convert_cbs() correctly runs with 10 latent factor Tangent", {
  cbs_out <- readRDS(testthat::test_path("test_data", "tumor_log2RCN_10latentFactors_TANGENT_CBS.rds"))
  sif_path <- testthat::test_path("test_data", "raw_inputs", "set_1", "sample_information.txt")
  tumor_path <- testthat::test_path("test_data", "raw_inputs", "set_1", "tumor_log2RCN.txt")

  cbs_reformatted <- convert_cbs(cbs_out, sif_path, tumor_path)
  cbs_reformatted_expected <- readRDS(testthat::test_path("test_data", "tumor_log2RCN_10latentFactors_TANGENT_CBS_REF.rds"))

  testthat::expect_equal(cbs_reformatted, cbs_reformatted_expected)
})
