test_that("run_tangent works", {
  # Set 1
  sif_path <- testthat::test_path("test_data", "raw_inputs",  "set_1", "sample_information.txt")
  tsig_path <- testthat::test_path("test_data", "raw_inputs", "set_1", "tumor_log2RCN.txt")
  nsig_path <- testthat::test_path("test_data", "raw_inputs", "set_1", "normal_log2RCN.txt")

  res <- run_tangent(sif_path, nsig_path, tsig_path, n_latent = 1, make_plots = FALSE)
  expect <- readr::read_delim(testthat::test_path(
    "test_data", "set_1", "TangentXY",
    "TangentXYnormalized_tumor_log2RCN_1.txt"
  ), progress=FALSE, show_col_types=FALSE) %>% as.data.frame()
  expect_equal(res, expect)

  res <- run_tangent(sif_path, nsig_path, tsig_path, n_latent = 5, make_plots = FALSE)
  expect <- readr::read_delim(testthat::test_path(
    "test_data", "set_1", "TangentXY",
    "TangentXYnormalized_tumor_log2RCN_5.txt"
  ), progress=FALSE, show_col_types=FALSE) %>% as.data.frame()
  expect_equal(res, expect)

  res <- run_tangent(sif_path, nsig_path, tsig_path, n_latent = 10, make_plots = FALSE)
  expect <- readr::read_delim(testthat::test_path(
    "test_data", "set_1", "TangentXY",
    "TangentXYnormalized_tumor_log2RCN_10.txt"
  ), progress=FALSE, show_col_types=FALSE) %>% as.data.frame()
  expect_equal(res, expect)

  # Set 2
  sif_path <- testthat::test_path("test_data", "raw_inputs",  "set_2", "sample_information.txt")
  tsig_path <- testthat::test_path("test_data", "raw_inputs", "set_2", "tumor_log2RCN.txt")
  nsig_path <- testthat::test_path("test_data", "raw_inputs", "set_2", "normal_log2RCN.txt")

  res <- run_tangent(sif_path, nsig_path, tsig_path, n_latent = 1, make_plots = FALSE)
  expect <- readr::read_delim(testthat::test_path(
    "test_data", "set_2", "TangentXY",
    "TangentXYnormalized_tumor_log2RCN_1.txt"
  ), progress=FALSE, show_col_types=FALSE) %>% as.data.frame()
  expect_equal(res, expect)

  res <- run_tangent(sif_path, nsig_path, tsig_path, n_latent = 5, make_plots = FALSE)
  expect <- readr::read_delim(testthat::test_path(
    "test_data", "set_2", "TangentXY",
    "TangentXYnormalized_tumor_log2RCN_5.txt"
  ), progress=FALSE, show_col_types=FALSE) %>% as.data.frame()
  expect_equal(res, expect)

  res <- run_tangent(sif_path, nsig_path, tsig_path, n_latent = 10, make_plots = FALSE)
  expect <- readr::read_delim(testthat::test_path(
    "test_data", "set_2", "TangentXY",
    "TangentXYnormalized_tumor_log2RCN_10.txt"
  ), progress=FALSE, show_col_types=FALSE) %>% as.data.frame()
  expect_equal(res, expect)

  # Check that scrambling row order doesn't affect results
  sif <- readr::read_delim(sif_path, progress=FALSE, show_col_types=FALSE) %>% dplyr::slice_sample(prop = 1)
  tsig <- readr::read_delim(tsig_path, progress=FALSE, show_col_types=FALSE) %>% dplyr::slice_sample(prop = 1)
  nsig <- readr::read_delim(nsig_path, progress=FALSE, show_col_types=FALSE) %>% dplyr::slice_sample(prop = 1)
  res <- run_tangent(sif, nsig, tsig, n_latent = 10, make_plots = FALSE)
  expect_equal(res, expect)
})
