test_that("transform_normals() correctly performs linear transformation on input normals", {
  # Set 1
  sif_path <- testthat::test_path("test_data", "raw_inputs",  "set_1", "sample_information.txt")
  pif_path <- testthat::test_path("test_data", "raw_inputs", "set_1", "probe_information.txt")
  tsig_path <- testthat::test_path("test_data", "raw_inputs", "set_1", "tumor_log2RCN.txt")

  tangent_res_1 <- readr::read_delim(testthat::test_path(
    "test_data", "set_1", "TangentXY",
    "TangentXYnormalized_tumor_log2RCN_1.txt"
  ), progress=FALSE, show_col_types=FALSE) %>% as.data.frame()
  tangent_res_5 <- readr::read_delim(testthat::test_path(
    "test_data", "set_1", "TangentXY",
    "TangentXYnormalized_tumor_log2RCN_5.txt"
  ), progress=FALSE, show_col_types=FALSE) %>% as.data.frame()
  tangent_res_10 <- readr::read_delim(testthat::test_path(
    "test_data", "set_1", "TangentXY",
    "TangentXYnormalized_tumor_log2RCN_10.txt"
  ), progress=FALSE, show_col_types=FALSE) %>% as.data.frame()
  tangent_res <- list(tangent_res_1, tangent_res_5, tangent_res_10)

  res <- plot_signal_noise(tangent_res, sif_path, pif_path, tsig_path, n_latent = c(1, 5, 10))
  expect <- readr::read_delim(testthat::test_path(
    "test_data", "set_1", "SignalNoise", "SignalNoise.txt"
  ), progress=FALSE, show_col_types=FALSE)
  expect_equal(res, as.data.frame(expect))

  # Set 2
  sif_path <- testthat::test_path("test_data", "raw_inputs",  "set_2", "sample_information.txt")
  pif_path <- testthat::test_path("test_data", "raw_inputs", "set_2", "probe_annotation.txt")
  tsig_path <- testthat::test_path("test_data", "raw_inputs", "set_2", "tumor_log2RCN.txt")

  tangent_res_1 <- readr::read_delim(testthat::test_path(
    "test_data", "set_2", "TangentXY",
    "TangentXYnormalized_tumor_log2RCN_1.txt"
  ), progress=FALSE, show_col_types=FALSE) %>% as.data.frame()
  tangent_res_5 <- readr::read_delim(testthat::test_path(
    "test_data", "set_2", "TangentXY",
    "TangentXYnormalized_tumor_log2RCN_5.txt"
  ), progress=FALSE, show_col_types=FALSE) %>% as.data.frame()
  tangent_res_10 <- readr::read_delim(testthat::test_path(
    "test_data", "set_2", "TangentXY",
    "TangentXYnormalized_tumor_log2RCN_10.txt"
  ), progress=FALSE, show_col_types=FALSE) %>% as.data.frame()
  tangent_res <- list(tangent_res_1, tangent_res_5, tangent_res_10)

  res <- plot_signal_noise(tangent_res, sif_path, pif_path, tsig_path, n_latent = c(1, 5, 10))
  expect <- readr::read_delim(testthat::test_path(
    "test_data", "set_2", "SignalNoise", "SignalNoise.txt"
  ), progress=FALSE, show_col_types=FALSE)
  expect_equal(res, as.data.frame(expect))
})
