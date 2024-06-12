test_that("plot_signal_noise() properly returns signal noise matrix on filepaths", {
  sif_path <- testthat::test_path("dummydata", "raw_inputs", "sample_information.txt")
  pif_path <- testthat::test_path("dummydata", "raw_inputs", "probe_information.txt")
  normal_path <- testthat::test_path("dummydata", "raw_inputs", "normal_log2RCN.txt")
  tumor_path <- testthat::test_path("dummydata", "raw_inputs", "tumor_log2RCN.txt")

  psn_out <- plot_signal_noise(sif_path, pif_path, normal_path, tumor_path, c(1,5,10))$signal_noise_df

  psn_expect_out <- readRDS(testthat::test_path("dummydata", "SignalNoise.rds"))
  testthat::expect_equal(psn_out, psn_expect_out)
})

test_that("plot_signal_noise() properly returns signal noise matrix on R objects", {
  sif <- readRDS(testthat::test_path("dummydata", "raw_inputs", "sample_information.rds"))
  pif <- readRDS(testthat::test_path("dummydata", "raw_inputs", "probe_information.rds"))
  n.df <- readRDS(testthat::test_path("dummydata", "raw_inputs", "normal_log2RCN.rds"))
  t.df <- readRDS(testthat::test_path("dummydata", "raw_inputs", "tumor_log2RCN.rds"))

  psn_out <- plot_signal_noise(sif, pif, n.df, t.df, c(1,5,10))$signal_noise_df

  psn_expect_out <- readRDS(testthat::test_path("dummydata", "SignalNoise.rds"))
  testthat::expect_equal(psn_out, psn_expect_out)
})
