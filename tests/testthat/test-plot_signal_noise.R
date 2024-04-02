test_that("plot_signal_noise() properly returns signal noise matrix", {
  sif_path <- testthat::test_path("dummydata", "raw_inputs", "sample_information.txt")
  pif_path <- testthat::test_path("dummydata", "raw_inputs", "probe_annotation.txt")
  normal_path <- testthat::test_path("dummydata", "raw_inputs", "normal_log2RCN.txt")
  tumor_path <- testthat::test_path("dummydata", "raw_inputs", "tumor_log2RCN.txt")

  psn_out <- plot_signal_noise(sif_path, pif_path, normal_path, tumor_path, c(1,5,10))$signal_noise_df

  psn_expect_out <- readRDS(testthat::test_path("dummydata", "SignalNoise.rds"))
  testthat::expect_equal(psn_out, psn_expect_out)
})


