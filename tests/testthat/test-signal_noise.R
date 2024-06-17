test_that("signal_noise() correctly runs on filepaths", {
  pif_path <- testthat::test_path("dummydata", "raw_inputs", "probe_information.txt")
  tumor_path <- testthat::test_path("dummydata", "raw_inputs", "tumor_log2RCN.txt")
  signal_noise_out <- signal_noise(tumor_path, pif_path)
  signal_noise_expectout <- readRDS(testthat::test_path("dummydata", "tumor_signal_prenormalization_signal_noise_out.rds"))
  testthat::expect_equal(signal_noise_out, signal_noise_expectout)
})

test_that("signal_noise() correctly runs on R objects", {
  pif <- readRDS(testthat::test_path("dummydata", "raw_inputs", "probe_information.rds"))
  t.df <- readRDS(testthat::test_path("dummydata", "raw_inputs", "tumor_log2RCN.rds"))
  signal_noise_out <- signal_noise(t.df, pif)
  signal_noise_expectout <- readRDS(testthat::test_path("dummydata", "tumor_signal_prenormalization_signal_noise_out.rds"))
  testthat::expect_equal(signal_noise_out, signal_noise_expectout)
})
