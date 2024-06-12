test_that("filter_blacklist() works on pif filepath", {
  pif_path <- testthat::test_path("dummydata", "raw_inputs", "probe_information.txt")
  filter_out <- filter_blacklist(pif_path, cores = 2)

  filter_expect_out <- readRDS(testthat::test_path("dummydata", "pif_blacklist_filtered.rds"))

  testthat::expect_equal(filter_out, filter_expect_out)
})

test_that("filter_blacklist() works on pif R object", {
  pif <- readRDS(testthat::test_path("dummydata", "raw_inputs", "probe_information.rds"))
  filter_out <- filter_blacklist(pif, cores = 2)

  filter_expect_out <- readRDS(testthat::test_path("dummydata", "pif_blacklist_filtered.rds"))

  testthat::expect_equal(filter_out, filter_expect_out)
})

test_that("filter_blacklist() works on signal filepath", {
  tumor_path <- testthat::test_path("dummydata", "raw_inputs", "tumor_log2RCN.txt")
  filter_out <- filter_blacklist(tumor_path, cores = 2)

  filter_expect_out <- readRDS(testthat::test_path("dummydata", "tdf_blacklist_filtered.rds"))

  testthat::expect_equal(filter_out, filter_expect_out)
})

test_that("filter_blacklist() works on signal R object", {
  t.df <- readRDS(testthat::test_path("dummydata", "raw_inputs", "tumor_log2RCN.rds"))
  filter_out <- filter_blacklist(t.df, cores = 2)

  filter_expect_out <- readRDS(testthat::test_path("dummydata", "tdf_blacklist_filtered.rds"))

  testthat::expect_equal(filter_out, filter_expect_out)
})
