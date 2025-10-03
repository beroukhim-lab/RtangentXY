test_that("run_svd() correctly performs SVD on input data", {
  # Set 1
  sif_path <- testthat::test_path("test_data", "raw_inputs",  "set_1", "sample_information.txt")
  tsig_path <- testthat::test_path("test_data", "raw_inputs", "set_1", "tumor_log2RCN.txt")
  nsig_path <- testthat::test_path("test_data", "raw_inputs", "set_1", "normal_log2RCN.txt")

  # Load data
  sif <- readr::read_delim(sif_path, progress=FALSE, show_col_types=FALSE)
  n.df <- read_and_format_input(nsig_path, locus_to_rownames = TRUE)
  t.df <- read_and_format_input(tsig_path, locus_to_rownames = TRUE)

  lt_res <- linear_transformation(sif, n.df, t.df)
  res <- run_svd(lt_res$n.df)
  svd_expect <- readRDS(testthat::test_path("test_data", "set_1", "SVD", "n.autox.svd.rds"))

  expect_equal(res$d, svd_expect$d)

  # Set 2
  sif_path <- testthat::test_path("test_data", "raw_inputs",  "set_2", "sample_information.txt")
  tsig_path <- testthat::test_path("test_data", "raw_inputs", "set_2", "tumor_log2RCN.txt")
  nsig_path <- testthat::test_path("test_data", "raw_inputs", "set_2", "normal_log2RCN.txt")

  # Load data
  sif <- readr::read_delim(sif_path, progress=FALSE, show_col_types=FALSE)
  n.df <- read_and_format_input(nsig_path, locus_to_rownames = TRUE)
  t.df <- read_and_format_input(tsig_path, locus_to_rownames = TRUE)

  lt_res <- linear_transformation(sif, n.df, t.df)
  res <- run_svd(lt_res$n.df)
  svd_expect <- readRDS(testthat::test_path("test_data", "set_2", "SVD", "n.autox.svd.rds"))

  expect_equal(res$d, svd_expect$d)
})
