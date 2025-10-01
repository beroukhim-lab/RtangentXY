test_that("linear_transformation works", {
  # Set 1
  sif_path <- testthat::test_path("test_data", "raw_inputs",  "set_1", "sample_information.txt")
  tsig_path <- testthat::test_path("test_data", "raw_inputs", "set_1", "tumor_log2RCN.txt")
  nsig_path <- testthat::test_path("test_data", "raw_inputs", "set_1", "normal_log2RCN.txt")

  # Load data
  sif <- readr::read_delim(sif_path, progress=FALSE, show_col_types=FALSE)
  n.df <- read_and_format_input(nsig_path, locus_to_rownames = TRUE)
  t.df <- read_and_format_input(tsig_path, locus_to_rownames = TRUE)

  res <- linear_transformation(sif, n.df, t.df)
  n_expect <- readr::read_delim(testthat::test_path(
    "test_data", "set_1",
    "LinearTransformation", "normal_log2RCN_linearTrans.txt"
  ), progress=FALSE, show_col_types=FALSE) %>% as.data.frame()
  t_expect <- readr::read_delim(testthat::test_path(
    "test_data", "set_1",
    "LinearTransformation", "tumor_log2RCN_linearTrans.txt"
  ), progress=FALSE, show_col_types=FALSE) %>% as.data.frame()
  expect_equal(res$n.df, n_expect)
  expect_equal(res$t.df, t_expect)


  # Set 2
  sif_path <- testthat::test_path("test_data", "raw_inputs",  "set_2", "sample_information.txt")
  tsig_path <- testthat::test_path("test_data", "raw_inputs", "set_2", "tumor_log2RCN.txt")
  nsig_path <- testthat::test_path("test_data", "raw_inputs", "set_2", "normal_log2RCN.txt")

  # Load data
  sif <- readr::read_delim(sif_path, progress=FALSE, show_col_types=FALSE)
  n.df <- read_and_format_input(nsig_path, locus_to_rownames = TRUE)
  t.df <- read_and_format_input(tsig_path, locus_to_rownames = TRUE)

  res <- linear_transformation(sif, n.df, t.df)
  n_expect <- readr::read_delim(testthat::test_path(
    "test_data", "set_2",
    "LinearTransformation", "normal_log2RCN_linearTrans.txt"
  ), progress=FALSE, show_col_types=FALSE) %>% as.data.frame()
  t_expect <- readr::read_delim(testthat::test_path(
    "test_data", "set_2",
    "LinearTransformation", "tumor_log2RCN_linearTrans.txt"
  ), progress=FALSE, show_col_types=FALSE) %>% as.data.frame()
  expect_equal(res$n.df, n_expect)
  expect_equal(res$t.df, t_expect)
})
