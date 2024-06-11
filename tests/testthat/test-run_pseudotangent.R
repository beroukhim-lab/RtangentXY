test_that("Running Pseudotangent with filepath inputs works", {
  sif_path <- testthat::test_path("dummydata", "raw_inputs", "sample_information.txt")
  normal_path <- testthat::test_path("dummydata", "raw_inputs", "normal_log2RCN.txt")
  tumor_path <- testthat::test_path("dummydata", "raw_inputs", "tumor_log2RCN.txt")

  pseudotangent_expect_out <- as.data.frame(readRDS(testthat::test_path("dummydata", "PseudoTangent_test_output.rds")))
  pseudotangent_out <- as.data.frame(run_pseudotangent(sif_path, normal_path, tumor_path, 5, num_partition = 2, n_latent_part = 2))

  testthat::expect_equal(pseudotangent_expect_out, pseudotangent_out)
})

test_that("Running Pseudotangent with dataframe inputs works", {
  sif <- readRDS(testthat::test_path("dummydata", "raw_inputs", "sample_information.rds"))
  n.df <- readRDS(testthat::test_path("dummydata", "raw_inputs", "normal_log2RCN.rds"))
  t.df <- readRDS(testthat::test_path("dummydata", "raw_inputs", "tumor_log2RCN.rds"))

  pseudotangent_expect_out <- as.data.frame(readRDS(testthat::test_path("dummydata", "PseudoTangent_test_output.rds")))
  pseudotangent_out <- as.data.frame(run_pseudotangent(sif, n.df, t.df, 5, num_partition = 2, n_latent_part = 2))

  testthat::expect_equal(pseudotangent_expect_out, pseudotangent_out)
})
