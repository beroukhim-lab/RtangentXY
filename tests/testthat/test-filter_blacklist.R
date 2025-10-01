test_that("filter_blacklist() works on pif filepath and small example", {
  pif_path <- testthat::test_path("test_data", "raw_inputs", "set_1", "probe_information.txt")
  filter_out <- filter_blacklist(pif_path)

  filter_expect_out <- readRDS(testthat::test_path("test_data", "pif_blacklist_filtered.rds"))
  expect_equal(filter_out, filter_expect_out)

  # Test that column order doesn't matter
  pif_path <- testthat::test_path("test_data", "raw_inputs", "set_2", "probe_annotation.txt") 
  filter_out <- filter_blacklist(pif_path) %>% dplyr::relocate('locus')
  expect_equal(filter_out, filter_expect_out)

  pif <- readr::read_tsv(pif_path, progress=FALSE, show_col_types=FALSE) %>%
    dplyr::slice(1:10)

  blacklist <- tibble::tibble(
    chr = c(1,   1,    1,    2,    "X",  "Y"),
    s =   c(2,   1500, 3650, 1230, 100,  400),
    e =   c(200, 2003, 4000, 2000, 2100, 505)
  )
  res <- filter_blacklist(pif, blacklist = blacklist, start_col = "s", end_col = "e")
  expect <- pif[c(2,3, 6:10), ]
  expect_equal(res, expect)

})

test_that("filter_blacklist() works on pif R object", {
  pif <- readRDS(testthat::test_path("test_data", "raw_inputs", "set_1", "probe_information.rds"))
  filter_out <- filter_blacklist(pif)

  filter_expect_out <- readRDS(testthat::test_path("test_data", "pif_blacklist_filtered.rds"))

  testthat::expect_equal(filter_out, filter_expect_out)
})
