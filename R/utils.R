### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: September 2025
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

# We need to add a global variable for ".", as this is used in dplyr
utils::globalVariables(c("."))

#' Read in a data frame or filepath with a `locus` column and arrange it in 
#' order of increasing genomic position (1-22, X, Y; within each chromosome, by start position).
#'
#' @param input A tibble or a filepath to a text file containing the data. Must have a `locus` column
#' in the format `"{chr}:{start}-{end}"` where `{chr}` is 1-22, X, or Y and where `{start}` and `{end}` are genomic
#' coordinates.
#' @param locus_to_rownames If `TRUE`, the `locus` column will be converted to rownames.
#'
#' @returns Input data frame arranged by genomic position
#' 
#' @keywords internal
read_and_format_input <- function(input, locus_to_rownames = FALSE) {
  # If input is a filepath, read it in
  if (inherits(input, "character")) {
    input <- readr::read_tsv(input, progress=FALSE, show_col_types=FALSE)
  }

  res <- arrange_by_genomic_position(input)

  if (locus_to_rownames) res <- res %>% tibble::column_to_rownames('locus')

  return(res)
}

#' Arrange a data frame by genomic position based on a `locus` column
#' @param df A tibble or data.frame with a `locus` column in the format `"{chr}:{start}-{end}"`
#' where `{chr}` is 1-22, X, or Y and where `{start}` and `{end}` are genomic coordinates.
#' @returns Input data frame arranged by genomic position
#' @keywords internal
arrange_by_genomic_position <- function(df) {
  df <- df %>%
    tidyr::separate_wider_delim(cols="locus", names=c('tmp_chr', 'tmp_pos'), delim=':', cols_remove = FALSE) %>%
    tidyr::separate_wider_delim(cols="tmp_pos", names=c('tmp_start', 'tmp_end'), delim='-', cols_remove = TRUE) %>%
    dplyr::mutate(tmp_chr = factor(.data$tmp_chr, levels = c(as.character(1:22), 'X', 'Y'))) %>%
    dplyr::arrange(.data$tmp_chr, as.numeric(.data$tmp_start)) %>%
    dplyr::select(-c('tmp_chr', 'tmp_start', 'tmp_end'))
  return(df)
}