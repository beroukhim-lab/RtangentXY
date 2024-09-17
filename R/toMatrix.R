### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: August 5, 2024
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Convert Read Counts to Matrix
#'
#' @description
#' This takes in the filepath or tibble of GATK DepthOfCoverage or CollectReadCounts
#' outputs, checks if they contain the same loci, and preprocesses them for the
#' Tangent pipeline. Input either all tumor signals or normal signals to be used
#' for the analysis.
#'
#' @param ... The tibble of or the filepath for the read counts
#'
#' @returns A pre-processed read counts matrix
#'
#' @import readr
#' @export

to_matrix <- function(...) {
  # Use this to check if all the loci are the same
  loci_col <- NULL

  for (i in list(...)) {
    # Check if it's a filepath or dataframe
    if (inherits(i, "character")) {
      cur_df <- readr::read_delim(i, progress=FALSE, show_col_types=FALSE)
    } else { cur_df <- i }

  }
}
