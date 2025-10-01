### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: September 2025
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Calculate SVD of linearly transformed normal signal matrix
#'
#' Calculate singular value decomposition (SVD) of the linearly transformed normal signal matrix
#' in preparation for input into the Tangent algorithm.
#'
#' @param nt.df A data frame containing the linearly transformed normal signal matrix (output of \code{\link{linear_transformation}})
#'
#' @returns A list containing the SVD components
# 
# @examples
# nt.df <- linear_transformation(example_sif, example_nsig_df, example_tsig_df)$n.df
# res <- run_svd(nt.df)
# str(res)
#'
#' @keywords internal
run_svd <- function(nt.df) {

  cat('\nRunning SVD ...\n')

  n.lt.df <- nt.df %>%
    tibble::column_to_rownames('locus')

  n.autox <- n.lt.df[!grepl('^Y', rownames(n.lt.df)), , drop = FALSE] %>%
    as.matrix()

  n.autox.svd <- svd(n.autox)
  colnames(n.autox.svd$u) <- colnames(n.autox)
  rownames(n.autox.svd$u) <- rownames(n.autox)
  rownames(n.autox.svd$v) <- colnames(n.autox)

  return(n.autox.svd)
}

