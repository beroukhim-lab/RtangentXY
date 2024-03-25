### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: March 25, 2024
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Conduct singular value decomposition (SVD) in preparation for input into the Tangent algorithm
#'
#' @param sif_filepath The filepath for the sample information file
#' @param nltdf_filepath The filepath for the normal signal matrix file that has been linearly transformed by transform_normals()
#'
#' @return A matrix with the male X chromosomes linearly transformed
#'
#' @import readr
#' @import dplyr
#' @export

run_svd <- function(sif_filepath, nltdf_filepath) {
  sif <- readr::read_delim(sif_filepath, progress=FALSE, show_col_types=FALSE)
  n.lt.df <- readr::read_delim(nltdf_filepath, progress=FALSE, show_col_types=FALSE) %>%
    column_to_rownames('locus')

  n.autox <- n.lt.df[!grepl('^Y', rownames(n.lt.df)),] %>%
    as.matrix()

  n.autox.svd <- svd(n.autox)
  colnames(n.autox.svd$u) <- colnames(n.autox)
  rownames(n.autox.svd$u) <- rownames(n.autox)
  rownames(n.autox.svd$v) <- colnames(n.autox)

  return(n.autox.svd)
}

