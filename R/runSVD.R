### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: March 25, 2024
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Conduct SVD
#'
#' @description
#' Conduct singular value decomposition (SVD) on the normal signals in
#' preparation for input into the Tangent algorithm
#'
#' @param sif_filepath The filepath for the sample information file
#' @param ndf_filepath The filepath for the normal signal matrix file
#'
#' @returns A matrix with the male X chromosomes linearly transformed
#'
#' @import readr
#' @import dplyr
#' @import tibble
#' @export

run_svd <- function(sif_filepath, ndf_filepath) {
  cat('\nRunning SVD ...\n')
  sif <- readr::read_delim(sif_filepath, progress=FALSE, show_col_types=FALSE)
  n.lt.df <- transform_normals(sif_filepath, ndf_filepath) %>%
    tibble::column_to_rownames('locus')

  n.autox <- n.lt.df[!grepl('^Y', rownames(n.lt.df)),] %>%
    as.matrix()

  n.autox.svd <- svd(n.autox)
  colnames(n.autox.svd$u) <- colnames(n.autox)
  rownames(n.autox.svd$u) <- rownames(n.autox)
  rownames(n.autox.svd$v) <- colnames(n.autox)

  return(n.autox.svd)
}

