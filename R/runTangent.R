### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: March 26, 2024
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' This code runs the TangentXY algorithm, which isolates the tumor signal using a set of normals.
#'
#' @param sif_filepath The filepath for the sample information file
#' @param ndf_filepath The filepath for the normal signal matrix file
#' @param tdf_filepath The filepath for the tumor signal matrix file
#' @param n_latent An integer representing the number of latent factors to reconstruct normal subspace
#'
#' @return A normalized tumor signal matrix
#'
#' @import readr
#' @import dplyr
#' @export

run_tangent <- function(sif_filepath, ndf_filepath, tdf_filepath, n_latent) {
  ndf_lt <- transform_normals(sif_filepath, ndf_filepath)
  n_svd <- run_svd(sif_filepath, ndf_filepath)

  return(1)
}

