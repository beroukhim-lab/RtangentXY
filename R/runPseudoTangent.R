### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: April 8, 2024
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Runs the Pseudotangent pipeline
#'
#' @description
#' Pseudotangent also isolates the tumor signals using the Tangent algorithm,
#' but it should be used when the set of normals are particularly non-representative.
#' NOTE: all signal file rows should be in order (ie. from 1-22, X, Y).
#'
#' @param sif_filepath The filepath for the sample information file
#' @param ndf_filepath The filepath for the normal signal matrix file
#' @param tdf_filepath The filepath for the tumor signal matrix file
#' @param n_latent_init An integer representing the number of latent factors to reconstruct normal subspace
#' @param n_latent_tsvd An integer representing the number of latent factors for each of the runs
#'
#' @returns A normalized tumor signal matrix
#'
#' @import readr
#' @import dplyr
#' @import tibble
#' @import pracma
#' @importFrom stats median
#' @export

run_pseudotangent <- function(sif_filepath, ndf_filepath, tdf_filepath, n_latent_init, n_latent_tsvd) {
  # Step 1: Run TangentXY on a small set of normals
  cat('\nRunning Tangent on initial set of normals...\n')
  step1_tangent_out <- run_tangent(sif_filepath, ndf_filepath, tdf_filepath, n_latent_init)

  # Step 2: Run CBS on the TangentXY outputs to get the tentative CN profile
  cat('\nRunning CBS on inital Tangent output...\n')
  return(1)
}
