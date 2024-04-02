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
#' @returns A normalized tumor signal matrix
#'
#' @import readr
#' @import dplyr
#' @import tibble
#' @import pracma
#' @importFrom stats median
#' @export

run_tangent <- function(sif_filepath, ndf_filepath, tdf_filepath, n_latent) {
  sif <- readr::read_delim(sif_filepath, progress=FALSE, show_col_types=FALSE)
  t.df <- readr::read_delim(tdf_filepath, progress=FALSE, show_col_types=FALSE) %>%
    tibble::column_to_rownames('locus')

  # Other functions from this package
  ndf_lt <- transform_normals(sif_filepath, ndf_filepath)
  n.autox.svd <- run_svd(sif_filepath, ndf_filepath)

  ## Tangent on autosomes and chrX
  cat('\nRunning Tangent on autosomes and chrX ...\n')

  ## Reconstruct a normal subspace with a specific number of latent factors
  num.lf <- n_latent
  if (num.lf == 0 | num.lf > length(n.autox.svd$d)) {
    stop('-l option need to be greater than 0, and less than or equal to the number of normal samples.')
  } else if (num.lf == 1) {
    N.autox <- as.matrix(n.autox.svd$u[, num.lf]) %*% n.autox.svd$d[num.lf] %*% t(n.autox.svd$v[, num.lf])
  } else {
    N.autox <- n.autox.svd$u[, 1:num.lf] %*% diag(n.autox.svd$d[1:num.lf]) %*% t(n.autox.svd$v[, 1:num.lf])
  }

  T.autox <- t.df[!grepl('^Y', rownames(t.df)),] %>%
    as.matrix()

  ## Get the origin in the normal subspace
  N.autox.means <- apply(N.autox, 1, mean)
  N.autox0 <- N.autox - N.autox.means
  T.autox0 <- T.autox - N.autox.means

  ## Run Tangent
  Npi.autox <- pracma::pinv(N.autox0)
  weights.autox <- Npi.autox %*% T.autox0
  proj.autox <- N.autox0 %*% weights.autox
  T.autox.norm <- T.autox0 - proj.autox

  ## Re-scaling after Tangent by median
  T.autox.norm.medians <- T.autox.norm[!grepl('X', rownames(T.autox.norm)),] %>%
    apply(., 2, median)

  T.autox.norm.rescaled <- t(t(T.autox.norm)- T.autox.norm.medians)
  cat('Done.\n')

  ## Tangent on male chrY
  cat('\nRunning Tangent on male chrY ...\n')
  n.df <- readr::read_delim(ndf_filepath, progress=FALSE, show_col_types=FALSE) %>%
    tibble::column_to_rownames('locus')

  male.normals <- sif %>%
    dplyr::filter(sample.id %in% colnames(n.df) & gender=='male') %>%
    dplyr::pull(sample.id)

  male.tumors <- sif %>%
    dplyr::filter(sample.id %in% colnames(t.df) & gender=='male') %>%
    dplyr::pull(sample.id)

  N.m <- n.df[, male.normals] %>%
    as.matrix()

  T.m <- t.df[, male.tumors] %>%
    as.matrix()

  ## Get the origin in the normal subspace
  N.m.means <- apply(N.m, 1, mean)
  N.m0 <- N.m - N.m.means
  T.m0 <- T.m - N.m.means

  ## Run Tangent
  Npi.m <- pracma::pinv(N.m0)
  weights.m <- Npi.m %*% T.m0
  proj.m <- N.m0 %*% weights.m
  T.m.norm <- T.m0 - proj.m

  ## Re-scaling after Tangent by median
  T.m.norm.medians <- T.m.norm[!grepl('^X|^Y', rownames(T.m.norm)),] %>%
    apply(., 2, median)

  T.m.norm.rescaled <- t(t(T.m.norm)- T.m.norm.medians)
  T.m.y <- T.m.norm.rescaled[grepl('^Y', rownames(T.m.norm.rescaled)), ]

  ## Adjust chrY so that it is relative to CN=2
  T.m.y.adj <- T.m.y - 1
  cat('Done.\n')

  ## Combine "autosomes & chrX" and "chrY"
  ## Replace NAs in female chrY for downstream analysis (e.g. Circular Binary Segmentation)
  T.norm <- as.data.frame(T.autox.norm.rescaled) %>%
    dplyr::bind_rows(as.data.frame(T.m.y.adj))

  rownames(T.norm) <- rownames(t.df)

  T.norm <- T.norm %>%
    as.data.frame() %>%
    tibble::rownames_to_column('locus')

  return(T.norm)
}
