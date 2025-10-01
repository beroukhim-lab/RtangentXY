### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: September 2025
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Run CBS segmentation on normalized tumor signals
#'
#' @param tnorm_df Normalized tumor signal from \code{\link{run_tangent}}.
#' @param alpha Passed to \code{\link[DNAcopy]{segment}}.
#' @param min_width Passed to \code{\link[DNAcopy]{segment}}.
#' @param n_cores Number of cores to use for parallel processing.
#'
#' @returns A named list with the outputs of the CBS algorithm. Each element corresponds to a sample.
#'
#' @export
run_cbs <- function(tnorm_df, alpha = 0.005, min_width = 3, n_cores = 1) {
  tnorm_df <- tnorm_df %>% tibble::column_to_rownames('locus')

  res <- parallel::mclapply(colnames(tnorm_df), function(x) {
    sample_data <- setNames(tnorm_df[, x], rownames(tnorm_df))
    sample.cbs(sample_data, alpha, min_width)
  }, mc.cores = n_cores)
  names(res) <- colnames(tnorm_df)

  return(res)
}

#' Run CBS segmentation on a single sample
#' 
#' Uses the start position as the probe location.
#' @keywords internal
sample.cbs <- function(sample, alpha, min.width) {
  chr <- unlist(lapply(strsplit(names(sample), ":"), `[[`, 1))
  start <- unlist(lapply(strsplit(names(sample), ":"), `[[`, 2))
  start <- as.numeric(unlist(lapply(strsplit(start, "-"), `[[`, 1)))

  cna <- DNAcopy::CNA(sample, chr, start, data.type = 'logratio')
  smoothed_cna <- DNAcopy::smooth.CNA(cna)
  segmented_cna <- DNAcopy::segment(smoothed_cna, verbose = 0, alpha = alpha, min.width = min.width)

  return(segmented_cna)
}
