### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: May 19, 2024
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Run CBS segmentation on normalized tumor signals
#'
#' @description
#' The function `run_cbs` expects the output of `runTangent()` or an input of
#' the same format. This function expects the locus column to be in the format
#' chr:start-end (e.g. 1:500-9501). There are more arguments in `DNAcopy::segment`
#' that aren't modified in this function and package. In this case, the `R/runCBS.R`
#' file can be modified to include these arguments.
#'
#' @param tnorm_df The `runTangent()` output tibble. It should be in the same format as the tumor signal matrix file
#' @param data_type This is the `data.type` argument for `DNAcopy::CNA()`. Default is logratio, but can be `'binary'` if LOH data
#' @param alpha This is the `alpha` argument for `DNAcopy::segment()`. Default is 0.005
#' @param min_width This is the `min.width` argument for `DNAcopy::segment()`. Default is 3
#'
#' @returns A named list with the outputs of the CBS algorithm. Each element is another list, returned by `DNAcopy::segment()`
#'
#' @import DNAcopy
#' @export

run_cbs <- function(tnorm_df, data_type = 'logratio', alpha = 0.005, min_width = 3) {
  # First transform this matrix so that we can parallelize: loci as column names and sample names as row names
  # Assumes locus is a column
  tnorm_data <- tnorm_df[ , -1, drop = FALSE]
  rownames(tnorm_data) <- tnorm_df[ , 1, drop = TRUE]
  tnorm_data <- t(tnorm_data)

  ret_list <- apply(tnorm_data, 1, row.cbs, data_type, alpha, min_width)
  return(ret_list)
}

# The following is a function called for each row in the dataframe
row.cbs <- function(r, data.type, alpha, min.width) {
  chr <- unlist(lapply(strsplit(names(r), ":"), `[[`, 1))
  start <- unlist(lapply(strsplit(names(r), ":"), `[[`, 2))
  start <- as.numeric(unlist(lapply(strsplit(start, "-"), `[[`, 1)))

  cna <- DNAcopy::CNA(r, chr, start, data.type = data.type)
  smoothed_cna <- DNAcopy::smooth.CNA(cna)
  segmented_cna <- DNAcopy::segment(smoothed_cna, verbose = 0, alpha = alpha, min.width = min.width)

  return(segmented_cna)
}
