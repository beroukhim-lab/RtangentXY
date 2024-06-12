### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: June 12, 2024
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Filter common germline CN variants
#'
#' @description
#' Pseudotangent also isolates the tumor signals using the Tangent algorithm,
#' but it should be used when the set of normals are particularly non-representative.
#' NOTE: all signal file rows should be in order (ie. from 1-22, X, Y).
#'
#' @param df The tibble of or the filepath for the file to filter. Must contain a `locus`
#' column of the following form: `"{chr}:{start}-{end}"` where `{chr}` is 1-22 or X,
#' and where `{start}` and `{end}` are all integers
#' @param version By default `hg38`. If the data in `df` is aligned to `hg19`, set this parameter to
#' `"hg19"`. A dataframe can also be supplied as a parameter for custom filtering. Must contain a `chr`
#' column and columns with start and end CNV loci
#' @param col_range By default `NULL`, as the flanking regions of the blacklist files are selected.
#' If a custom version dataframe is supplied (or if the flanking regions are not desired), a numeric
#' vector can be supplied in this field of the form `c({start}, {end})`, where `{start}` and `{end}`
#' are integers corresponding to the desired start and end columns in `version`
#' @param cores The number of cores to parallelize with
#'
#' @returns A normalized tumor signal matrix
#'
#' @import readr
#' @importFrom GenomicRanges GRanges countOverlaps
#' @import parallel
#' @export

filter_blacklist <- function(df, version = "hg38", col_range = NULL, cores = 2) {

  # Load data
  if (inherits(df, "character")) {
    df <- readr::read_delim(df, progress=FALSE, show_col_types=FALSE)
  } else { df <- df }

  # Load in the blacklist file
  if (version == "hg38") { blacklist_df <- RtangentXY::hg38_blacklist }
  else if (version == "hg19") { blacklist_df <- RtangentXY::hg19_blacklist }
  else if (is.data.frame(version)) { blacklist_df <- version }
  else { stop("ERROR: version parameter not recognized") }

  # Make sure dataframe as a locus column
  if (!("locus" %in% colnames(df))) { stop("ERROR: input dataframe does not have a locus column") }

  # Create GRanges object from blacklist table
  if (is.null(col_range)) {
    start <- as.numeric(as.vector(unlist(blacklist_df[, 5, drop = FALSE])))
    end <- as.numeric(as.vector(unlist(blacklist_df[, 6, drop = FALSE])))
  } else if (length(col_range) == 2) {
    start <- as.numeric(as.vector(unlist(blacklist_df[, col_range[1], drop = FALSE])))
    end <- as.numeric(as.vector(unlist(blacklist_df[, col_range[2], drop = FALSE])))
  } else { stop("ERROR: col_range parameter not recognized") }

  lox_ranges <- paste0(paste0("chr", blacklist_df$chr), ":", as.character(start), "-", as.character(end))
  blacklist_gr <- GenomicRanges::GRanges(lox_ranges)

  # Find whether each row overlaps
  nonoverlapping_rows <- parallel::mclapply(df$locus, is.not.overlapping, blacklist_gr, mc.cores = cores)
  nonoverlapping_rows <- unlist(nonoverlapping_rows)

  filtered_df <- df[nonoverlapping_rows, ]

  return(filtered_df)
}

is.not.overlapping <- function(row, blacklist_gr) {
  # Input: a row from a dataframe
  # Returns TRUE if this region overlaps
  suppressWarnings({
    row_gr <- GenomicRanges::GRanges(paste0("chr", row))
    overlaps <- GenomicRanges::countOverlaps(row_gr, blacklist_gr)
  })
  if (overlaps == 0) { return(TRUE) }
  else { return(FALSE) }
}
