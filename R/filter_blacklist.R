### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: September 2025
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Filter out common germline copy-number variants
#'
#' @param probe_info Tibble or filepath to a text file containing probe information.
#' Should have a "locus" column of the following form: `"{chr}:{start}-{end}"` where `{chr}` is 1-22, X, or Y
#' and where `{start}` and `{end}` are genomic coordinates
#' @param blacklist By default `hg38`. If the data in `probe_info` is aligned to `hg19`, set this parameter to
#' `"hg19"`. A dataframe or filepath can also be supplied as a parameter for custom filtering.
#' Blacklist should have a `chr` column with chromosome names (1-22, X, or Y) and columns for start and end coordinates.
#' @param start_col The name of the column in the blacklist file to use as the start coordinate for filtering
#' @param end_col The name of the column in the blacklist file to use as the end coordinate for filtering
#'
#' @returns Input `probe_info` dataframe with rows overlapping the blacklist removed
#'
#' @export
filter_blacklist <- function(probe_info,
                            blacklist = "hg38",
                            start_col = "flanking.start",
                            end_col = "flanking.end") {

  # Load data
  if (inherits(probe_info, "character")) {
    probe_info <- readr::read_delim(probe_info, progress=FALSE, show_col_types=FALSE)
  }

  # Load blacklist file
  if (!inherits(blacklist, "data.frame")) {
    if (blacklist == "hg38") { blacklist <- RtangentXY::hg38_blacklist }
    else if (blacklist == "hg19") { blacklist <- RtangentXY::hg19_blacklist }
    else if (inherits(blacklist, "character")) { blacklist <- readr::read_delim(blacklist, progress=FALSE, show_col_types=FALSE) }
    else { stop("ERROR: blacklist parameter not recognized") }
  }

  # Make sure probe information table has a locus column
  if (!("locus" %in% colnames(probe_info))) { stop("ERROR: input dataframe does not have a locus column") }

  blacklist_gr <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(paste0("chr", blacklist$chr)),
    ranges = IRanges::IRanges(start = blacklist[[start_col]], end = blacklist[[end_col]])
  )

  probe_gr <- GenomicRanges::GRanges(paste0("chr", probe_info$locus))
  overlaps <- GenomicRanges::findOverlaps(probe_gr, blacklist_gr)
  remove_idx <- unique(S4Vectors::queryHits(overlaps))

  res <- probe_info[-remove_idx, ]

  return(res)
}