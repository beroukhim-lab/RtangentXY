### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: September 2025
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Converts CBS output to signal format
#'
#' @description
#' Converts output of \code{\link{run_cbs}} into the format of the
#' signal files (normal and tumor) used elsewhere in this package. Fills female Y chromosome values with NAs.
#'
#' @inheritParams run_tangent
#' @param cbs_out The output from \code{\link{run_cbs}} (i.e. list of outputs from \code{\link[DNAcopy]{segment}}).
#'
#' @returns Tumor signal matrix in the same format as the input tumor signal matrix.
#'
#' @keywords internal
convert_cbs <- function(cbs_out, sif_df, tsig_df) {
  # Load data
  # We're going to use tumor signal file format to create the output cbs
  if (inherits(sif_df, "character")) {
    sif <- readr::read_delim(sif_df, progress=FALSE, show_col_types=FALSE)
  } else { sif <- sif_df }
  tumor_signal_data <- read_and_format_input(tsig_df, locus_to_rownames = TRUE)

  cbs_df_list <- list()

  for (samplename in colnames(tumor_signal_data)) {
    cur_sample_out <- cbs_out[[samplename]]

    # t is the table containing the actual output from the total cbs output for the current sample
    t <- cur_sample_out$output %>%
      dplyr::mutate(chrom = factor(.data$chrom, levels = c(as.character(1:22), 'X', 'Y')))
    t <- t[with(t, order(chrom, loc.start)),]

    # The vector containing the segment means for the current sample
    seg_mean_vector <- unlist(mapply(rep, t$seg.mean, t$num.mark))

    # Fill missing values with NAs (e.g. female samples don't have Y values)
    if (sif[sif$sample.id == samplename, ]$sex == "female") {
      diff <- nrow(cur_sample_out$data) - length(seg_mean_vector)
      na_vector <- rep(NA, diff)
      seg_mean_vector <- c(seg_mean_vector, na_vector)
    }

    cbs_df_list[[samplename]] <- seg_mean_vector
  }

  final_cbs_df <- do.call(cbind, cbs_df_list)

  # All signal files rows should be ordered (1-22, X, Y)
  rownames(final_cbs_df) <- rownames(tumor_signal_data)
  return (final_cbs_df)
}
