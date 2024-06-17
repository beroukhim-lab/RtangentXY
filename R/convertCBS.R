### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: May 21, 2024
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Converts CBS output to signal format
#'
#' @description
#' This takes the output from the `runCBS()` function in this package (ie. the
#' output in the DNAcopy pakage format) and converts it into the format of the
#' signal files (normal and tumor) used elsewhere in this package. This makes
#' computation between the CBS outputs and signal files simple.
#'
#' @param cbs_out The output from `runCBS()` (i.e. list of outputs from `DNAcopy::segment()`)
#' @param sif_df The tibble of or the filepath for the sample information file
#' @param tsig_df The tibble of or the filepath for the tumor signal matrix file
#'
#' @returns A normalized tumor signal matrix
#'
#' @import readr
#' @import dplyr
#' @import tibble
#' @export

convert_cbs <- function(cbs_out, sif_df, tsig_df) {

  # Note: when the CBS was conducted, we only took into account the start (location)
  # Load data
  # We're going to use tumor signal file format to create the output cbs
  if (inherits(sif_df, "character")) {
    sif <- readr::read_delim(sif_df, progress=FALSE, show_col_types=FALSE)
  } else { sif <- sif_df }

  if (inherits(tsig_df, "character")) {
    tumor_signal_data <- readr::read_delim(tsig_df, progress=FALSE, show_col_types=FALSE) %>%
      tibble::column_to_rownames('locus')
  } else {
    if ('locus' %in% colnames(tsig_df)) { tumor_signal_data <- tsig_df %>% tibble::column_to_rownames('locus') }
    else { tumor_signal_data <- tsig_df }
  }

  cbs_df_list <- list()

  for (samplename in colnames(tumor_signal_data)) {
    cur_sample_out <- cbs_out[[samplename]]

    # t is the table containing the actual output from the total cbs output for the current sample
    t <- cur_sample_out$output %>%
      dplyr::mutate(chrom = factor(chrom, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                     "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                                                     "21", "22", "X", "Y")))
    t <- t[with(t, order(chrom, loc.start)),]

    # The vector containing the segment means for the current sample
    seg_mean_vector <- unlist(mapply(rep, t$seg.mean, t$num.mark))

    if (sif[sif$sample.id == samplename, ]$gender == "female") {
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
