### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: September 2025
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Save segmentation file
#'
#' @description
#' Save each sample's segmentation file as a tab-delimited text file.
#'
#' @param seg_list The output of \code{\link{run_cbs}}, i.e. a list of DNAcopy objects
#' @param output_dir Directory where the segmentation files should be saved
#' @param filenames A character vector of the names to save each sample by. This should be the same length and order as `seg_list`
#'
#' @returns 1 if the file was successfully saved
#'
#' @export
save_segmentation <- function(seg_list, output_dir, filenames = NULL) {
  if (is.null(filenames)) {
    for (i in 1:length(seg_list)) {
      cur_sample_seg <- tibble::as_tibble(seg_list[[i]]$output)
      cur_sample <- names(seg_list)[i]

      cur_filepath <- file.path(output_dir, paste0(cur_sample, ".seg.txt"))
      readr::write_delim(cur_sample_seg, file = cur_filepath, delim = "\t", quote = "none")
    }
  } else {
    for (i in 1:length(seg_list)) {
      if (length(seg_list) != length(filenames)) {
        stop("ERROR: seg_list and filenames must have the same length")
      }
      cur_sample_seg <- tibble::as_tibble(seg_list[[i]]$output)
      cur_sample <- filenames[i]

      cur_filepath <- file.path(output_dir, paste0(cur_sample, ".seg.txt"))
      readr::write_delim(cur_sample_seg, file = cur_filepath, delim = "\t", quote = "none")
    }
  }

  return(1)
}
