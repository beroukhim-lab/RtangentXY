### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: September 18, 2024
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Save segmentation file
#'
#' @description
#' Save the output of CBS to a given filepath as a segmentation file.
#'
#' @param seg_list The output of `run_cbs()`, i.e. a list of DNAcopy objects
#' @param dirpath The path of the directory where the segmentation files should be saved
#'
#' @return 1 if the file was successfully saved
#'
#' @import DNAcopy
#' @import tibble
#' @import readr
#' @export

save_segmentation <- function(seg_list, dirpath) {
  for (i in 1:length(seg_list)) {
    cur_sample_seg <- tibble::as_tibble(seg_list[[i]]$output)
    cur_sample <- names(seg_list)[i]

    cur_filepath <- file.path(dirpath, paste0(cur_sample, ".seg.text"))
    readr::write_delim(cur_sample_seg, file = cur_filepath, delim = "\t", quote = "none")
  }
  return(1)
}
