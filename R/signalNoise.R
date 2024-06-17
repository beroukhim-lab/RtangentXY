### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: June 17, 2024
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Calculate the signal and noise profile of each sample
#'
#' @description
#' This function calculates the signal and noise profile of each sample. In this
#' case, the signal is defined as the standard deviation of the median signal of
#' each arm (grouped into autosomes, X, and Y chromosomes), and the noise is
#' defined as the median of the absolute values of the vector of consecutive
#' differences between each signal.
#'
#' @param signal_df The tibble of or the filepath for the tumor signal file (pre- or post-normalization)
#' @param pif_df The tibble of or the filepath for the probe information file
#'
#' @returns A dataframe of the signal and noise profiles for each sample
#'
#' @import readr
#' @import tibble
#' @import dplyr
#' @import parallel
#' @importFrom stats setNames median sd
#' @export

signal_noise <- function(signal_df, pif_df) {

  # Load data
  if (inherits(signal_df, "character")) {
    sig.df <- readr::read_delim(signal_df, progress=FALSE, show_col_types=FALSE) %>%
      tibble::column_to_rownames('locus')
  } else {
    if ('locus' %in% colnames(signal_df)) { sig.df <- signal_df %>% tibble::column_to_rownames('locus') }
    else { sig.df <- signal_df }
  }

  if (inherits(pif_df, "character")) {
    pif <- readr::read_delim(pif_df, progress=FALSE, show_col_types=FALSE)
  } else { pif <- pif_df }

  #signal.noise.list <- parallel::mclapply(sig.df, calc.signal.noise, pif, mc.cores = 2)
  signal.noise.list <- lapply(sig.df, calc.signal.noise, pif)
  signal.noise.df <- signal.noise.list %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(sample.id=names(signal.noise.list))
  return(signal.noise.df)
}

calc.signal.noise <- function(list, pif) {
  data <- list %>%
    as.data.frame() %>%
    stats::setNames('signal') %>%
    dplyr::bind_cols(pif)

  signal.auto <- data %>%
    dplyr::filter(!chr %in% c('X', 'Y')) %>%
    dplyr::group_by(chr, arm) %>%
    dplyr::summarize(arm.median=median(signal)) %>%
    dplyr::ungroup() %>%
    dplyr::summarize(signal=stats::sd(arm.median)) %>%
    dplyr::pull(signal)

  noise.auto <- data %>%
    dplyr::filter(!chr %in% c('X', 'Y')) %>%
    dplyr::pull(signal) %>%
    diff() %>%
    abs() %>%
    stats::median()

  noise.x <- data %>%
    dplyr::filter(chr=='X') %>%
    dplyr::pull(signal) %>%
    diff() %>%
    abs() %>%
    stats::median()

  noise.y <- data %>%
    dplyr::filter(chr=='Y') %>%
    dplyr::pull(signal) %>%
    diff() %>%
    abs() %>%
    stats::median()

  result.df <- data.frame(signal.auto=signal.auto,
                          noise.auto=noise.auto,
                          sn.auto=signal.auto/noise.auto,
                          noise.chrx=noise.x,
                          noise.chry=noise.y)

  return(result.df)
}
