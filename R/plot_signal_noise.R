### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: September 2025
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Plot signal and noise before and after TangentXY normalization
#'
#' @description
#' This function calculates and plots the signal and noise profiles of tumor samples before and
#' after TangentXY normalization with varying numbers of latent factors. The signal is defined
#' as the standard deviation of the median signal of each chromosome arm (grouped into autosomes, X, and Y chromosomes),
#' and the noise is defined as the median of the absolute values of the vector of
#' consecutive differences between each signal.
#'
#' @param tnorm A normalized tumor signal matrix output by \link{run_tangent}, or a list of such matrices.
#' @inheritParams run_tangent
#' @param pif_df Tibble or filepath to a text file containing probe information.
#' Should have a 'chr' column with chromosome names (1-22, X, Y) and an 'arm' column with arm names (p or q).
#' Should also have a 'locus' column with locus names in the format `"{chr}:{start}-{end}"` where `{chr}` is 1-22, X, or Y
#' and where `{start}` and `{end}` are genomic coordinates.
#' @param n_latent A numeric vector of the numbers of latent factors used in TangentXY normalization.
#' Should correspond to the elements in `tnorm`.
#' @param output_dir Directory to save the plot. If `NULL`, the plot will be printed to the screen.
#' @param n_cores Number of cores to use for parallel processing.
#' 
#' @examples 
#' n_latent <- c(5, 10)
#' tangent_res <- lapply(n_latent, function(nlf) {
#'  run_tangent(example_sif, example_nsig_df, example_tsig_df, nlf, make_plots = FALSE)
#' })
#' res <- plot_signal_noise(tangent_res, example_sif, example_pif,
#'                          example_tsig_df, n_latent = n_latent)
#'
#' @returns (Invisibly) A tibble containing the signal and noise profiles for each sample before and after normalization.
#'
#' @export
plot_signal_noise <- function(tnorm, sif_df, pif_df, tsig_df, n_latent, output_dir = NULL, n_cores = 1) {
  # Ensure tnorm is a list
  if (!inherits(tnorm, "list")) {
    tnorm <- list(tnorm)
  }
  # Load data
  if (inherits(sif_df, "character")) {
    sif <- readr::read_delim(sif_df, progress=FALSE, show_col_types=FALSE)
  } else { sif <- sif_df }
  # Read in signal data and arrange by genomic position
  t.df <- read_and_format_input(tsig_df, locus_to_rownames = TRUE)
  pif <- read_and_format_input(pif_df, locus_to_rownames = FALSE)

  # Associate each tangent output with its corresponding number of latent factors
  names(tnorm) <- paste("nlf", n_latent, sep="_")

  # Sort latent factors in increasing order for plotting
  num.lf <- sort(n_latent)

  options(dplyr.summarise.inform = FALSE)

  ## Pre-normalization signal
  cat('Calculating signal and noise in pre-normalization data...\n')

  t0.signal.noise.list <- parallel::mclapply(t.df, calc.signal.noise, pif = pif, mc.cores= n_cores)
  # t0.signal.noise.list <- lapply(t.df, calc.signal.noise)
  t0.signal.noise.df <- t0.signal.noise.list %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(sample.id=names(t0.signal.noise.list)) %>%
    dplyr::mutate(lf='Pre-normalization')

  ## After TangentXY
  for (i in 1:length(num.lf)) {
    num.lf.i <- num.lf[i]
    if (num.lf.i==1) {
      cat(paste0('Calculating signal and noise in normalized data with ', num.lf.i, ' latent factor...\n'))
    } else {
      cat(paste0('Calculating signal and noise in normalized data with ', num.lf.i, ' latent factors...\n'))
    }

    t.df.i <- tnorm[[paste("nlf", num.lf.i, sep="_")]] %>%
      tibble::column_to_rownames('locus')

    t.i.signal.noise.list <- parallel::mclapply(t.df.i, calc.signal.noise, pif = pif, mc.cores=n_cores)
    t.i.signal.noise.df <- t.i.signal.noise.list %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(sample.id=names(t.i.signal.noise.list)) %>%
      dplyr::mutate(lf=as.character(num.lf.i))

    if (i==1) {
      signal.noise.df <- t0.signal.noise.df %>% dplyr::bind_rows(t.i.signal.noise.df)
    } else {
      signal.noise.df <- signal.noise.df %>% dplyr::bind_rows(t.i.signal.noise.df)
    }
  }

  ## Plot
  signal.noise.df.l <- signal.noise.df %>%
    tidyr::pivot_longer(cols=dplyr::matches('^signal\\.|^noise\\.|^sn\\.'), names_pattern='(^signal\\.|^noise\\.|^sn\\.)(.*)$', names_to=c('metric', 'chr')) %>%
    dplyr::mutate(lf=dplyr::case_when(.data$chr=='chry' & .data$lf!='Pre-normalization' ~ 'TangentXY', TRUE ~ .data$lf)) %>%
    dplyr::distinct(.data$sample.id, .data$lf, .data$metric, .data$chr, .keep_all=TRUE) %>%
    dplyr::mutate(lf=factor(.$lf, levels=c('Pre-normalization', num.lf, 'TangentXY'))) %>%
    dplyr::mutate(metric=dplyr::case_when(.data$metric=='signal.' ~ 'Signal',
                            .data$metric=='noise.' ~ 'Noise',
                            .data$metric=='sn.' ~ 'SN')) %>%
    dplyr::mutate(metric=factor(.$metric, levels=c('Signal', 'Noise', 'SN'))) %>%
    dplyr::mutate(chr=dplyr::case_when(.data$chr=='auto' ~ 'Auto',
                          .data$chr=='chrx' ~ 'ChrX',
                          .data$chr=='chry' ~ 'ChrY')) %>%
    dplyr::mutate(chr=factor(.$chr, levels=c('Auto', 'ChrX', 'ChrY'))) %>%
    dplyr::left_join(sif, by='sample.id') %>%
    dplyr::filter(!(.data$chr=='ChrY' & .data$sex=='female')) %>%
    dplyr::filter(!is.na(.data$value)) %>%
    dplyr::mutate(sex=stringr::str_to_title(.data$sex))
  
  plt <- ggplot2::ggplot(signal.noise.df.l, ggplot2::aes(x=.data$lf, y=.data$value)) +
    ggplot2::geom_violin() +
    ggplot2::geom_point(ggplot2::aes(col=.data$sex), position=ggbeeswarm::position_beeswarm()) +
    ggplot2::ylim(0, NA) +
    ggh4x::facet_grid2(.data$metric ~ .data$chr, scales='free', space='free_x', independent='y') +
    ggplot2::labs(y='Value', col='Sex') +
    ggplot2::theme_bw(base_size=30) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, vjust =1, hjust=1)) +
    ggplot2::theme(axis.title.x=ggplot2::element_blank())
  
  if (is.null(output_dir)) {
    print(plt)
  } else {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    ggplot2::ggsave(plt, file=file.path(output_dir, 'SignalNoise.jpg'), width=24, height=20)
  }

  return(invisible(signal.noise.df))
}

#' Calculate the signal and noise profile of a sample
#'
#' This function calculates the signal and noise profile of a sample. The signal is defined
#' as the standard deviation of the median signal of each arm (grouped into autosomes, X, and Y chromosomes),
#' and the noise is defined as the median of the absolute values of the vector of
#' consecutive differences between each signal.
#' 
#' @param list A numeric vector representing the signal values for a sample
#' @param pif Tibble containing probe information.
#' Should have a 'chr' column with chromosome names (1-22, X, Y) and an 'arm' column with arm names (p or q)
#' @param pif A data frame containing the positional information for each locus, including chromosome
#'
#' @keywords internal
calc.signal.noise <- function(list, pif) {
  data <- list %>%
    as.data.frame() %>%
    stats::setNames('signal') %>%
    dplyr::bind_cols(pif)

  signal.auto <- data %>%
    dplyr::filter(!.data$chr %in% c('X', 'Y')) %>%
    dplyr::group_by(.data$chr, .data$arm) %>%
    dplyr::summarize(arm.median=stats::median(.data$signal)) %>%
    dplyr::ungroup() %>%
    dplyr::summarize(signal=stats::sd(.data$arm.median)) %>%
    dplyr::pull(.data$signal)

  noise.auto <- data %>% 
    dplyr::filter(!.data$chr %in% c('X', 'Y')) %>%
    dplyr::pull(.data$signal) %>%
    diff() %>%
    abs() %>%
    stats::median()

  noise.x <- data %>% 
    dplyr::filter(.data$chr=='X') %>%
    dplyr::pull(.data$signal) %>%
    diff() %>%
    abs() %>%
    stats::median()

  noise.y <- data %>%
    dplyr::filter(.data$chr=='Y') %>%
    dplyr::pull(.data$signal) %>%
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