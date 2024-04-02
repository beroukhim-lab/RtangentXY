### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: April 2, 2024
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' For each of the given number of latent factors, this function runs `run_tangent()` and plots the signal and noise data for each latent factor.
#'
#' @param sif_filepath The filepath for the sample information file
#' @param pif_filepath The filepath for the probe information file
#' @param ndf_filepath The filepath for the normal signal matrix file
#' @param tdf_filepath The filepath for the tumor signal matrix file
#' @param n_latents A vector of latent factors to reconstruct the normal subspace
#' @param cores The number of parallel cores to use when calculating signal to noise. Default is set to 2 cores
#'
#' @returns A matrix of the signal and noise data that is plotted with this function as well as the `run_tangent()` outputs. These are stored in a list under `signal_noise_df` and `tangent_out_df` respectively.
#'
#' @import readr
#' @import tibble
#' @import dplyr
#' @import ggplot2
#' @import parallel
#' @import stringr
#' @importFrom tidyselect matches
#' @importFrom stats setNames median sd
#' @importFrom ggbeeswarm position_beeswarm
#' @importFrom ggh4x facet_grid2
#'
#' @export

plot_signal_noise <- function(sif_filepath, pif_filepath, ndf_filepath, tdf_filepath, n_latents, cores = 2) {
  # Load data
  sif <- readr::read_delim(sif_filepath, progress=FALSE, show_col_types=FALSE)
  pif <- readr::read_delim(pif_filepath, progress=FALSE, show_col_types=FALSE)
  t.df <- readr::read_delim(tdf_filepath, progress=FALSE, show_col_types=FALSE) %>%
    tibble::column_to_rownames('locus')

  num.lf <- sort(n_latents)

  male.tumors <- sif %>%
    dplyr::filter(sample.id %in% colnames(t.df) & gender=='male') %>%
    dplyr::pull(sample.id)

  cat('Calculating signal and noise in pre-normalization data...\n')
  t0.signal.noise.list <- parallel::mclapply(t.df, calc.signal.noise, pif, mc.cores = cores)
  t0.signal.noise.df <- t0.signal.noise.list %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(sample.id=names(t0.signal.noise.list)) %>%
    dplyr::mutate(lf='Pre-normalization')

  ## Post-TangentXY
  for (i in 1:length(num.lf)) {
    num.lf.i <- num.lf[i]
    if (num.lf.i==1) {
      cat(paste0('Calculating signal and noise in normalized data with ', num.lf.i, ' latent factor...\n'))
    } else {
      cat(paste0('Calculating signal and noise in normalized data with ', num.lf.i, ' latent factors...\n'))
    }

    # Run Tangent XY
    t.df.i <- run_tangent(sif_filepath, ndf_filepath, tdf_filepath, num.lf.i) %>%
      tibble::column_to_rownames('locus')

    t.i.signal.noise.list <- parallel::mclapply(t.df.i, calc.signal.noise, pif, mc.cores = cores)
    t.i.signal.noise.df <- t.i.signal.noise.list %>%
      bind_rows() %>%
      mutate(sample.id=names(t.i.signal.noise.list)) %>%
      mutate(lf=as.character(num.lf.i))

    if (i==1) {
      signal.noise.df <- t0.signal.noise.df %>% bind_rows(t.i.signal.noise.df)
      i.names <- c(toString(num.lf.i))
      t.df.combined <- list(t.df.i)
    } else {
      signal.noise.df <- signal.noise.df %>% bind_rows(t.i.signal.noise.df)
      i.names <- append(i.names, toString(num.lf.i))
      t.df.combined[[length(t.df.combined) + 1]] <- t.df.i
    }
  }
  t.df.combined <- setNames(t.df.combined, i.names)

  # Plotting signal noise
  signal.noise.df.l <- signal.noise.df %>%
    tidyr::pivot_longer(cols=tidyselect::matches('^signal\\.|^noise\\.|^sn\\.'), names_pattern='(^signal\\.|^noise\\.|^sn\\.)(.*)$', names_to=c('metric', 'chr')) %>%
    dplyr::mutate(lf=dplyr::case_when(chr=='chry' & lf!='Pre-normalization' ~ 'TangentXY', TRUE ~ lf)) %>%
    dplyr::distinct(sample.id, lf, metric, chr, .keep_all=TRUE) %>%
    dplyr::mutate(lf=factor(.$lf, levels=c('Pre-normalization', num.lf, 'TangentXY'))) %>%
    dplyr::mutate(metric=dplyr::case_when(metric=='signal.' ~ 'Signal',
                            metric=='noise.' ~ 'Noise',
                            metric=='sn.' ~ 'SN')) %>%
    dplyr::mutate(metric=factor(.$metric, levels=c('Signal', 'Noise', 'SN'))) %>%
    dplyr::mutate(chr=dplyr::case_when(chr=='auto' ~ 'Auto',
                         chr=='chrx' ~ 'ChrX',
                         chr=='chry' ~ 'ChrY')) %>%
    dplyr::mutate(chr=factor(.$chr, levels=c('Auto', 'ChrX', 'ChrY'))) %>%
    dplyr::left_join(sif, by='sample.id') %>%
    dplyr::filter(!(chr=='ChrY' & gender=='female')) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::mutate(gender=stringr::str_to_title(gender))

  g <- ggplot2::ggplot(signal.noise.df.l, ggplot2::aes(x=lf, y=value)) +
    ggplot2::geom_violin() +
    ggplot2::geom_point(ggplot2::aes(col=gender), position=ggbeeswarm::position_beeswarm()) +
    ggplot2::ylim(0, NA) +
    ggh4x::facet_grid2(metric ~ chr, scales='free', space='free_x', independent='y') +
    ggplot2::labs(y='Value', col='Gender') +
    ggplot2::theme_bw(base_size=30) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, vjust =1, hjust=1)) +
    ggplot2::theme(axis.title.x=ggplot2::element_blank())

  print(g)

  # Format return list
  ret.list <- list(signal_noise_df = signal.noise.df,
                   tangent_out_df = t.df.combined)
  return(ret.list)
}

## Make function for calculating signal and noise
## The following function is used in plot_signal_noise()
calc.signal.noise <- function(list, pif) {
  data <- list %>%
    as.data.frame() %>%
    setNames('signal') %>%
    dplyr::bind_cols(pif)

  signal.auto <- data %>%
    dplyr::filter(!chr %in% c('X', 'Y')) %>%
    dplyr::group_by(chr, arm) %>%
    dplyr::summarize(arm.median=median(signal)) %>%
    dplyr::ungroup() %>%
    dplyr::summarize(signal=sd(arm.median)) %>%
    dplyr::pull(signal)

  noise.auto <- data %>%
    dplyr::filter(!chr %in% c('X', 'Y')) %>%
    dplyr::pull(signal) %>%
    diff() %>%
    abs() %>%
    median()

  noise.x <- data %>%
    dplyr::filter(chr=='X') %>%
    dplyr::pull(signal) %>%
    diff() %>%
    abs() %>%
    median()

  noise.y <- data %>%
    dplyr::filter(chr=='Y') %>%
    dplyr::pull(signal) %>%
    diff() %>%
    abs() %>%
    median()

  result.df <- data.frame(signal.auto=signal.auto,
                          noise.auto=noise.auto,
                          sn.auto=signal.auto/noise.auto,
                          noise.chrx=noise.x,
                          noise.chry=noise.y)

  return(result.df)
}
