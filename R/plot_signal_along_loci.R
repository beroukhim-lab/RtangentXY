### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: September 2025
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Plot signal along genome before and after TangentXY normalization
#' 
#' @inheritParams run_tangent
#' @param sample_id Sample ID to plot
#' @param pif_df Tibble or filepath to a text file containing probe information.
#' Should have a 'chr' column with chromosome names (1-22, X, Y) and an 'arm' column with arm names (p or q)
#' @param tnorm A normalized tumor signal matrix output by \code{\link{run_tangent}}, or a list of such matrices
#' @param n_latent A numeric vector of the numbers of latent factors used in TangentXY normalization.
#' Should correspond to the elements in `tnorm`.
#' @param output_dir Directory to save the plot. If `NULL`, the plot will be printed to the screen.
#' 
#' @examples
#' n_latent <- c(5, 10)
#' tangent_res <- lapply(n_latent, function(nlf) {
#'   run_tangent(example_sif, example_nsig_df, example_tsig_df, nlf, make_plots = FALSE)
#' })
#' plot_signal_along_loci("tumor.female2", example_pif,
#'                        example_tsig_df, tangent_res, n_latent = n_latent)
#'
#' @returns (Invisibly) A ggplot object containing the signal profiles for the sample before and after normalization.
#'
#' @export
plot_signal_along_loci <- function(sample_id, pif_df, tsig_df, tnorm, n_latent, output_dir = NULL) {
  options(dplyr.summarise.inform = FALSE)
  if (!inherits(tnorm, "list")) {
    tnorm <- list(tnorm)
  }

  # Read in signal data and arrange by genomic position
  t.df <- read_and_format_input(tsig_df, locus_to_rownames = TRUE)
  pif <- read_and_format_input(pif_df, locus_to_rownames = FALSE)

  # Associate each tangent output with its corresponding number of latent factors
  names(tnorm) <- paste("nlf", n_latent, sep="_")

  # Sort latent factors in increasing order for plotting
  num.lf <- sort(n_latent)
  sample.of.interest <- sample_id

  t0 <- t.df %>%
      dplyr::select(dplyr::all_of(sample.of.interest)) %>%
      stats::setNames('signal') %>%
      tibble::rownames_to_column('locus') %>%
      dplyr::mutate(index=1:dplyr::n()) %>%
      dplyr::mutate(ln='Pre-normalization')

  for (i in 1:length(num.lf)) {
    num.lf.i <- num.lf[i]
    t <- tnorm[[paste("nlf", num.lf.i, sep="_")]] %>%
      tibble::column_to_rownames('locus') %>%
      dplyr::select(dplyr::all_of(sample.of.interest)) %>%
      stats::setNames('signal') %>%
      tibble::rownames_to_column('locus') %>%
      dplyr::mutate(index=1:dplyr::n()) %>%
      dplyr::mutate(ln=as.character(num.lf.i))
    if (i==1) {
      t.df <- dplyr::bind_rows(t0, t)
    } else {
      t.df <- dplyr::bind_rows(t.df, t)
    }
  }

  data.soi <- t.df %>%
    dplyr::filter(!is.na(.data$signal)) %>%
    dplyr::left_join(pif, by='locus') %>%
    dplyr::mutate(chr.class=dplyr::case_when(.data$chr %in% c(1,3,5,7,9,11,13,15,17,19,21,'X') ~ 'odd', TRUE ~ 'even')) %>%
    dplyr::group_by(.data$chr, .data$arm, .data$ln) %>%
    dplyr::mutate(median=stats::median(.data$signal)) %>%
    dplyr::ungroup()
  data.soi <- data.soi %>%
    dplyr::mutate(chr=factor(.data$chr, levels=.data$chr %>% unique() %>% gtools::mixedsort())) %>%
    dplyr::arrange(.data$chr, .data$start) %>%
    dplyr::mutate(ln=factor(.data$ln, levels=c('Pre-normalization', as.character(num.lf))))

  plt <- ggplot2::ggplot(data.soi, ggplot2::aes(x=.data$index, y=.data$signal)) +
    ggplot2::geom_hline(yintercept=0, linetype='dashed') +
    ggplot2::geom_point(ggplot2::aes(col=.data$chr.class), show.legend=FALSE) +
    ggplot2::geom_line(ggplot2::aes(y=.data$median), col='red') +
    ggplot2::scale_color_manual(values=c('odd'='black', 'even'='green')) +
    ggplot2::facet_wrap(~ln, nrow=1) +
    ggplot2::labs(y='Signal') +
    ggplot2::theme_bw(base_size=20)
  
  if (is.null(output_dir)) {
    print(plt)
  } else {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive=TRUE)
    }
    ggplot2::ggsave(plt,
      file = file.path(output_dir, paste0(sample.of.interest, "_signal.jpg")),
      width = 8 * (length(num.lf) + 1), height = 2 * (length(num.lf) + 1)
    )
  }

  return(invisible(plt))
}