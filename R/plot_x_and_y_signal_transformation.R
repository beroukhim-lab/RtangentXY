### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: September 2025
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Plot Chromosome X and Y signals before and after linear transformation
#'
#' @param sif Tibble or data.frame containing sample metadata.
#' @param n.df Tibble or data.frame containing the normal signal matrix. Rownames
#' should be locus names in the format `"{chr}:{start}-{end}"` where `{chr}` is 1-22, X, or Y
#' and where `{start}` and `{end}` are genomic coordinates
#' @param t.df Tibble or data.frame containing the tumor signal matrix. Rownames
#' should match `n.df`.
#' @param n.df.transformed Tibble or data.frame containing the normal signal matrix after linear transformation.
#' Output of \link{linear_transformation}.
#' @param t.df.transformed Tibble or data.frame containing the tumor signal matrix after linear transformation.
#' Output of \link{linear_transformation}.
#' @param output_dir Directory to save the plots. If `NULL`, the plots will be printed to the screen.
#' 
#' @returns Nothing. Generates and saves or prints plots.
#'
#' @keywords internal
plot_x_and_y_signal_transformation <- function(
  sif, n.df, t.df,
  n.df.transformed, t.df.transformed,
  output_dir = NULL
) {
  # Chromosome X
  signaln.n.x.before <- n.df[grepl('^X', rownames(n.df)), , drop = FALSE] %>%
    tidyr::pivot_longer(names_to='sample.id', values_to='signal', cols=dplyr::everything()) %>%
    dplyr::left_join(sif, by='sample.id') %>%
    dplyr::mutate(transformation='Before transformation') %>%
    dplyr::mutate(type='Normal')

  signal.n.x.after <- n.df.transformed %>%
    dplyr::filter(grepl('^X', .data$locus)) %>%
    tibble::column_to_rownames('locus') %>%
    tidyr::pivot_longer(names_to='sample.id', values_to='signal', cols=dplyr::everything()) %>%
    dplyr::left_join(sif, by='sample.id') %>%
    dplyr::mutate(transformation='After transformation') %>%
    dplyr::mutate(type='Normal')

  signal.t.x.before <- t.df[grepl('^X', rownames(t.df)), , drop = FALSE] %>%
    tidyr::pivot_longer(names_to='sample.id', values_to='signal', cols=dplyr::everything()) %>%
    dplyr::left_join(sif, by='sample.id') %>%
    dplyr::mutate(transformation='Before transformation') %>%
    dplyr::mutate(type='Tumor')

  signal.t.x.after <- t.df.transformed %>%
    dplyr::filter(grepl('^X', .data$locus)) %>%
    tibble::column_to_rownames('locus') %>%
    tidyr::pivot_longer(names_to='sample.id', values_to='signal', cols=dplyr::everything()) %>%
    dplyr::left_join(sif, by='sample.id') %>%
    dplyr::mutate(transformation='After transformation') %>%
    dplyr::mutate(type='Tumor')

  signal.x <- signaln.n.x.before %>%
    dplyr::bind_rows(signal.n.x.after) %>%
    dplyr::bind_rows(signal.t.x.before) %>%
    dplyr::bind_rows(signal.t.x.after)
  signal.x <- signal.x %>%
    dplyr::mutate(transformation=factor(.data$transformation, levels=c('Before transformation', 'After transformation'))) %>%
    dplyr::mutate(gender=stringr::str_to_title(.data$sex))

  plt <- ggplot2::ggplot(signal.x, ggplot2::aes(x=.data$signal, group=.data$sample.id)) +
    ggplot2::geom_density(ggplot2::aes(fill=.data$sex), alpha=0.25) +
    ggplot2::geom_vline(xintercept=0, col='red', linetype='dashed') +
    ggplot2::geom_vline(xintercept=-1, col='blue', linetype='dashed') +
    ggplot2::coord_flip() +
    ggplot2::facet_grid(.data$type~.data$transformation) +
    ggplot2::labs(title='Chromosome X signal distribution', fill='Sex') +
    ggplot2::theme_bw(base_size=20)

  if (is.null(output_dir)) {
    print(plt)
  } else {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    ggplot2::ggsave(plt, file=file.path(output_dir, 'ChrX_SignalDistribution.jpg'), width=12, height=8)
  }

  # Chromosome Y
  signaln.n.y.before <- n.df[grepl('^Y', rownames(n.df)), , drop = FALSE] %>%
    tidyr::pivot_longer(names_to='sample.id', values_to='signal', cols=dplyr::everything()) %>%
    dplyr::left_join(sif, by='sample.id') %>%
    dplyr::mutate(transformation='Before transformation') %>%
    dplyr::mutate(type='Normal')

  signal.n.y.after <- n.df.transformed %>%
    dplyr::filter(grepl('^Y', .data$locus)) %>%
    tibble::column_to_rownames('locus') %>%
    tidyr::pivot_longer(names_to='sample.id', values_to='signal', cols=dplyr::everything()) %>%
    dplyr::left_join(sif, by='sample.id') %>%
    dplyr::mutate(transformation='After transformation') %>%
    dplyr::mutate(type='Normal')

  signal.t.y.before <- t.df[grepl('^Y', rownames(t.df)), , drop = FALSE] %>%
    tidyr::pivot_longer(names_to='sample.id', values_to='signal', cols=dplyr::everything()) %>%
    dplyr::left_join(sif, by='sample.id') %>%
    dplyr::mutate(transformation='Before transformation') %>%
    dplyr::mutate(type='Tumor')

  signal.t.y.after <- t.df.transformed %>%
    dplyr::filter(grepl('^Y', .data$locus)) %>%
    tibble::column_to_rownames('locus') %>%
    tidyr::pivot_longer(names_to='sample.id', values_to='signal', cols=dplyr::everything()) %>%
    dplyr::left_join(sif, by='sample.id') %>%
    dplyr::mutate(transformation='After transformation') %>%
    dplyr::mutate(type='Tumor')

  signal.y <- signaln.n.y.before %>%
    dplyr::bind_rows(signal.n.y.after) %>%
    dplyr::bind_rows(signal.t.y.before) %>%
    dplyr::bind_rows(signal.t.y.after)
  signal.y <- signal.y %>%
    dplyr::mutate(transformation=factor(.data$transformation, levels=c('Before transformation', 'After transformation'))) %>%
    dplyr::mutate(gender=stringr::str_to_title(.data$sex))

  plt <- ggplot2::ggplot(signal.y %>% dplyr::filter(.data$sex=='male'), ggplot2::aes(x=.data$signal, group=.data$sample.id)) +
    ggplot2::geom_density(ggplot2::aes(fill=.data$sex), alpha=0.25) +
    ggplot2::geom_vline(xintercept=0, col='red', linetype='dashed') +
    ggplot2::geom_vline(xintercept=-1, col='blue', linetype='dashed') +
    ggplot2::coord_flip() +
    ggplot2::facet_grid(.data$type~.data$transformation) +
    ggplot2::scale_fill_manual(values=c('male'='#00BFC4')) +
    ggplot2::labs(title='Chromosome Y signal distribution', fill='Sex') +
    ggplot2::theme_bw(base_size=20)
  
  if (is.null(output_dir)) {
    print(plt)
  } else {
    ggplot2::ggsave(plt, file=file.path(output_dir, 'ChrY_SignalDistribution.jpg'), width=12, height=8)
  }
}
