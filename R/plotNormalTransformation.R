### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: March 29, 2024
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Plot Chromosome X transformation
#'
#' @description
#' Plots the signal distribution of Chromosome X after the linear transformation
#' done by `transform_normals()`
#'
#' @param sif_filepath The filepath for the sample information file used to create `ndf_transformed`
#' @param ndf_filepath The filepath for the normal signal matrix file used to create `ndf_transformed`
#' @param ndf_transformed The returned transformed dataframe from `transform_normals()`
#'
#' @returns A dataframe of the plotted data
#'
#' @import readr
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import ggplot2
#' @import stringr
#'
#' @export

plot_normal_transformation <- function(sif_filepath, ndf_filepath, ndf_transformed) {
  # Load data
  sif <- readr::read_delim(sif_filepath, progress=FALSE, show_col_types=FALSE)
  n.df <- readr::read_delim(ndf_filepath, progress=FALSE, show_col_types=FALSE) %>%
    tibble::column_to_rownames('locus')
  n.df.transformed <- ndf_transformed

  signal.x.before <- n.df[grepl('^X', rownames(n.df)),] %>%
    tidyr::pivot_longer(names_to='sample.id', values_to='signal', cols=everything()) %>%
    dplyr::left_join(sif, by='sample.id') %>%
    dplyr::mutate(transformation='Before transformation')

  signal.x.after <- n.df.transformed %>%
    dplyr::filter(grepl('^X', locus)) %>%
    tibble::column_to_rownames('locus') %>%
    tidyr::pivot_longer(names_to='sample.id', values_to='signal', cols=everything()) %>%
    dplyr::left_join(sif, by='sample.id') %>%
    dplyr::mutate(transformation='After transformation')

  signal.x <- dplyr::bind_rows(signal.x.before, signal.x.after) %>%
    dplyr::mutate(transformation=factor(.$transformation, levels=c('Before transformation', 'After transformation'))) %>%
    dplyr::mutate(gender=stringr::str_to_title(gender))

  g <- ggplot2::ggplot(signal.x, ggplot2::aes(x=signal, group=sample.id)) +
    ggplot2::geom_density(ggplot2::aes(fill=gender), alpha=0.25) +
    ggplot2::geom_vline(xintercept=0, col='red', linetype='dashed') +
    ggplot2::geom_vline(xintercept=-1, col='blue', linetype='dashed') +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~transformation, nrow=1) +
    ggplot2::labs(title='Normal samples chrX signal distribution', fill='Gender') +
    ggplot2::theme_bw(base_size=20)

  print(g)

  return(as.data.frame(signal.x))
}
