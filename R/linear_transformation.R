### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: September 2025
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Apply linear transformation to male chromosome X and Y signals
#'
#' Applies linear transformation to male normal sample chromosome X and Y signals.
#' Also applies linear transformation to male tumor sample chromosome X signals.
#' 
#' @param sif Tibble or data.frame containing sample metadata.
#' @param n.df Tibble or data.frame containing the normal signal matrix. Rownames
#' should be locus names in the format `"{chr}:{start}-{end}"` where `{chr}` is 1-22, X, or Y
#' and where `{start}` and `{end}` are genomic coordinates
#' @param t.df Tibble or data.frame containing the tumor signal matrix. Rownames
#' should match `n.df`.
#'
#' @return A list with two data frames, both having a `locus` column:
#' \describe{
#'  \item{n.df: }{Normal signal matrix with male chromosome X and Y signals linearly transformed}
#'  \item{t.df: }{Tumor signal matrix with male chromosome X signals linearly transformed}
#' }
# @examples
# res <- linear_transformation(example_sif,
#  example_nsig_df %>% tibble::column_to_rownames('locus'),
#  example_tsig_df %>% tibble::column_to_rownames('locus'))
# head(res$n.df)
# head(res$t.df)
#'
#' @keywords internal
linear_transformation <- function(sif, n.df, t.df) {

  cat('\nApplying linear transformation ...\n')
  res <- list()

  ## Linear transformation on male chrX in normal samples
  female.normal.samples <- sif %>%
    dplyr::filter(.data$sex=='female' & .data$type=='normal') %>%
    dplyr::pull(.data$sample.id)

  male.normal.samples <- sif %>%
    dplyr::filter(.data$sex=='male' & .data$type=='normal') %>%
    dplyr::pull(.data$sample.id)

  female.x.mean <- n.df[grepl('^X', rownames(n.df)), , drop = FALSE] %>%
    dplyr::select(dplyr::all_of(female.normal.samples)) %>%
    as.matrix() %>%
    mean()

  male.x.mean <- n.df[grepl('^X', rownames(n.df)), , drop = FALSE] %>%
    dplyr::select(dplyr::all_of(male.normal.samples)) %>%
    as.matrix() %>%
    mean()

  n.df.x.transformed <- n.df[grepl('^X', rownames(n.df)), , drop = FALSE] %>%
    tibble::rownames_to_column('locus') %>%
    tidyr::separate(col="locus", into=c('chr', 'pos'), sep=':') %>%
    tidyr::pivot_longer(names_to='sample.id', values_to='signal', cols=-c('chr', 'pos')) %>%
    dplyr::mutate(signal=dplyr::case_when(.data$sample.id %in% male.normal.samples ~ .data$signal - male.x.mean + female.x.mean,
                            TRUE ~ .data$signal)) %>%
    tidyr::pivot_wider(names_from='sample.id', values_from='signal') %>%
    tidyr::unite(col="locus", c('chr', 'pos'), sep=':') %>%
    tibble::column_to_rownames('locus')

  ## Linear transformation on male chrY in normal samples
  male.y.mean <- n.df[grepl('^Y', rownames(n.df)), , drop = FALSE] %>%
    dplyr::select(dplyr::all_of(male.normal.samples)) %>%
    as.matrix() %>%
    mean()

  male.y.mean.mode <- n.df[grepl('^Y', rownames(n.df)), , drop = FALSE] %>%
    dplyr::select(dplyr::all_of(male.normal.samples)) %>%
    as.matrix() %>%
    apply(., 2, mean) %>%
    stats::density() %>%
    {.$x[which.max(.$y)]}

  n.df.y.transformed <- n.df[grepl('^Y', rownames(n.df)), , drop = FALSE] %>%
    tibble::rownames_to_column('locus') %>%
    tidyr::separate(col="locus", into=c('chr', 'pos'), sep=':') %>%
    tidyr::pivot_longer(names_to='sample.id', values_to='signal', cols=-c('chr', 'pos')) %>%
    dplyr::group_by(.data$sample.id) %>%
    dplyr::mutate(sample_mean=mean(.data$signal)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(signal=dplyr::case_when(.data$sample.id %in% male.normal.samples ~ .data$signal - .data$sample_mean + male.y.mean.mode,
                                        TRUE ~ .data$signal)) %>%
    dplyr::select("chr", "pos", "sample.id", "signal") %>%
    tidyr::pivot_wider(names_from='sample.id', values_from='signal') %>%
    tidyr::unite(col="locus", c('chr', 'pos'), sep=':') %>%
    tibble::column_to_rownames('locus')

  n.df.transformed <- n.df[!grepl('^X|^Y', rownames(n.df)), , drop = FALSE] %>%
    dplyr::bind_rows(n.df.x.transformed) %>%
    dplyr::bind_rows(n.df.y.transformed) %>%
    tibble::rownames_to_column('locus')

  res$n.df <- n.df.transformed

  ## Linear transformation on male chrX in tumor samples
  male.tumor.samples <- sif %>%
    dplyr::filter(.data$sex=='male' & .data$type=='tumor') %>%
    dplyr::pull(.data$sample.id)

  t.df.x.transformed <- t.df[grepl('^X|^Y', rownames(t.df)), , drop = FALSE] %>%
    tibble::rownames_to_column('locus') %>%
    tidyr::separate(col="locus", into=c('chr', 'pos'), sep=':') %>%
    tidyr::pivot_longer(names_to='sample.id', values_to='signal', cols=-c('chr', 'pos')) %>%
    dplyr::mutate(signal=dplyr::case_when(.data$sample.id %in% male.tumor.samples & .data$chr=='X' ~ .data$signal - male.x.mean + female.x.mean,
                            TRUE ~ .data$signal)) %>%
    tidyr::pivot_wider(names_from='sample.id', values_from='signal') %>%
    tidyr::unite(col="locus", c('chr', 'pos'), sep=':') %>%
    tibble::column_to_rownames('locus')

  t.df.transformed <- t.df[!grepl('^X|^Y', rownames(t.df)), , drop = FALSE] %>%
    dplyr::bind_rows(t.df.x.transformed) %>%
    tibble::rownames_to_column('locus')

  res$t.df <- t.df.transformed

  return(res)
}
