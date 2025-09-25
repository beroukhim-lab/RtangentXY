### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: March 19, 2024
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Transform Chromosome X signals
#'
#' @description
#' Conduct linear transformation on the CN signals of the male X chromosomes.
#'
#' @param sif_df The tibble of or the filepath for the sample information file
#' @param nsig_df The tibble of or the filepath for the normal signal matrix file
#'
#' @return A matrix with the male X chromosomes linearly transformed
#'
#' @import readr
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @importFrom stats sd
#' @export

linear_transformation <- function(sif_df, nsig_df, tsig_df) {

  cat('\nApplying linear transformation ...\n')
  res <- list()

  ## Load data
  if (inherits(sif_df, "character")) {
    sif <- readr::read_delim(sif_df, progress=FALSE, show_col_types=FALSE)
  } else { sif <- sif_df }

  # Normal
  if (inherits(nsig_df, "character")) {
    n.df <- readr::read_delim(nsig_df, progress=FALSE, show_col_types=FALSE) %>%
      tibble::column_to_rownames('locus')
  } else {
    if ('locus' %in% colnames(nsig_df)) { n.df <- nsig_df %>% tibble::column_to_rownames('locus') }
    else { n.df <- nsig_df }
  }

  # Tumor
  if (inherits(tsig_df, "character")) {
    t.df <- readr::read_delim(tsig_df, progress=FALSE, show_col_types=FALSE) %>%
      tibble::column_to_rownames('locus')
  } else {
    if ('locus' %in% colnames(tsig_df)) { t.df <- tsig_df %>% tibble::column_to_rownames('locus') }
    else { t.df <- tsig_df }
  }

  # Consider the case where there are only female samples
  # Does it also work when there are only male samples?

  ## Linear transformation on male chrX in normal samples
  female.normal.samples <- sif %>%
    dplyr::filter(gender=='female' & type=='normal') %>%
    dplyr::pull(sample.id)

  male.normal.samples <- sif %>%
    dplyr::filter(gender=='male' & type=='normal') %>%
    dplyr::pull(sample.id)

  female.x.mean <- n.df[grepl('^X', rownames(n.df)), , drop = FALSE] %>%
    dplyr::select(all_of(female.normal.samples)) %>%
    as.matrix() %>%
    mean()

  male.x.mean <- n.df[grepl('^X', rownames(n.df)), , drop = FALSE] %>%
    dplyr::select(all_of(male.normal.samples)) %>%
    as.matrix() %>%
    mean()

  n.df.x.transformed <- n.df[grepl('^X', rownames(n.df)), , drop = FALSE] %>%
    rownames_to_column('locus') %>%
    separate(col=locus, into=c('chr', 'pos'), sep=':') %>%
    pivot_longer(names_to='sample.id', values_to='signal', cols=-c('chr', 'pos')) %>%
    mutate(signal=case_when(sample.id %in% male.normal.samples ~ signal - male.x.mean + female.x.mean,
                            TRUE ~ signal)) %>%
    pivot_wider(names_from='sample.id', values_from='signal') %>%
    unite(col=locus, c('chr', 'pos'), sep=':') %>%
    column_to_rownames('locus')

  ## Linear transformation on male chrY in normal samples
  male.y.mean <- n.df[grepl('^Y', rownames(n.df)), , drop = FALSE] %>%
    select(all_of(male.normal.samples)) %>%
    as.matrix() %>%
    mean()

  male.y.mean.mode <- n.df[grepl('^Y', rownames(n.df)), , drop = FALSE] %>%
    select(all_of(male.normal.samples)) %>%
    as.matrix() %>%
    apply(., 2, mean) %>%
    density() %>%
    {.$x[which.max(.$y)]}

  n.df.y.transformed <- n.df[grepl('^Y', rownames(n.df)), , drop = FALSE] %>%
    rownames_to_column('locus') %>%
    separate(col=locus, into=c('chr', 'pos'), sep=':') %>%
    pivot_longer(names_to='sample.id', values_to='signal', cols=-c('chr', 'pos')) %>%
    group_by(sample.id) %>%
    mutate(sample_mean=mean(signal), sample_sd=sd(signal)) %>%
    ungroup() %>%
    mutate(signal=case_when(sample.id %in% male.normal.samples ~ signal - sample_mean + male.y.mean.mode,
                                        TRUE ~ signal)) %>%
    select(chr, pos, sample.id, signal) %>%
    pivot_wider(names_from='sample.id', values_from='signal') %>%
    unite(col=locus, c('chr', 'pos'), sep=':') %>%
    column_to_rownames('locus')

  n.df.transformed <- n.df[!grepl('^X|^Y', rownames(n.df)), , drop = FALSE] %>%
    bind_rows(n.df.x.transformed) %>%
    bind_rows(n.df.y.transformed) %>%
    rownames_to_column('locus')

  res$n.df <- n.df.transformed

  ## Linear transformation on male chrX in tumor samples
  male.tumor.samples <- sif %>%
    filter(gender=='male' & type=='tumor') %>%
    pull(sample.id)

  t.df.x.transformed <- t.df[grepl('^X|^Y', rownames(t.df)), , drop = FALSE] %>%
    rownames_to_column('locus') %>%
    separate(col=locus, into=c('chr', 'pos'), sep=':') %>%
    pivot_longer(names_to='sample.id', values_to='signal', cols=-c('chr', 'pos')) %>%
    mutate(signal=case_when(sample.id %in% male.tumor.samples & chr=='X' ~ signal - male.x.mean + female.x.mean,
                            TRUE ~ signal)) %>%
    pivot_wider(names_from='sample.id', values_from='signal') %>%
    unite(col=locus, c('chr', 'pos'), sep=':') %>%
    column_to_rownames('locus')

  t.df.transformed <- t.df[!grepl('^X|^Y', rownames(t.df)), , drop = FALSE] %>%
    bind_rows(t.df.x.transformed) %>%
    rownames_to_column('locus')

  res$t.df <- t.df.transformed

  return(res)
}
