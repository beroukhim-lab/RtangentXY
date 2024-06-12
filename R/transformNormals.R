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

transform_normals <- function(sif_df, nsig_df) {

  cat('\nTransforming normals ...\n')

  # Load data
  if (inherits(sif_df, "character")) {
    sif <- readr::read_delim(sif_df, progress=FALSE, show_col_types=FALSE)
  } else { sif <- sif_df }

  if (inherits(nsig_df, "character")) {
    n.df <- readr::read_delim(nsig_df, progress=FALSE, show_col_types=FALSE) %>%
      tibble::column_to_rownames('locus')
  } else {
    if ('locus' %in% colnames(nsig_df)) { n.df <- nsig_df %>% tibble::column_to_rownames('locus') }
    else { n.df <- nsig_df }
  }

  ## Linear transformation only on male chrX in normal samples
  female.normal.samples <- sif %>%
    dplyr::filter(gender=='female' & type=='normal') %>%
    dplyr::pull(sample.id)
  male.normal.samples <- sif %>%
    dplyr::filter(gender=='male' & type=='normal') %>%
    dplyr::pull(sample.id)

  female.x.mean <- n.df[grepl('^X', rownames(n.df)),] %>%
    dplyr::select(all_of(female.normal.samples)) %>%
    as.matrix() %>%
    mean()
  male.x.mean <- n.df[grepl('^X', rownames(n.df)),] %>%
    dplyr::select(all_of(male.normal.samples)) %>%
    as.matrix() %>%
    mean()
  female.x.sd <- n.df[grepl('^X', rownames(n.df)),] %>%
    dplyr::select(all_of(female.normal.samples)) %>%
    as.matrix() %>%
    sd()
  male.x.sd <- n.df[grepl('X', rownames(n.df)),] %>%
    dplyr::select(all_of(male.normal.samples)) %>%
    as.matrix() %>%
    sd()

  n.df.x.transformed <- n.df[grepl('^X|^Y', rownames(n.df)), ] %>%
    tibble::rownames_to_column('locus') %>%
    tidyr::separate(col=locus, into=c('chr', 'pos'), sep=':') %>%
    tidyr::pivot_longer(names_to='sample.id', values_to='signal', cols=-c('chr', 'pos')) %>%
    dplyr::mutate(signal=case_when(sample.id %in% male.normal.samples & chr=='X' ~ ((signal - male.x.mean)/male.x.sd) * female.x.sd + female.x.mean,
                            TRUE ~ signal)) %>%
    tidyr::pivot_wider(names_from='sample.id', values_from='signal') %>%
    tidyr::unite(col=locus, c('chr', 'pos'), sep=':') %>%
    tibble::column_to_rownames('locus')

  n.df.transformed <- n.df[!grepl('^X|^Y', rownames(n.df)), ] %>%
    dplyr::bind_rows(n.df.x.transformed) %>%
    tibble::rownames_to_column('locus')

  return(n.df.transformed)
}
