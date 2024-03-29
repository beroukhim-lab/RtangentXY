### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: March 29, 2024
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Plots the importance of latent factors  `run_svd()`
#'
#' @param n_svd The returned matrix decomposition of `run_svd()`
#' @param top_100 A boolean. Set to TRUE if interested in the top 100 latent factors. If set to TRUE and if there are fewer than 100 latent factors, an error will be thrown
#'
#' @return A boolean of whether the plot was successfully generated.
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom scales pretty_breaks
#' @importFrom stats setNames
#'
#' @export

plot_latent_factors <- function(n_svd, top_100 = FALSE) {
  n.autox.svd <- n_svd

  d.df <- n.autox.svd$d %>%
    as.data.frame() %>%
    setNames('d') %>%
    dplyr::mutate(n=1:n())

  if (top_100) {
    if (nrow(d.df) >= 100) {
      g <- ggplot2::ggplot(d.df, aes(x=n, y=d)) +
        ggplot2::geom_line() +
        ggplot2::scale_x_continuous(breaks=scales::pretty_breaks()) +
        ggplot2::labs(title='Importance of latent factors', y='r (Importance)', col='# of normals') +
        ggplot2::theme_bw(base_size=30) +
        ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              axis.title.x=element_blank())
      print(g)
      return (TRUE)
    } else {
      cat("ERROR: Fewer than 100 latent factors\n")
      return (FALSE)
    }

  } else{
    g <- ggplot2::ggplot(d.df, aes(x=n, y=d)) +
      ggplot2::geom_line() +
      ggplot2::scale_x_continuous(breaks=scales::pretty_breaks()) +
      ggplot2::labs(title='Importance of latent factors', y='r (Importance)', col='# of normals') +
      ggplot2::theme_bw(base_size=30) +
      ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            axis.title.x=element_blank())
    print(g)
    return (TRUE)
  }
}
