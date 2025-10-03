### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: September 2025
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Plot importance of latent factors
#' 
#' @inheritParams run_tangent
#' @param output_dir Directory to save the plots. If `NULL`, the plots will be printed to the screen.
#'
#' @returns (Invisibly) A list of ggplot objects. Component `all` contains the plot of all latent factors.
#' Component `top_100` contains the plot of the top 100 latent factors (if there are at least 100 latent factors).
#' 
#' @export
plot_latent_factor_importance <- function(sif_df, nsig_df, tsig_df, output_dir = NULL) {
  # Load data
  if (inherits(sif_df, "character")) {
    sif <- readr::read_delim(sif_df, progress=FALSE, show_col_types=FALSE)
  } else { sif <- sif_df }
  # Read in signal data and arrange by genomic position
  n.df <- read_and_format_input(nsig_df, locus_to_rownames = TRUE)
  t.df <- read_and_format_input(tsig_df, locus_to_rownames = TRUE)


  # Apply linear transformation
  linear.transformation <- linear_transformation(sif, n.df, t.df)
  n.df <- linear.transformation$n.df
  
  # Run SVD
  n.autox.svd <- run_svd(n.df)

  # Plot latent factors
  plot_latent_factors(n.autox.svd, output_dir)
}

#' Plot importance of latent factors
#'
#' Plot singular values to visualize the importance of latent factors
#'
#' @param n.autox.svd A list containing the SVD components (output of \link{run_svd})
#'
#' @returns A list containing ggplot objects. If `output_dir` is specified, the plots will also be saved as .jpg files in the specified directory.
#'
#' @keywords internal
plot_latent_factors <- function(n.autox.svd, output_dir = NULL) {

  d.df <- n.autox.svd$d %>%
    as.data.frame() %>%
    stats::setNames('d') %>%
    dplyr::mutate(n = 1:dplyr::n())
  
  res <- list()

  # Plot all latent factors
  res$all <- ggplot2::ggplot(d.df, ggplot2::aes(x = .data$n, y = .data$d)) +
    ggplot2::geom_line() +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
    ggplot2::labs(title = 'Importance of latent factors', y = 'r (Importance)') +
    ggplot2::theme_bw(base_size = 30) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title.x = ggplot2::element_blank())

  if (is.null(output_dir)) {
    print(res$all)
  } else {
    ggplot2::ggsave(res$all, file = file.path(output_dir, 'LatentFactors_Importance.jpg'), width = 16, height = 8)
  }

  # Plot top 100 latent factors
  if (nrow(d.df) >= 100) {
    res$top_100 <- ggplot2::ggplot(d.df %>% utils::head(100), ggplot2::aes(x = .data$n, y = .data$d)) +
        ggplot2::geom_line() +
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
        ggplot2::labs(title = 'Importance of top 100 latent factors', y = 'r (Importance)') +
        ggplot2::theme_bw(base_size = 30) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title.x = ggplot2::element_blank()
        )
    if (is.null(output_dir)) {
      print(res$top_100)
    } else {
      ggplot2::ggsave(res$top_100, file = file.path(output_dir, 'LatentFactors_Importance_Top100.jpg'), width = 16, height = 8)
    }
  }

  return(invisible(res))
}
