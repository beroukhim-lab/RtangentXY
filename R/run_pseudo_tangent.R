### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: September 2025
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' Run the PseudoTangentXY algorithm
#'
#' @description
#' Pseudotangent also isolates the tumor signals using the Tangent algorithm,
#' but it should be used when the set of normals are particularly non-representative.
#'
#' @inheritParams run_tangent
#' @param cbs_a The alpha parameter for the function \code{\link{run_cbs}}
#' @param cbs_mw The minimum width parameter for the function \code{\link{run_cbs}}
#' @param cbs_ncores The number of cores to use for the function \code{\link{run_cbs}}
#' @param num_partitions The number of partitions to create in the pseudotangent pipeline. This
#' must be less than the number of tumors
#' @param partition_seed The seed to set for reproducibility for the random partitioning
#' @param n_latent_init The number of latent factors to reconstruct the initial normal subspace
#' @param n_latent_part The number of latent factors for each of the partition runs. This should not exceed the minimum number of normal/pseudonormal samples across all partitions.
#'
#' @returns A normalized tumor signal matrix
#'
#' @export
run_pseudotangent <- function(sif_df, nsig_df, tsig_df, n_latent_init, num_partitions, n_latent_part,
                              cbs_a = 0.005, cbs_mw = 3, cbs_ncores = 1,
                              partition_seed = 37, make_plots = FALSE, output_dir = NULL) {

  # Check to make sure the input number of partitions is not 1
  if (num_partitions == 1) {
    stop("ERROR: Must choose more than 1 partition")
  }

  # Load data
  if (inherits(sif_df, "character")) {
    sif <- readr::read_delim(sif_df, progress=FALSE, show_col_types=FALSE)
  } else { sif <- sif_df }
  # Read in signal data and arrange by genomic position
  n.df <- read_and_format_input(nsig_df, locus_to_rownames = FALSE)
  t.df <- read_and_format_input(tsig_df, locus_to_rownames = FALSE)
  

  # Check to make sure the input number of partitions is not more than the number of tumors
  if (num_partitions > length(colnames(t.df))) {
    stop("ERROR: number of partitions exceeds the number of tumors")
  }

  # Step 1: Run TangentXY on a small set of normals
  cat('\nRunning Tangent on initial set of normals...')
  if (is.null(output_dir)) {
    initial_output_dir <- NULL
  } else {
    initial_output_dir <- file.path(output_dir, "initial_pseudotangent_run")
    if (!dir.exists(initial_output_dir)) {
      dir.create(initial_output_dir, recursive = TRUE)
    }
  }
  step1_tangent_out <- run_tangent(sif, n.df, t.df, n_latent_init, make_plots = make_plots, output_dir = initial_output_dir)

  # Step 2: Run CBS on the TangentXY outputs to get the tentative CN profile
  cat('\nRunning CBS on inital Tangent output...\n')
  step2_run_cbs <- run_cbs(step1_tangent_out, alpha = cbs_a, min_width = cbs_mw, n_cores = cbs_ncores)

  # Step 3: Subtract this tentatitve CN profile
  cat('\nSubtracting tentative CN profile...\n')
  # NOTE: Check to make sure the NAs when subtracted can actually be input into Tangent
  step2_run_cbs_df <- convert_cbs(step2_run_cbs, sif, t.df)
  t.df <- t.df %>% tibble::column_to_rownames('locus')
  step3_pseudo_normal <- t.df - step2_run_cbs_df

  # Step 4: Partition samples and run Tangent
  # Only partition the tumors
  tumor_sif <- sif[sif$type == "tumor", , drop = FALSE]
  sif_partitions <- partition.dataframe(tumor_sif, num_partitions, partition_seed)

  # For each partition, create complement pseudonormals and rename with pseudonormal
  # This is in the same list format as the sif_partitions above
  pseudonormal_partitions <- lapply(sif_partitions, create.complement.pseudonormal, step3_pseudo_normal, sif)

  # For each partition, create the matrix of tumor signals
  tumor_partitions <- lapply(sif_partitions, create.df.tumors, t.df)

  # Loop through and create files individually
  # Also use this loop to modify sif_partitions to have pseudonormal profiles
  cat('\nRunning Tangent on partitions...\n')
  step4_tangent_partitions <- list()
  for (i in 1:length(sif_partitions)) {
    cat(paste0('\nParition ', toString(i), '\n'))
    # Set up all of the files necessary for Tangent
    cur_mod_sif <- add.pseudonormal.sif(sif_partitions[[i]], pseudonormal_partitions[[i]], sif)
    cur_pseudonormals <- pseudonormal_partitions[[i]] %>%
      tibble::rownames_to_column('locus')
    cur_tumors <- tumor_partitions[[i]] %>%
      tibble::rownames_to_column('locus')

    if (is.null(output_dir)) {
      cur_output_dir <- NULL
    } else {
      cur_output_dir <- file.path(output_dir, paste0("partition_", i))
      if (!dir.exists(cur_output_dir)) {
        dir.create(cur_output_dir, recursive = TRUE)
      }
    }
    cur_tangent_out <- run_tangent(cur_mod_sif, cur_pseudonormals, cur_tumors, n_latent_part, make_plots = make_plots, output_dir = cur_output_dir)
    step4_tangent_partitions[[i]] <- cur_tangent_out
  }

  # Step 5: Combine all of the partitions
  cat('\nCombining all partitions...\n')
  step5_combined_output <- partitions.format.join(step4_tangent_partitions, t.df)
  step5_combined_output <- data.frame(locus = row.names(t.df), step5_combined_output)

  return(step5_combined_output)

  # # Step 6: Run CBS on combined output
  # cat('\nRunning CBS on final output...\n')
  # step6_final_cbs_output <- run_cbs(step5_combined_output, data_type = cbs_dt, alpha = cbs_a, min_width = cbs_mw)
  # cat('Done.\n')
  #
  # return(step6_final_cbs_output)
}

partition.dataframe <- function(df, num_partitions, seed) {
  # Set the seed for reproducibility
  set.seed(seed)

  run_count <- 0
  male_count_PASS <- FALSE

  # Shuffle the rows of the dataframe
  shuffled_indices <- sample(nrow(df))
  shuffled_df <- df[shuffled_indices, , drop = FALSE]

  # Calculate the size of each partition
  partition_sizes <- rep(floor(nrow(df) / num_partitions), num_partitions)
  remainder <- nrow(df) %% num_partitions

  # Distribute the remainder across the first few partitions
  if (remainder > 0) {
    partition_sizes[1:remainder] <- partition_sizes[1:remainder] + 1
  }

  # Split the dataframe into partitions
  partitions <- vector("list", num_partitions)
  start <- 1
  for (i in 1:num_partitions) {
    end <- start + partition_sizes[i] - 1
    partitions[[i]] <- shuffled_df[start:end, , drop = FALSE]
    start <- end + 1
  }
  return(partitions)
}

create.complement.pseudonormal <- function(cur_tumor_set, pseudonormal_df, sif) {
  # Remove the tumors listed in listy[[i]] from the pseudonormal_df
  # We assume that the pseudonormal profiles in pseudonormal_df are named with the tumor names
  subset_df <- pseudonormal_df
  for (id in cur_tumor_set$sample.id) {
    subset_df <- subset_df %>% dplyr::select(-dplyr::all_of(id))
  }

  # Now we need to rename these tumor names with the respective normal name
  cur_names <- colnames(subset_df)
  new_names <- paste0(cur_names, "_pseudonormal")
  colnames(subset_df) <- new_names

  return(subset_df)
}

create.df.tumors <- function(cur_tumor_set, t.df) {
  # Combine all of the tumors in the set into one table with all signals
  tumors_to_select <- cur_tumor_set$sample.id
  subset_df <- t.df %>% dplyr::select(dplyr::all_of(tumors_to_select))
  return(subset_df)
}

add.pseudonormal.sif <- function(tumor_sif, pseudonormal_df, sif) {
  final_sif <- tumor_sif

  for (pseudonormal in colnames(pseudonormal_df)) {
    original_tumor <- gsub("_pseudonormal", "", pseudonormal)
    original_tumor_row <- sif[sif$sample.id == original_tumor, , drop = FALSE]
    cur_patient <- original_tumor_row$patient.id
    cur_sex <- original_tumor_row$sex

    new_row <- list(sample.id = pseudonormal,
                    patient.id = cur_patient,
                    sex = cur_sex,
                    type = "normal")
    final_sif[nrow(final_sif) + 1, ] <- new_row
  }

  return(final_sif)
}

locus.join <- function(df1, df2) {
  # Inner join on "locus" column
  merged <- dplyr::inner_join(df1, df2, by = "locus")
  # Reorder columns to keep "locus" first
  merged <- merged %>% dplyr::select("locus", tidyselect::everything())
  return(merged)
}

partitions.format.join <- function(listy, t.df) {
  # Combine all of the inidividual tangent outputs and format it like t.df
  combined_df <- purrr::reduce(listy, locus.join)
  combined_df <- combined_df[, colnames(t.df), drop = FALSE]
  cat('\nChecking if the columns of the pseudotangent output match the input tumor signal columns...\n')
  cat(paste0(toString(all(colnames(combined_df) == colnames(t.df)))), '\n')
  return(combined_df)
}
