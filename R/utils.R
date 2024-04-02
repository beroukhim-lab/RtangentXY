### Author: Siyun Lee siyun@broadinstitute.org, Kei Enomoto kenomoto@broadinstitute.org, Greg Raskind graskind@g.harvard.edu
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: March 19, 2024
### License: GNU GPL >= 2, Copyright (C) 2024 Dana-Farber Cancer Institute
### Dependencies: Check DESCRIPTION file
### See https://github.com/beroukhim-lab/RtangentXY for more information

#' @import utils

# Note, all of the functions in this file are used internally

# We need to add global variables, these refer to the raw inputs column names
utils::globalVariables(c("gender", "locus", "sample.id", "patient.id", "type", "start", "end", "arm"))

# We need to add global variables used to generate plots
utils::globalVariables(c("signal", "d", "chr", "arm.median", "lf", "metric", "value"))

# We need to add a global variable for ".", as this is used in dplyr
utils::globalVariables(c("."))
