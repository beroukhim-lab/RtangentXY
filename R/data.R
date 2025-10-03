#' Germline Copy Number Variants Blacklist Data Version hg38
#'
#' An updated CNV blacklist, used for GISTIC, derived from a SNP6 array analysis
#' with Cancer Gene Census loci removed
#'
#' @format A data frame with 7,285 rows and 6 columns:
#' \describe{
#'   \item{id}{Unique ID for CNV}
#'   \item{chr}{Chromosome of CNV}
#'   \item{start}{CNV start location}
#'   \item{end}{CNV end location}
#'   \item{flanking.start}{Start location of region flanking CNV}
#'   \item{flanking.end}{End location of region flanking CNV}
#' }
#' @source <https://github.com/beroukhim-lab/gistic_reference/tree/main>
"hg38_blacklist"

#' Germline Copy Number Variants Blacklist Data Version hg19
#'
#' A CNV blacklist, used for GISTIC, derived from a SNP6 array analysis
#'
#' @format A data frame with 7,454 rows and 6 columns:
#' \describe{
#'   \item{id}{Unique ID for CNV}
#'   \item{chr}{Chromosome of CNV}
#'   \item{start}{CNV start location}
#'   \item{end}{CNV end location}
#'   \item{flanking.start}{Start location of region flanking CNV}
#'   \item{flanking.end}{End location of region flanking CNV}
#' }
#' @source <ftp://ftp.broadinstitute.org/pub/GISTIC2.0/hg19_support/CNV.hg19.bypos.111213.txt>
"hg19_blacklist"

#' Example Sample Information File
#'
#' A tibble containing sample information for a small test dataset
#' 
#' @format A tibble with 4 columns:
#'  \describe{
#'    \item{sample.id}{Sample ID}
#'    \item{patient.id}{Patient ID}
#'    \item{sex}{One of "male" or "female"}
#'    \item{type}{One of "normal" or "tumor"}
#' }
"example_sif"

#' Example Probe Information File
#' 
#' A tibble containing probe information for a small test dataset
#' 
#' @format A tibble with 5 columns:
#'  \describe{
#'    \item{locus}{Locus in the format `"{chr}:{start}-{end}"` where `{chr}` is 1-22, X, or Y}
#'    \item{chr}{Chromosome of the probe}
#'    \item{start}{Start position of the probe}
#'    \item{end}{End position of the probe}
#'    \item{arm}{Arm of the chromosome (p or q)}
#' }
"example_pif"

#' Example normal signal matrix file
#'
#' A tibble containing log2 relative copy number values for normal samples in a small test dataset
#' 
#' @format A tibble with 4,800 rows and 11 columns:
#'  \describe{
#'    \item{locus}{Locus in the format `"{chr}:{start}-{end}"` where `{chr}` is 1-22, X, or Y}
#'    \item{sample.id}{One column for each normal sample containing log2(RCN) values}
#' }
"example_nsig_df"

#' Example tumor signal matrix file
#' 
#' A tibble containing log2 relative copy number values for tumor samples in a small test dataset
#' 
#' @format Analagous to `example_nsig_df`
"example_tsig_df"