#' Germline Copy Number Variants Blacklist Data Version hg38
#'
#' An updated CNV blacklist, used for GISTIC, derived from a SNP6 array analysis
#' with Cancer Gene Census loci removed
#'
#' @format ## `hg38_blacklist`
#' A data frame with 7,285 rows and 6 columns:
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
#' @format ## `hg19_blacklist`
#' A data frame with 7,285 rows and 6 columns:
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
