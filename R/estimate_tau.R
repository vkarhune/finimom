#' Estimate tau based on sample size and z score
#'
#' @param n Sample size in GWAS summary statistics
#' @param pval P-value threshold for detectable effect size (z-score used if both p-value and z-score provided)
#' @param z Z-score for detectable effect size (z-score used if both p-value and z-score provided)
#' @param q Quantile for the proposed Z score
#'
#' @return
#' @export
#'
#' @examples
estimate_tau <- function(n = NULL, q = 0.05, pval = 0.001, z = NULL){
  if(!(is.null(pval)) & !(is.null(z))) {
    warning(sprintf("Both p-value and z-score provided; using z=%s to calculate tau\n", as.character(z)))
  } else if(!(is.null(pval)) & is.null(z)) {
    z <- qnorm(1-pval/2)
  }

  bq <- z/sqrt(n)
  qgamma(1 - q, shape = 0.5, scale = 1)*bq^2
}
