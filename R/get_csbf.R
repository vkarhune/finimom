#' Title
#'
#' @param beta
#' @param se
#' @param tau
#' @param r
#' @param level
#'
#' @return
#' @export
#'
#' @examples
get_csbf <- function(beta, se, tau, r, level = 0.95){
  lbf <- finimom:::bf(beta = beta, se = se, tau = tau, r = r)

  pp <- exp(lbf)/sum(exp(lbf))

  credset <- order(pp, decreasing = T)[which(cumsum(sort(pp, decreasing = TRUE)/sum(pp)) < level)]

  if(length(credset) == 0) credset <- which.max(pp)

  credset

}
