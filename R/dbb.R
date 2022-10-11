# beta-binomial prior for model size:
#' Title
#'
#' @param x Function argument.
#' @param p Parameter p.
#' @param a Parameter a.
#' @param b Parameter b.
#'
#' @return Log-density value.
#' @export
#'
#' @examples
dbb <- function(x, p, a, b, log = TRUE){
  out <- lchoose(p, x) + lbeta(x + a, p - x + b) - lbeta(a, b)
  if(!(log)) { out <- exp(out) }
  return(out)
}
