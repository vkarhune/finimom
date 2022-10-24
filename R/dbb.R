# beta-binomial prior for model size:
#' Beta-binomial density.
#'
#' @param x Function argument.
#' @param p Parameter p.
#' @param a Parameter a.
#' @param b Parameter b.
#' @param log Whether log-density is returned -- defaults to TRUE.
#'
#' @return (Log-)density value.
#' @export
#'
#' @examples
dbb <- function(x, p, a, b, log = TRUE){
  out <- lchoose(p, x) + lbeta(x + a, p - x + b) - lbeta(a, b)
  if(!(log)) { out <- exp(out) }
  return(out)
}
