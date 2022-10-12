#' Title
#'
#' @param x
#' @param r
#' @param tau
#' @param phi
#'
#' @return
#' @export
#'
#' @examples
#'
dimom <- Vectorize(function(x, r, tau, phi = 1){
  if(x == 0){
    0
  } else {
    ((tau*phi)^(r/2)/gamma(r/2))*abs(x)^(-(r+1))*exp(-tau*phi/x^2)
  }
})
