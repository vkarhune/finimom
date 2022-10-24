#' Log-marginal likelihood for the data.
#'
#' @param x Vector of effect sizes.
#' @param se Vector of standard errors.
#' @param tau Parameter tau.
#' @param psi Parameter psi (not used).
#' @param r Parameter r.
#' @param k Model size.
#' @param R LD matrix.
#' @param gval Value of g.
#'
#' @return Value.
#' @export
#'
#' @examples
logmyb <- function(x, se, tau, psi, r, k, R, gval){

  if(k == 1){
    S <- se
    invS <- 1/S
  } else {
    S <- diag(se)
    invS <- solve(S)
  }


  hessian <- diag(6*tau*psi/x^4 - (r + 1)/x^2) + invS %*% R %*% invS

  logdeth <- determinant(hessian, logarithm = TRUE)$modulus

  -gval + 0.5*k*log(2*pi) - 0.5*logdeth

}
