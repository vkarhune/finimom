create_Cmatrix <- function(rho){

  if(length(rho) == 1){
    edrho <- eigen(matrix(c(1, rho, rho, 1), nrow = 2))
  } else {
    edrho <- eigen(rho)
  }
  S <- edrho$vectors %*% diag(edrho$values^(-0.5)) %*% t(edrho$vectors)

  out <- stats::cov2cor(S)
  return(out)
  #return(S)
}
