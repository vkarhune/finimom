create_Cmatrix <- function(rho){
 edrho <- eigen(matrix(c(1, rho, rho, 1), nrow = 2))
 S <- edrho$vectors %*% diag(edrho$values^(-0.5)) %*% t(edrho$vectors)
 out <- stats::cov2cor(S)
 return(out)
}
