g <- function(x, beta, se, tau, psi, r, k, R){

  if(k == 1){
    S <- se
    invS <- 1/S
    #invR <- 1/R
  } else {
    S <- diag(se)
    invS <- solve(S)
    #invR <- solve(R)
  }

  z <- as.vector(beta/se)

  -k*(0.5*r*sum(log(tau*psi)) - sum(log(gamma(0.5*r)))) + 0.5*(r + 1)*sum(log(x^2)) +
    # 0.5*n*determinant(2*pi * S %*% R %*% S, logarithm = TRUE)$modulus +
    sum(tau*psi/(x^2)) + 0.5*( - t(z) %*% invS %*% x - t(x) %*% invS %*% z + # + t(z) %*% invR %*% z
                                 t(x) %*% invS %*% R %*% invS %*% x) # - logmk(lprior, k)
}
