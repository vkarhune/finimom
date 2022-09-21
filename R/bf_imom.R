bf_imom <- function(beta, se, tau, r){

  opt <- optimise(f = bfg, interval = c(-1, 1), beta = beta, se = se, tau = tau, r = r)
  # loglap <- opt$objective + 0.5*log(2*pi) - 0.5*determinant(bfhessian(beta, se, tau, r))
  # -opt$objective + 0.5*log(2*pi) - 0.5*log(bfhessian(beta, se, tau, r))
  0.5*r*log(tau) - lgamma(r/2) - opt$objective + 0.5*log(2*pi) - 0.5*log(bfhessian(beta, se, tau, r))

}

bf <- Vectorize(bf_imom, vectorize.args = c("beta", "se"))

bfg <- function(x, beta, se, tau, r){
  log(x^2)*(r + 1)/2 - beta*x/(se^2) + x^2/(2*(se^2)) + tau/(x^2)
}

bfhessian <- function(beta, se, tau, r){
  -(r + 1)/beta^2 + 1/se^2 + 6*tau/beta^4
}

