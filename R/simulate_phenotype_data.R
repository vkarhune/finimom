simulate_phenotype_data <- function(X, N, p, causals, R2, seed,
                                    meanbeta = 0, sdbeta = 1){

  set.seed(seed)

  # meanbeta <- 0
  # sdbeta <- 1

  betas <- vector("numeric", length = p)
  betas[causals] <- rnorm(length(causals), mean = meanbeta, sd = sdbeta)
  xtb <- X %*% betas

  # R2 <- 0.005

  sigmayy2 <- var(xtb)*(1/R2 - 1)

  YY <- xtb + rnorm(nrow(X), mean = 0, sd = sqrt(sigmayy2))

  # scale back to sd(Y) = 1
  Y <- YY/sd(YY)
  betas <- betas/sd(YY)
  # sigma2 <- sigmayy2/var(YY)

  return(list(Y, betas))
}
