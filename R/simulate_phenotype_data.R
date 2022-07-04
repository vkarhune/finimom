#' Title
#'
#' @param X Genotype matrix.
#' @param N Number of individuals.
#' @param p Number of genetic variants.
#' @param causals Indices of causal variants.
#' @param R2 Variance explained by causal variants.
#' @param seed Random seed.
#' @param meanbeta Mean of the Gaussian where effect sizes are drawn from.
#' @param sdbeta Standard deviation of the Gaussian where the effect sizes are drawn from.
#'
#' @return List of the phenotype and the betas, after scaling to var(Y) = 1.
#' @export
#'
#' @examples
simulate_phenotype_data <- function(X, N, p, causals, R2, seed,
                                    meanbeta = 0, sdbeta = 1){

  set.seed(seed)

  # meanbeta <- 0
  # sdbeta <- 1

  betas <- vector("numeric", length = p)
  betas[causals] <- stats::rnorm(length(causals), mean = meanbeta, sd = sdbeta)
  xtb <- X %*% betas

  # R2 <- 0.005

  sigmayy2 <- stats::var(xtb)*(1/R2 - 1)

  YY <- xtb + stats::rnorm(nrow(X), mean = 0, sd = sqrt(sigmayy2))

  # scale back to sd(Y) = 1
  Y <- YY/stats::sd(YY)
  betas <- betas/stats::sd(YY)
  # sigma2 <- sigmayy2/var(YY)

  return(list(Y, betas))
}
