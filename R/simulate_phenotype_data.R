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
#' @param minpower
#' @param alpha
#'
#' @return List of the phenotype and the betas, after scaling to var(Y) = 1.
#' @export
#'
#' @examples
simulate_phenotype_data <- function(X, N, p, causals, R2, seed,
                                    meanbeta = 0, sdbeta = 1, minpower = NULL, alpha = NULL){

  set.seed(seed)

  # meanbeta <- 0
  # sdbeta <- 1

  betas <- vector("numeric", length = p)

  if(!(is.null(minpower))){
    power <- minpower
    # alpha <- 0.05
    bq <- (qnorm(1-alpha/2) + qnorm(power))/sqrt(N)
  }

  powercheck <- FALSE
  while(!powercheck){


  betas[causals] <- stats::rnorm(length(causals), mean = meanbeta, sd = sdbeta)


  xtb <- X %*% betas

  # R2 <- 0.005

  sigmayy2 <- stats::var(xtb)*(1/R2 - 1)

  YY <- xtb + stats::rnorm(nrow(X), mean = 0, sd = sqrt(sigmayy2))

  betas <- betas/stats::sd(YY)

  if(is.null(minpower) | all(abs(betas[causals]) > bq)){
    powercheck <- TRUE
  }
  }

  # scale back to sd(Y) = 1
  Y <- YY/stats::sd(YY)
  # sigma2 <- sigmayy2/var(YY)

  return(list(Y, betas))
}
