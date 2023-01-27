

#' Fine-mapping using an inverse-moment prior.
#'
#' @param beta Effect size estimates.
#' @param se Standard errors.
#' @param eaf Effect allele frequencies.
#' @param R LD-matrix.
#' @param cs Are credible sets returned?
#' @param cs_num How many credible sets? The default is to use the posterior mode.
#' @param cs_level Credible set level.
#' @param maxsize Maximum model size.
#' @param tau0 Prior parameter tau.
#' @param r0 Prior parameter r.
#' @param niter Number of iterations.
#' @param burnin Number of burn-in iterations.
#' @param seed Random seed.
#' @param excl.burnin Should burn-in be excluded?
#' @param a0 Hyperparameter a for the model size prior.
#' @param b0 Hyperparameter b for the model size prior (but preferably use parameter u).
#' @param inds0 Indices for the starting model.
#' @param standardize Should the effect sizes be standardised? Defaults to TRUE.
#' @param verbose Verbose output.
#' @param clump Should clumping be done for extremely highly-correlated variants?
#' @param clump_r2 Clumping threshold for extremely highly-correlated variants.
#' @param check_ld Should LD discrepancy check be performed?
#' @param pip Are posterior inclusion probabilities returned?
#' @param u Hyperparameter for model size prior. Defaults to 1.5 for in-sample LD matrix and 1.75 for out-of-sample LD matrix.
#' @param insampleLD Is in-sample LD used?
#' @param ala Whether Approximate Laplace should be used?
#'
#' @return List.
#' @export
#'
#' @examples
finimom <- function(beta, se, eaf, R,
                    cs = TRUE,
                    cs_num = NULL,
                    cs_level = 0.95,
                    pip = FALSE,
                    maxsize = 10, tau0 = 0.0083, r0 = 1, niter = 12500, burnin = 2500, seed = 456, excl.burnin = TRUE,
                    a0 = 1, b0 = NULL, u = NULL, inds0 = NULL, standardize = TRUE,
                    verbose = TRUE,
                    insampleLD = NULL,
                    clump = TRUE, clump_r2 = 0.99^2, check_ld = FALSE,
                    ala = NULL){

  # all checks here
  if(is.null(beta)) { stop("Effect sizes required") }
  p <- length(beta)

  if(length(se) != p){ stop("beta and se are of different length") }

  if(!(is.null(eaf)) & (length(eaf) != p)){ stop("beta and eaf are of different length") }

  R <- as.matrix(R)

  if(dim(R)[1] != length(beta) | dim(R)[2] != length(beta)) { stop("Dimension of R does not match with beta") }

  if(is.null(insampleLD)) { stop("Determine whether in-sample LD is used") }

  if(is.null(eaf) & standardize){
    warning("No allele frequencies provided - changing to 'standardize = FALSE'")
    standardize <- FALSE
  }

  if(!(is.null(b0)) & !(is.null(u))) { warning("Both b0 and u given - using p^u for beta-binomial hyperparameter") }

  # define input params to workhorse function

  if(!(is.null(u)) | is.null(b0)){
    b0 <- p^u
  }

  if(insampleLD){
    if(is.null(check_ld)){ check_ld <- FALSE }
    if(is.null(u)) {
      u <- 1.5
      b0 <- p^u
    }
  } else {
    if(is.null(check_ld)){ check_ld <- TRUE }
    if(is.null(u)) {
      u <- 1.75
      b0 <- p^u
    }
  }


if(is.null(ala) | ala == FALSE){
  samples <- posterior_samples(
    beta = beta, se = se, eaf = eaf, R = R,
    maxsize = maxsize, tau0 = tau0, r0 = r0, p = p,
    niter = niter, burnin = burnin, seed = seed, excl.burnin = TRUE,
    a0 = a0, b0 = b0, inds0 = inds0, standardize = standardize,
    verbose = verbose, clump = clump, clump_r2 = clump_r2, check_ld = check_ld
    )
} else {
  samples0 <- finimom::posterior(
    dat = list(beta = beta, se = se, LDmat = R),
    tau = rep(tau0, length(beta)),
    r = r0, p = p, u = u,
    lpriorval = lprior,
    niter = niter,
    maxsize = maxsize)
  samples <- samples0[c("betavecmat","modsize","modindices","value")]
  samples[["betavecmat"]] <- t(samples[["betavecmat"]])
  samples[["modsize"]] <- as.numeric(samples[["modsize"]])
  samples[["value"]] <- as.numeric(samples[["value"]])
  # names(samples) <- NULL
  samples <- c(samples, samples0[5])
  samples[[3]] <- sapply(samples[[3]], function(x){gsub(" ", ",", x)})

  if(excl.burnin){
    samples <- list(samples[[1]][(burnin + 1):niter,],
                    samples[[2]][(burnin + 1):niter],
                    samples[[3]][(burnin + 1):niter],
                    samples[[4]][(burnin + 1):niter])
  }
}

  ppcs <- prop.table(table(samples[[2]]))

  cs_best <- as.numeric(names(which.max(ppcs)))

  if(is.null(cs_num)) {
    cs_num <- cs_best
  } else {
    if(cs_num != cs_best) { warning("The requested number of credible sets is not the same as the posterior mode") }
  }

  out <- list("samples" = samples)

  if(cs){
    sets <- get_credible_sets(samples = samples, num_signals = cs_num, level = cs_level)
    out <- c(out, "sets" = list(sets))
  }

  if(pip){
    pip <- get_pips(samples = samples)
    out <- c(out, "pip" = list(pip[,2]))
  }

  # out <- list(samples, sets)

  return(out)

}
