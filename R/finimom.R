

#' Fine-mapping using an inverse-moment prior.
#'
#' @param beta Effect size estimates.
#' @param se Standard errors.
#' @param eaf Effect allele frequencies.
#' @param R LD-matrix.
#' @param n GWAS sample size. If provided, then tau is calculated based on n.
#' @param cs Are credible sets returned?
#' @param cs_num How many credible sets? The default is to use the posterior mode.
#' @param cs_level Credible set level.
#' @param maxsize Maximum model size.
#' @param tau Prior parameter tau.
#' @param r Prior parameter r.
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
#' @param purity Credible set purity, defined as the minimum absolute correlation between the variants in a credible set.
#'
#' @return List.
#' @export
#'
#' @examples
finimom <- function(beta, se, eaf, R,
                    n = NULL,
                    cs = TRUE,
                    cs_num = NULL,
                    cs_level = 0.95,
                    pip = FALSE,
                    maxsize = 10, tau = NULL, r = 1, niter = 12500, burnin = 2500, seed = 456, excl.burnin = TRUE,
                    a0 = 1, b0 = NULL, u = NULL, inds0 = NULL, standardize = TRUE,
                    verbose = TRUE,
                    insampleLD = NULL,
                    clump = TRUE, clump_r2 = 0.99^2, check_ld = FALSE,
                    ala = NULL, purity = NULL){

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

  if(is.null(n) & is.null(tau)) { stop("Please provide n to calculate tau, or tau directly.") }

  if(!(is.null(n)) & is.null(tau)) {
    cat("Calculating tau based on the sample size.\n")

    tau <- estimate_tau(n = n)
  }

  if(!(is.null(b0)) & !(is.null(u))) { warning("Both b0 and u given - using p^u for beta-binomial hyperparameter") }

  default_u_insampleLD <- 2
  default_u_refLD <- 2.25

  # define input params to workhorse function

  if(!(is.null(u)) | is.null(b0)){
    b0 <- p^u
  }

  if(insampleLD){
    if(is.null(check_ld)){ check_ld <- FALSE }
    if(is.null(u)) {
      u <- default_u_insampleLD
      b0 <- p^u
    }
  } else {
    if(is.null(check_ld)){ check_ld <- TRUE }
    if(is.null(u)) {
      u <- default_u_refLD
      b0 <- p^u
    }
  }

if(is.null(ala)) ala <- FALSE

samples <- posterior_samples(
  beta = beta, se = se, eaf = eaf, R = R,
  maxsize = maxsize, tau0 = tau, r0 = r, p = p,
  niter = niter, burnin = burnin, seed = seed, excl.burnin = TRUE,
  a0 = a0, b0 = b0, inds0 = inds0, standardize = standardize,
  verbose = verbose, clump = clump, clump_r2 = clump_r2, check_ld = check_ld,
  ala = ala
)

  ppcs <- prop.table(table(samples[[2]]))

  cs_best <- as.numeric(names(which.max(ppcs)))

  if(is.null(cs_num)) {
    cs_num <- cs_best
  } else {
    if(cs_num != cs_best) { warning("The requested number of credible sets is not the same as the posterior mode") }
  }

  out <- list("samples" = samples, "signals" = ppcs)

  if(cs){
    sets <- get_credible_sets(samples = samples, num_signals = cs_num, level = cs_level, purity = purity, R = R)
    out <- c(out, "sets" = list(sets))
  }

  if(pip){
    pip <- get_pips(samples = samples)
    out <- c(out, "pip" = list(pip[,2]))
  }

  # out <- list(samples, sets)

  return(out)

}
