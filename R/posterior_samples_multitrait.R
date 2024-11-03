#' Samples from the posterior distribution
#'
#' @param beta Vector of effect sizes.
#' @param se Vector of standard errors.
#' @param eaf Vector of effect allele frequencies.
#' @param R LD matrix.
#' @param maxsize The maximum number of causal variants.
#' @param tau0 Parameter tau.
#' @param r0 Parameter $r$.
#' @param niter Number of iterations.
#' @param burnin Number of burn-in samples.
#' @param p Number of variants.
#' @param seed Random seed.
#' @param excl.burnin Should the burn-in be excluded?
#' @param a0 Hyperparameter a for the model size prior.
#' @param b0 Hyperparameter b for the model size prior.
#' @param inds0 Initial model indices (not used).
#' @param standardize Should the effect sizes be standardised? Defaults to TRUE.
#' @param verbose Verbose output.
#' @param clump Whether to clump extremely highly correlated variants.
#' @param clump_r2 Clumping threshold for extremely highly correlated variants.
#' @param check_ld Should the Z-scores be checked for LD discrepancy?
#' @param ala Whether Approximate Laplace should be used?
#'
#' @return List.
#' @export
#'
#' @examples
posterior_samples_multitrait <- function(
    beta, se, eaf, R, k, omega, u, n,
    maxsize, tau0, r0, niter, burnin, p, seed = 456, excl.burnin = TRUE,
    a0 = 1, b0 = NULL, inds0 = NULL, standardize = TRUE,
    verbose = TRUE, clump = TRUE, clump_r2 = 0.99^2, check_ld = FALSE,
    ala = NULL){



  if(standardize){
    beta <- lapply(beta, function(x) x*sqrt(2*eaf*(1-eaf)))
    se <- lapply(se, function(x) x*sqrt(2*eaf*(1-eaf)))
  }

  z <- lapply(seq_len(k), function(x) beta[[x]]/se[[x]])

  if(k == 2){
  z1c <- z[[1]][apply(do.call("cbind", lapply(z, abs)), 1, max) < 2]
  z2c <- z[[2]][apply(do.call("cbind", lapply(z, abs)), 1, max) < 2]

  rho_null <- cor(z1c, z2c)
  Cmat <- matrix(c(1, -rho_null, -rho_null, 1), nrow = 2)

  } else {
    Cmat <- diag(k)
  }

  if(clump){
    cat(sprintf("Clumping variants at r2=%.3g\n", clump_r2))

    zclump <- apply(do.call("cbind", lapply(z, abs)), 1, max)

    keeplist_cleaned <- clump_variants(R = R, clump_r2 = clump_r2, z = zclump)

    if(check_ld){
      ldcheck <- lapply(keeplist_cleaned, check_ld_disc, z = zclump,
                        # ldcheck <- lapply(keeplist_cleaned[1:5], check_ld_disc, z = z,
                        Chi2_quantile = 0.5,
                        LDm = R, clump_r2 = clump_r2)

      keeplist_cleaned <- lapply(rapply(ldcheck, enquote, how = "unlist"), eval)


    }

    keepinds <- sapply(keeplist_cleaned, "[", 1)



    beta <- lapply(beta, function(x) x[keepinds])
    se <- lapply(se, function(x) x[keepinds])
    z <- lapply(z, function(x) x[keepinds])
    p <- length(keepinds)

    if(is.null(omega)){
      omega <- rep(1/p, times = length(unlist(beta)))
    }

    R <- R[keepinds, keepinds]

    if(!(is.matrix(R))){
      cat(sprintf("Only one set of variants for fine-mapping - exiting\n"))
      return(NULL)
    }

    # sum(abs(R) > sqrt(clump_r2)) == length(keepinds)


  }

  if(maxsize > p){
    cat(sprintf("Note: maximum model size set to %i\n", p))
    maxsize <- p
  }

  msprior <- "complexity"

  if(msprior %in% "complexity") {
    lprior <- sapply(seq_len(maxsize), dbb, p = p, a = a0, b = p^u)
    lprior <- log(exp(lprior)/sum(exp(lprior)))
    lprior <- rep(lprior, k)
  }

  if(is.null(ala)) ala <- TRUE
  set.seed(seed)

  dat <- list(
    beta = unlist(beta),
    se = unlist(se),
    LDmat = Cmat %x% R,
    npheno = k
  )

  if(length(n) == 1){ n <- rep(n, k) }
  taus <- sapply(seq_along(n), function(x) estimate_tau(n = n[x]))


    if(verbose){ cat(sprintf("Sampling from the posterior...\n")) }
    prc <- proc.time()

    out <- posteriormv(dat = dat,
                       tau = c(rep(taus, each = p)),
                       maxsize = maxsize, r = 1, p = length(dat[["beta"]]),
                       niter = niter,
                       lpriorval = lprior,
                       approx = 1, k = 2,
                       omega = omega)

    cat(sprintf("\n%i iterations done in %.2f seconds\n", niter, (proc.time() - prc)[[3]]))

  if(excl.burnin){
    out <- list(out[[1]][(burnin + 1):niter,],
                out[[2]][,(burnin + 1):niter],
                out[[3]][(burnin + 1):niter])
  }

  if(clump){
    out <- c(out,
             list(keeplist_cleaned))
  } else {
    out <- c(out, list(seq_len(p)))
  }

  return(out)

}
