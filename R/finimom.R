

#' Title
#'
#' @param beta
#' @param se
#' @param eaf
#' @param R
#' @param cs
#' @param cs_num
#' @param cs_level
#' @param maxsize
#' @param tau0
#' @param r0
#' @param niter
#' @param burnin
#' @param seed
#' @param excl.burnin
#' @param a0
#' @param b0
#' @param inds0
#' @param standardize
#' @param verbose
#' @param clump
#' @param clump_r2
#' @param check_ld
#' @param pip
#' @param u
#' @param insampleLD
#'
#' @return
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
                    clump = TRUE, clump_r2 = 0.99^2, check_ld = FALSE){

  # all checks here
  if(is.na(beta)) { stop("Effect sizes required") }
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
    if(is.null(u)) { u <- 1.5 }
  } else {
    if(is.null(check_ld)){ check_ld <- TRUE }
    if(is.null(u)) { u <- 1.75 }
  }



  samples <- posterior_samples(
    beta = beta, se = se, eaf = eaf, R = R,
    maxsize = maxsize, tau0 = tau0, r0 = r0, p = p,
    niter = niter, burnin = burnin, seed = seed, excl.burnin = TRUE,
    a0 = a0, b0 = b0, inds0 = inds0, standardize = standardize,
    verbose = verbose, clump = clump, clump_r2 = clump_r2, check_ld = check_ld
    )

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
