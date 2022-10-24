#' Extract a subchain from exsisting posterior samples, excluding a burn-in.
#'
#' @param posterior_sample Object from posterior_samples().
#' @param burnin The number of initial iterations to be excluded.
#' @param niter The number of iterations kept.
#'
#' @return Posterior samples with burn-in excluded.
#' @export
#'
#' @examples
extract_subchain <- function(posterior_sample, burnin = NULL, niter = NULL){
  chainstart <- burnin + 1
  chainstop <- burnin + niter

  list(posterior_sample[[1]][chainstart:chainstop, ],
       posterior_sample[[2]][chainstart:chainstop],
       posterior_sample[[3]][chainstart:chainstop],
       posterior_sample[[4]][chainstart:chainstop],
       posterior_sample[[5]])
}
