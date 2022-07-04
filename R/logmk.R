#' Title
#'
#' @param lprior Output of dbb.
#' @param k Model size.
#'
#' @return Log-value of model size prior.
#' @export
#'
#' @examples
logmk <- function(lprior, k){ lprior[k] }
