# sampling from vectors of length 1
# see "?sample"
#' Title
#'
#' @param x Vector.
#' @param ... Other arguments to sample.
#'
#' @return Value.
#' @export
#'
#' @examples
sample2 <- function(x, ...){ x[sample.int(length(x), ...)] }
