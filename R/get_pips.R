#' Obtain posterior inclusion probabilities from the posterior samples.
#'
#' @param samples Posterior samples.
#'
#' @return Posterior inclusion probabilities
#' @export
#'
#' @examples
get_pips <- function(samples){

    if(is.na(samples[[1]])) { stop(
      "Make sure you are using the correct object.\nIf you are sure this is the correct object, use the main function with option 'pip = TRUE'.")
      }



  pip_clusters <- (colSums(samples[[1]] != 0)/nrow(samples[[1]]))

  out <- do.call("rbind", lapply(seq_along(pip_clusters), function(x) {
    cbind(samples[[5]][[x]], pip_clusters[x]/length(samples[[5]][[x]]))
  }))

  out <- out[order(out[,1]),]
}
