h2_cap <- function(betalist, LDmat, n_eigen, h2leads, verbose = verbose){

  cond <- FALSE

  while(!cond){

    eR <- RSpectra::eigs_sym(LDmat, k = n_eigen, which = "LM")

    kneedles <- tryCatch(
      kneedle_points(betalist = betalist, eR = eR, h2leadlist = h2leads,
                     n_eigen = n_eigen, verbose = verbose),
      error = function(error) { return(NULL) })


    if(!(is.null(kneedles))){
      cond <- TRUE
    } else {
      n_eigen <- n_eigen*2
    }
  }

  if(verbose){ cat(sprintf("Eigendecomposition nu = %i\n", n_eigen)) }

  return(kneedles)
}

