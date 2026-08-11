kneedle_points <- function(betalist, eR, h2leadlist, n_eigen = n_eigen, verbose = verbose){

  out <- lapply(seq_along(betalist), function(i){

    betavec <- betalist[[i]]
    h2lead <- h2leadlist[[i]]

    quadterm <- quadratic(b = betavec, eR = eR, q = n_eigen)

    kn <- kneedler(y = quadterm, h2lead = h2lead, Sens = 1)
    outk <- kn$knee_y

    if(verbose){
      cat(sprintf("Trait %i:\n", i))
      cat(sprintf("Kneedle point: %s\n",
                  formatC(signif(out, digits = 3), digits = 3, format = "fg", flag = "#")))
      cat(sprintf("Lead variant R2: %s\n",
                  formatC(signif(h2lead, digits = 3), digits = 3, format = "fg", flag = "#")))
    }

    return(outk)

  })

  return(out)

}
