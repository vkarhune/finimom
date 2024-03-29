#' Check discrepancy of the linkage equilibrium and test statistics
#'
#' @param z Test statistic.
#' @param Chi2_quantile Quantile for declaring a mismatch. The default is the median of a Chi-squared distribution with one degree of freedom.
#' @param LDm LD matrix.
#' @param indices Variant indices.
#' @param clump_r2 Clumping threshold for extremely highly correlated variants.
#'
#' @return List of indices within each LD clump, corrected for the LD discrepancy.
#' @export
#'
#' @examples
check_ld_disc <- function(indices, z, Chi2_quantile = NULL, LDm, clump_r2){

  if(length(indices) == 1){
    ld_disc_out <- indices
  } else {


  z <- z[indices]
  LDmat <- LDm[indices, indices]

  if(abs(z[1]) < 3){
    ld_disc_out <- indices
  } else {


  maxz <- 1 # this is sorted in "clump_variants.R"

  # observed
  observed <- z

  # expected
  expected <- z[maxz]*LDmat[maxz,]

  # test statistic
  # teststat <- ((observed - expected)/LDmat[maxz,])^2

  # teststat <- (observed - expected)^2

  # X - Y \sim N(mux - muy, sigmax2 + sigmay2 - 2*cor(X,y)*sigmax*sigmay)
  # Z = ((X - Y) - mu[x-y])/sigma[x-y] \sim N(0, 1)
  # https://srabbani.com/bivariate.pdf
  # Hormozdiari et al.
  # teststat <- (observed - expected)/sqrt(1 + 1 - 2*LDmat[maxz,])
  # this reduces to 0 if r = 1 - smth to protect this
  # solution: use clump_r2 threshold
  # rr <- ifelse(LDmat[maxz,] > sqrt(clump_r2), sqrt(clump_r2), LDmat[maxz,])
  teststat <- (observed - expected)^2/(2*(1 - sqrt(clump_r2)))



  if(is.null(Chi2_quantile)) { Chi2_quantile <- 0.5 }

  mismatch <- which(teststat > qchisq(Chi2_quantile, df = 1))

  if(length(mismatch) > 0){

    # projected_r <- mean(c(LDmat[mismatch,maxz], sqrt(observed[mismatch]/observed[maxz])*sign(LDmat[mismatch,maxz])))

    ratio <- observed[mismatch]/observed[maxz]
    sqratio <- ifelse(ratio < 0, abs(ratio), ratio)

    projected_r <- sqrt(sqratio)*sign(LDmat[mismatch,maxz])
    projected_r <- ifelse(ratio < 0, -1*projected_r, projected_r)

    LDmat[mismatch, maxz] <- projected_r
    LDmat[maxz, mismatch] <- projected_r

    ld_disc_out <- clump_variants(R = LDmat, clump_r2 = clump_r2, z = z)

    ld_disc_out <- lapply(ld_disc_out, function(x) indices[x])

    # ld_disc_out <- list(ld_disc_out, c(indices[maxz], indices[mismatch], projected_r))

  } else {
    ld_disc_out <- indices
  }

  }
}



  return(ld_disc_out)
}
