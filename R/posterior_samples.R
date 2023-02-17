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
posterior_samples <- function(
    beta, se, eaf, R, maxsize, tau0, r0, niter, burnin, p, seed = 456, excl.burnin = TRUE,
    a0 = 1, b0 = NULL, inds0 = NULL, standardize = TRUE,
    verbose = TRUE, clump = TRUE, clump_r2 = 0.99^2, check_ld = FALSE,
    ala = NULL){


  if(0){
    R <- LDmat
    beta <- summarystats[,1]
    se <- summarystats[,2]
    eaf <- eafs
    maxsize <- 10
    tau0 <- 0.0083
    r0 <- 1
    niter <- 2500
    burnin <- 500
    #p <- nrow(d)
    seed <- 456
    #n <- 1000
    p <- length(beta)
    u <- 1.5
    a0 <- 1
    b0 <- p^u
    inds0 <- NULL
    clump_r2 <- 0.99^2
  }

  if(is.null(b0)) { b0 <- p^1.5 }

  if(standardize){
    beta <- beta*sqrt(2*eaf*(1-eaf))
    se <- se*sqrt(2*eaf*(1-eaf))
  }

  z <- beta/se

  if(clump){
    cat(sprintf("Clumping variants at r2=%.3g\n", clump_r2))

    keeplist_cleaned <- clump_variants(R = R, clump_r2 = clump_r2, z = z)

    if(check_ld){
      ldcheck <- lapply(keeplist_cleaned, check_ld_disc, z = z,
                        # ldcheck <- lapply(keeplist_cleaned[1:5], check_ld_disc, z = z,
                        Chi2_quantile = 0.5,
                        LDm = R, clump_r2 = clump_r2)

      keeplist_cleaned <- lapply(rapply(ldcheck, enquote, how = "unlist"), eval)


    }

    keepinds <- sapply(keeplist_cleaned, "[", 1)



    beta <- beta[keepinds]
    se <- se[keepinds]
    z <- z[keepinds]
    p <- length(keepinds)

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
    lprior <- sapply(seq_len(maxsize), dbb, p = p, a = a0, b = b0)
    lprior <- log(exp(lprior)/sum(exp(lprior)))
  }

  # initialise
  if(is.null(inds0)){
    inds <- which.max(abs(z))
  } else {
    inds <- inds0
  }

  if(is.null(ala)) ala <- FALSE

  if(TRUE){
    set.seed(seed)

    #if(verbose){ cat(sprintf("Sampling from the posterior using approximate Laplace...\n")) }
    if(verbose){ cat(sprintf("Sampling from the posterior...\n")) }
    prc <- proc.time()

    out <- posterior(dat = list(beta = beta,
                                se = se,
                                LDmat = R),
                     tau = rep(tau0, length(beta)),
                     r = r0,
                     p = p,
                     lpriorval = lprior,
                     niter = niter,
                     maxsize = maxsize,
                     approx = ala*1)

    out <- out[c("betavecmat","modsize","modindices","value")]
    out[[1]] <- t(out[[1]])
    out[[2]] <- as.numeric(out[[2]])
    out[[3]] <- sapply(out[[3]], function(x){gsub(" ", ",", x)})
    out[[4]] <- as.numeric(out[[4]])

    cat(sprintf("\n%i iterations done in %.2f seconds\n", niter, (proc.time() - prc)[[3]]))

  } else {



  opt <- optim(as.matrix(beta[inds], ncol = 1), g,
                      beta = as.matrix(beta[inds], ncol = 1),
                      se = se[inds],
                      R = R[inds,inds,drop = F],
                      psi = 1,
                      k = length(inds), tau = tau0, r = r0,
                      # method = "BFGS")
                      method = "Nelder-Mead", control = list(warn.1d.NelderMead = FALSE))

  # lml <- logmyb(as.matrix(beta[inds], ncol = 1),
  lml <- logmyb(as.matrix(opt$par, ncol = 1),
                se[inds],
                tau = tau0, psi = 1, r = r0, k = length(inds),
                R = R[inds,inds,drop = F], gval = opt$value)
  lp <- lml + logmk(lprior, length(inds))

  modelsize <- length(inds)

  betavec <- rep(0, p)

  betavec[inds] <- opt$par

  set.seed(seed)

  out <- list(matrix(NA, nrow = niter, ncol = p),
              vector("numeric", length = niter),
              vector("character", length = niter),
              vector("numeric", length = niter)
  )

  pb <- txtProgressBar(min = 1, max = niter, style = 3)

  if(verbose){ cat(sprintf("Sampling from the posterior...\n")) }

  # loop
  prc <- proc.time()
  for(i in 1:niter){
    # for(i in 1:27){

    # i <- 1

    # i <- i + 1


      # randomly pick add = 1, delete = 0, swap = 2
      if(modelsize == 1){
        add <- sample(1:2, size = 1)
      } else if(modelsize == maxsize){
        if(modelsize == p){
          add <- 0
        } else {
          add <- sample(c(0, 2), size = 1)
        }
      } else {
        add <- sample(0:2, size = 1)
      }

      if(add == 1){
        # swapindex <- sample(which(betavec == 0), 1)
        # probs <- abs(z^2)[which(betavec == 0)]/sum(abs(z^2)[which(betavec == 0)])
        xtr <- beta - R %*% betavec
        probs <- (xtr[which(betavec == 0)]^2)/sum((xtr[which(betavec == 0)])^2)
        swapindex <- sample2(which(betavec == 0), size = 1, prob = probs)
        indsprop <- sort(c(inds, swapindex))
      } else if(add == 2){
        swapindex <- sample2(which(betavec != 0), 1)
        probs <- R[swapindex,which(betavec == 0)]^2/sum(R[swapindex,which(betavec == 0)]^2)
        # probs <- abs(R[swapindex,which(betavec == 0)])/sum(abs(R[swapindex,which(betavec == 0)]))
        addindex <- sample2(which(betavec == 0), size = 1, prob = probs)
        indsprop <- sort(c(addindex, setdiff(inds, swapindex)))
      } else {
        swapindex <- sample2(which(betavec != 0), 1)
        indsprop <- setdiff(which(betavec != 0), swapindex)
      }

      cond <- qr(R[indsprop,indsprop])$rank == length(indsprop)


    # inds <- sort(c(inds, newinds))

      if(!(cond)){
        lpnew <- Inf
      } else {

        initpars <- as.matrix(beta[indsprop], ncol = 1)
        if(any(initpars == 0)){ initpars[initpars == 0] <- 1e-7 }

    opt <- optim(initpars,
                        g,
                        beta = as.matrix(beta[indsprop], ncol = 1),
                        se = se[indsprop],
                        R = R[c(indsprop),c(indsprop),drop = F],
                        psi = 1,
                        k = length(indsprop), tau = tau0, r = r0,
                        # method = "BFGS")
                        method = "Nelder-Mead", control = list(warn.1d.NelderMead = FALSE))


    b <- as.matrix(beta[indsprop], ncol = 1)
    if(any(b == 0)){ b[b == 0] <- 1e-7 }

    #lmlnew <- logmyb(b,
    lmlnew <- logmyb(as.matrix(opt$par, ncol = 1),
                     se[indsprop],
                     tau = tau0, psi = 1, r = r0, k = length(indsprop),
                     R = R[c(indsprop),c(indsprop),drop = F],
                     gval = opt$value)
    lpnew <- lmlnew + logmk(lprior, length(indsprop))

    betaprop <- rep(0, p)
    betaprop[indsprop] <- opt$par
    modelsizeprop <- length(indsprop)

    # lp
    # lpnew

    # log likelihood ratio
    # llr <- lpnew - lp
    # barker <- exp(lpnew)/(exp(lp) + exp(lpnew))
    barker <- 1/(1+exp(lp - lpnew))

      }

    if(is.infinite(lpnew)) { barker <- -1 }
    # if(opt$convergence %in% c(1, 10)) { barker <- -1 }

    u <- runif(1)

    if(u < barker){
      betavec <- betaprop
      modelsize <- modelsizeprop
      lml <- lmlnew
      lp <- lpnew
      inds <- indsprop
    }

    out[[1]][i,] <- betavec
    out[[2]][i] <- modelsize
    out[[3]][i] <- paste0(inds, collapse = ",")
    out[[4]][i] <- lml

    if(verbose){
      #if(i %% 100 == 0) { cat(sprintf("%i\n", i)) }
      setTxtProgressBar(pb, i)
    }

  }
  close(pb)
  cat(sprintf("\n%i iterations done in %.2f seconds\n", niter, (proc.time() - prc)[[3]]))
  }

  if(excl.burnin){
    out <- list(out[[1]][(burnin + 1):niter,],
                out[[2]][(burnin + 1):niter],
                out[[3]][(burnin + 1):niter],
                out[[4]][(burnin + 1):niter])
  }

  if(clump){
    out <- c(out,
             list(keeplist_cleaned))
  } else {
    out <- c(out, list(seq_len(p)))
  }

  return(out)

}
