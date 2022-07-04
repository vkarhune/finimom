generate_sumstats <- function(phenotype, X, adj = NULL){

  out <- t(vapply(seq_len(ncol(X)), function(i){

    xx <- X[,i]

    if(!(is.null(adj))) { xx <- cbind(xx, adj)}

    modelfit <- lm.fit(x = cbind(1, xx), y = phenotype)

    B <- coef(modelfit)[2]
    SE <- tryCatch( sqrt(diag((sum(resid(modelfit)^2)/modelfit$df.residual)*chol2inv(modelfit$qr$qr)))[2],
                    error=function(error){return(NA)})
    logP <- tryCatch({log(2) + pt(abs(B/SE), df = modelfit$df.residual, lower.tail = FALSE, log.p = TRUE)},
                     error=function(error){return(NA)})
    P <- ifelse(logP < log(.Machine$double.xmin), 1e-300, exp(logP))

    return(c(B, SE, logP, P))

  }, FUN.VALUE=numeric(4), USE.NAMES = FALSE))

}
