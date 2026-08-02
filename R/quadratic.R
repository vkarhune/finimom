quadratic <- function(b, Rsvd, q){
  sapply(seq_len(q), function(i){
  
  if(i == 1){
    lm1 <- as.matrix(1/Rsvd$values[1])
  } else {
    lm1 <- diag(1/Rsvd$values[1:i])
  }
  
  evec <- Rsvd$vectors[,1:i,drop = F]
  
  brb <- t(b) %*% evec %*% lm1 %*% t(evec) %*% b

  brb
  
})

}
