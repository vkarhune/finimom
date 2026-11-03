truncbb <- function(x, p, a, maxsize, k){
  lprobs <- sapply(seq_len(maxsize), function(y) lchoose(p, y) + lbeta(y + a, p - y + x) - lbeta(a, x) )
  probs <- exp(lprobs - max(lprobs)) / sum(exp(lprobs - max(lprobs)))
  out <- sum(probs*seq_len(maxsize)) - k
  return(out)
}
