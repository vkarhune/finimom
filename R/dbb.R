# beta-binomial prior for model size:
dbb <- function(x, p, a, b){
  lchoose(p, x) + lbeta(x + a, p - x + b) - lbeta(a, b)
}
