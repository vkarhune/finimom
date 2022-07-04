# sampling from vectors of length 1
# see "?sample"
sample2 <- function(x, ...) x[sample.int(length(x), ...)]
