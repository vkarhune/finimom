select_causals <- function(numcausals, LDmat, mincorr = 0.5, maxcorr = 0.8, minmutualcorr = 0.5){

# numcausals <- sample(c(1, 2, 4, 8), size = 1)
causals <- numeric(length = numcausals)

cond <- FALSE

while(!cond){
  causalindex <- sample(seq_len(nrow(LDmat)), size = 1)
  if(sum(abs(LDmat[causalindex,]) > mincorr & abs(LDmat[causalindex,]) < maxcorr) > 1){
    cond <- TRUE
  }
}

# causalindex

if(numcausals == 1) causals <- causalindex

# select a second causal
if(numcausals > 1){
  chooseinds <- setdiff(which(abs(LDmat[causalindex,]) > mincorr &
                                abs(LDmat[causalindex,]) < maxcorr), causalindex)
  causal2 <- sample(chooseinds, size = 1)
  causals[1:2] <- c(causalindex, causal2)
}

# select more causals
if(numcausals > 2){
  cond <- FALSE

  while(!cond){
    chooseinds2 <- which(apply(LDmat, 1, function(x) sum(abs(x) > minmutualcorr)) > 1)

    causal3 <- sample(chooseinds2, size = numcausals - 2)
    # causal3 <- sample(nrow(dd), size = numcausals - 2) # alternative
    tm <- LDmat[c(causalindex, causal2, causal3),c(causalindex, causal2, causal3)]
    if(max(abs(tm[lower.tri(tm)])) < maxcorr){
      cond <- TRUE
    }
  }
  causals[3:numcausals] <- causal3
}

sort(causals)
}
