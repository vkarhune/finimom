get_multitrait_credible_sets <- function(input, num_signals, level, purity, R){

  p <- nrow(R)

  ctable <- sort(table(input[[3]][apply(input[[2]], 2, function(x) all(x == num_signals))]),
                 #ctable <- sort(table(input1[[3]][apply(input1[[2]], 2, function(x) all(x == num_signals))]),
                 #ctable <- sort(table(input2[[3]][apply(input2[[2]], 2, function(x) all(x == num_signals))]),
                 decreasing = T)
  csind <- (function(.) min(which(. > level)))(cumsum(prop.table(ctable)))
  mods <- lapply(names(ctable)[1:csind], function(x) {
    list(unlist(strsplit(x, ",")), ctable[x])
  })

  sets1 <- lapply(seq_along(num_signals), function(i){
    # sets2 <- lapply(seq_along(num_signals), function(i){
    # numsig <- num_signals[1]
    # numsig <- num_signals[2]
    # i <- 1

    numsig <- num_signals[i]
    if(i == 1){
      modss <- lapply(mods, function(x) as.numeric(x[[1]][seq_len(numsig)]))
    } else {
      # modss <- lapply(mods, function(x) as.numeric(x[[1]][seq_len(numsig) + num_signals[1:(i - 1)]]))
      modss <- lapply(mods, function(x) as.numeric(x[[1]][seq_len(numsig) + sum(num_signals[1:(i - 1)])]))
    }

    cs_indices <- lapply(seq_len(numsig), function(x) {
      #signal <- mods[[1]][[1]][-x]
      signal <- modss[[1]][-x]
      #unlist(lapply(sapply(mods, "[", 1), function(mod) {
      unlist(lapply(modss, function(mod) {
        if (all(signal %in% mod)) {
          setdiff(mod, signal)
        }
        else {
          NULL
        }
      }))
    })

    cs_indices <- lapply(cs_indices, unique)

    allmods <- lapply(split(expand.grid(cs_indices, stringsAsFactors = F),
                            seq_len(nrow(expand.grid(cs_indices)))), as.character)

    all_cs <- lapply(seq_len(numsig), function(x) {
      unique(unlist(lapply(lapply(allmods, function(basemod) {
        lapply(seq_len(numsig), function(z) {
          signal <- basemod[-z]
          #unlist(lapply(sapply(mods, "[", 1), function(mod) {
          unlist(lapply(modss, function(mod) {
            if (all(signal %in% mod)) {
              setdiff(mod, signal)
            }
            else {
              NULL
            }
          }))
        })
      }), "[", x)))
    })



    if (!(is.null(purity))) {
      credible_sets <- lapply(all_cs, function(cs) {
        #indices <- unlist(samples[[5]][as.numeric(cs)])
        indices <- as.numeric(cs)
        if(i > 1){ indices <- indices - p*(i - 1) }
        Rcs <- R[indices, indices, drop = F]
        Rcs[upper.tri(Rcs)] <- 1
        arrinds <- unique(which(abs(Rcs) < purity, arr.ind = T)[,
                                                                1])
        if (length(arrinds) == 0) {
          sort(indices)
        }
        else {
          sort(indices[-arrinds])
        }
      })
    }

    credible_sets

  })

}
