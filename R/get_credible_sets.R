#' Obtain credible sets from the posterior samples.
#'
#' @param samples Posterior samples from posterior_samples() function.
#' @param num_signals Number of signals.
#' @param level Probability level, default = 0.95.
#'
#' @return List of credible sets.
#' @export
#'
#' @examples
get_credible_sets <- function(samples, num_signals, level = 0.95){

  ctable <- sort(table(samples[[3]][samples[[2]] %in% num_signals]), decreasing = T)
  csind <- ctable |> prop.table() |> cumsum() |> (\(.) min(which(. > level)))()

  mods <- lapply(names(ctable)[1:csind], function(x) {
    list(unlist(strsplit(x, ",")), ctable[x])
  })



  cs_indices <- lapply(seq_len(num_signals), function(x){
    signal <- mods[[1]][[1]][-x]

    unlist(lapply(sapply(mods, "[", 1), function(mod){
      if(all(signal %in% mod)){
        setdiff(mod, signal)
      } else { NULL }
    }))

  })

  allmods <- lapply(split(expand.grid(cs_indices, stringsAsFactors = F),
                          seq_len(nrow(expand.grid(cs_indices)))), as.character)

  all_cs <- lapply(seq_len(num_signals), function(x) {
    unique(unlist(lapply(lapply(allmods, function(basemod){
      lapply(seq_len(num_signals), function(x){
        signal <- basemod[-x]

        # loop over allmods
        unlist(lapply(sapply(mods, "[", 1), function(mod){
          if(all(signal %in% mod)){
            setdiff(mod, signal)
          } else { NULL }
        }))
      })

  }), "[", x))) })



  credible_sets <- lapply(all_cs, function(x){
    sort(unlist(samples[[5]][as.numeric(x)]))
  })


}
