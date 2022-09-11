#' Title
#'
#' @param R
#' @param clump_r2
#' @param z
#'
#' @return
#' @export
#'
#' @examples
clump_variants <- function(R, clump_r2, z){


  ldlist <- lapply(seq_len(nrow(R)), function(i){
    c(i, setdiff(which(abs(R[i,]) > sqrt(clump_r2)), i))
  })

  ldlist_sort <- ldlist[order(-abs(z))]




  #
  keeplist_sort <- vector("list", nrow(R))

  keeplist_sort[[1]] <- ldlist_sort[[1]]

  for(i in 2:length(ldlist_sort)){
    if(ldlist_sort[[i]][[1]] %in% unlist(keeplist_sort[1:(i-1)])){
      NULL } else {
        keeplist_sort[[i]] <- setdiff(ldlist_sort[[i]], unlist(keeplist_sort[1:(i-1)]))
    }
  }



  keeplist_cleaned <- keeplist_sort[sapply(keeplist_sort,
                                           function(x) length(x) > 0)]

}
