
#' This function convert count matrix to log2CPM matrix
#' @param count The count matrix
#' @return The log2CPM matrix
#' @author Junqiang Wang
#' @export
Log2CPM_OLD<-function(count){

  count <- apply(count,2,function(x){
    y <- 1e6*x/sum(x) + 1
    z <- log2(y)
    return(z)
  })

  count<-round(count, digits=2)

return(count)

}



#' This function convert count matrix to log2CPM matrix
#' @param count The count matrix
#' @return The log2CPM matrix
#' @author Junqiang Wang
#' @export
Log2CPM<-function(count){

  count <- t(t(count) / (colSums(count) / 1e6))
  count<-log2(count+1)

  count<-round(count, digits=2)

  return(count)

}







