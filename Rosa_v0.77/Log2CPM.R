#' This function converts count matrix to log2CPM matrix
#' @param count The count matrix
#' @return The log2CPM matrix
#' @export
Log2CPM<-function(count){

  count <- t(t(count) / (colSums(count) / 1e6))
  count<-log2(count+1)

  count<-round(count, digits=2)

  return(count)

}







