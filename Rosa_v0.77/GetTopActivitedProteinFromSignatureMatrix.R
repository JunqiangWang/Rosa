#' @description  This function selects top activated proteins from signatures
#' @param  seu The seurat object
#' @param group.by The column name in the seu@meta.data
#' @return Signatures for each group
#' @export
GetTopActivitedProteinFromSignatureMatrix<-function(mat,
                                                    top.n=10){

  out<-data.frame()

  for (i in 1:ncol(mat)){
   signature.i<-mat[,i]
   signature.i<-sort(signature.i, decreasing=TRUE)
   signature.i<-signature.i[1:top.n]
   signature.i<-data.frame(protein=names(signature.i), signature=unname(signature.i), group=rep(colnames(mat)[i], top.n))

   out<-rbind(out, signature.i)
  }

  return(out)

}
