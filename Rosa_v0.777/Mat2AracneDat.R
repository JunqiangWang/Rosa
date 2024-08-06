

#' This function generate the arance mat files from matrix
#' @param  log2CPM if on, use log2CPM transformation
#' @param  write.output if TRUE, output will write to files
#' @export
Mat2AracneDat<-function( mat,
                      dat.name=NULL,
                      log2CPM=FALSE,
                      write.output=TRUE,
                      subsampling=FALSE,
                      r.cutoff=0,
                      subsampling.size=NULL){

  if(is.null(dat.name)){
  stop('please name the dat.name')
  }


dir_main<-getwd()

dir_sub <- paste0(dir_main, "/in_mat_aracne")


if(file.exists(dir_sub)==FALSE){

 dir.create(file.path(dir_sub))
}


file.out<-paste0(dir_sub, "/", dat.name, ".dat")

require(dplyr)

  if(log2CPM==TRUE){

  mat<-Log2CPM(as.matrix(mat))
  }


  if (subsampling==TRUE){

id<-sample(1:ncol(mat), size=min(subsampling.size, ncol(mat)))

mat<-mat[,id]
}

exprs<-mat

exprs<-round(exprs, 2)

# remove genes do not express in 95% samples

nonzeors<-apply(exprs, 1, FUN=function(x){
  length(which(x != 0))
  }
  )

r<-nonzeors/ncol(exprs)

exprs<-exprs[which(r>r.cutoff), ] %>% as.matrix()

exprs<-round(exprs, digits=2)

colnames.exprs<-c("gene", colnames(exprs))
exprs<-cbind(rownames(exprs), exprs)
exprs<-rbind(colnames.exprs, exprs)
exprs<-as.data.frame(exprs)

rownames(exprs)<-NULL
colnames(exprs)<-NULL

message("dim of the output exprs is: ", dim(exprs)[1], " ", dim(exprs)[2])


if(write.output==TRUE){

write.table(exprs, file=file.out, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)
}

message("writing is finished!")
}


