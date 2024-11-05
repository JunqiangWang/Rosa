#' Generate signature from the specified assay and slot
#' @description  This script generates protein activity signatures by mean
#' @param  seu The Seurat object
#' @param group.by The column name in the meta.data of the Seurat object
#' @param assay assay used
#' @param slot slot used
#' @return Signatures for each group
#' @export
GenerateSignatureByMean<-function(seu,
                                  assay="SCT",
                                  slot="data",
                                  group.by=NULL){


  tmp<-match(group.by, colnames(seu@meta.data))
  groups<-as.vector(sort(unique(seu@meta.data[,tmp])))

  dat <- GetAssayData(object = seu, assay = assay, slot = slot)
  out<-matrix(0, nrow(dat), length(groups))

  rownames(out)<-rownames(dat)
  colnames(out)<-groups

  
  for (i in 1:length(groups)){
    reg_pattern <- FetchData(object = seu, vars = group.by)
    seu.i<-seu[, which(x=reg_pattern==groups[i])]
    
    dat<-as.matrix(GetAssayData(object = seu.i, assay = assay, slot = slot))
    dat<-rowSums(dat)/ncol(dat)

    out[,i]<-dat
  }

  out<-round(out, digits=2)
  return(out)
}




