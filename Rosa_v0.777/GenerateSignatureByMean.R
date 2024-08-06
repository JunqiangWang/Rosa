#' Generate signature from the specified assay and slot
#'
#' @description  This script generate protein activity signature by mean
#'
#' @param  seu The Seurat object
#' @param group.by the column name of the meta.data in the Seurat object
#' @param assay assay used
#' @param slot slot used
#'
#' @return Signatures for each group
#' @examples
#' \dontrun{
#' # Assuming you have a Seurat object named seu
#' signatures <- GenerateSignatureByMean(seu, assay = "RNA", slot = "data", group.by = "cell_type")
#' head(signatures)
#' }
#' @export
GenerateSignatureByMean<-function(seu,
                                  assay="SCT",
                                  slot="data",
                                  group.by=NULL){


  tmp<-match(group.by, colnames(seu@meta.data))
  groups<-as.vector(sort(unique(seu@meta.data[,tmp])))

  dat <- GetAssayData(object = seu, assay = assay, slot = slot)
  out<-matrix(0, nrow(dat), length(groups))

  # Extract data from the specified assay and slot

  rownames(out)<-rownames(dat)
  colnames(out)<-groups

  #
  for (i in 1:length(groups)){
    reg_pattern <- FetchData(object = seu, vars = group.by)
    seu.i<-seu[, which(x=reg_pattern==groups[i])]
    #
    dat<-as.matrix(GetAssayData(object = seu.i, assay = assay, slot = slot))
    dat<-rowSums(dat)/ncol(dat)
    # convert to RPKM
    #exp<-log2(1e6*exp/sum(exp)+1)
    out[,i]<-dat
  }

  out<-round(out, digits=2)
  return(out)
}




