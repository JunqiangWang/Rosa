#' This function performs Stouffer integration of the PA in each cluster by rows.
#' @param  seu The seurat object
#' @group.by The colname of the metadata used
#' @export
GenerateStoufferSignatureFromSeu<-function(seu,
                                          assay="Activity",
                                          slot="data",
                                          group.by=NULL
                                         ){

  require(Seurat)
  require(dplyr)

  Idents(seu) <- group.by

  clusters<-seu@meta.data %>%
    dplyr::select(!!sym(group.by)) %>%
    pull() %>%
    unique() %>%
    as.character()

  clusters

   mat<-GetAssayData(object = seu[[assay]], slot = slot) %>% as.matrix()

   nes.stouffer<-matrix(0, nrow=nrow(mat), ncol=length(clusters))
   rownames(nes.stouffer)<-rownames(mat)
   colnames(nes.stouffer)<-clusters


for (i in 1:length(clusters)){

   cluster.i<-clusters[i]

   seu.i <- subset(seu, idents =cluster.i)

   mat<-GetAssayData(object = seu.i[[assay]], slot = slot) %>% as.matrix()

   nes.stouffer.i<-rowSums(mat)/sqrt(ncol(mat))

   nes.stouffer[,i]<-nes.stouffer.i

}

 return(nes.stouffer)

}



