

#' This function do Clustering Analysis. Different distance methods are included, Euclidean, Viper Similairty and Correlation distance. Please choose the assay abd the slot to do the analysis.
#' @param  seu The seu expression object
#' res=0.6 for PA
#' @return The seu results
#' @author Junqiang Wang
#' @export
SeuPAClustering<-function(seu,
                            dims=50,
                          resolution=0.8,
                          features.sd.cutoff=0,
                           n.features=500,
                          assay=c("Activity","RNA"),
                          slot=c("data", "counts", "scaled.data"),
                          distance=c("viper_similarity","default",  "cor"),
                          nn=50,
                          vp.similarity.method="greater"
                          ){

  require(Seurat)
  require(ggplot2)

  if (is.null(seu@assays$Activity)==TRUE){

  message("Activity assay is NULL, please add activity assay")

 # seu<-AddActivityAssay(seu, vp_mat=as.matrix(seu@assays$RNA@counts))

  }

 # seu@active.assay <- 'Activity'

  DefaultAssay(seu)<-'Activity'

  message ("clustering using seurat default method...")

  seu@assays$Activity@scale.data<- seu@assays$Activity@data %>% as.matrix() %>% scale()

# Feature selection by sd

tmp.sd<-apply(as.matrix(seu@assays$Activity@data), 1, sd)

tmp.sd<-sort(tmp.sd, decreasing=TRUE)

features.1<-names(tmp.sd[1:n.features])

features.2<-tmp.sd[tmp.sd>features.sd.cutoff] %>% names()

features<-intersect(features.1, features.2)


#---------------------------
# Feature slection by MAD
#----------------------------

  # as.matrix(seu@assays$Activity@scale.data) -> Data
  #
  # mads<-apply(Data,1,mad) %>% sort( decreasing = TRUE)
  #
  # mads<-mads[1:min(n.HVF, length(mads))]
  #
  # fs<-names(mads)
  #
  # features<-fs


  seu <- RunPCA(seu, features = features)
  seu <- FindNeighbors(seu)
  seu <- FindClusters(seu, resolution = resolution)
  seu <- RunUMAP(seu, dims=1:dims)
  #seu <- RunTSNE(seu, dims=1:dims)

  seu@meta.data$Clusters_PA_EuclideanDist<-as.factor(as.numeric(seu@meta.data$seurat_clusters))


  if(distance=="viper_similarity"){


  message ("Clustering using viper similarity as distance...")

  vpmat<-as.matrix(seu@assays$Activity@data)

  vp_dist <- as.dist(viperSimilarity(vpmat, nn=nn, method=vp.similarity.method)) %>% as.matrix()

  seu@graphs <-FindNeighbors(vp_dist, distance.matrix=TRUE)

  seu <- FindClusters(seu, resolution =resolution, graph.name="snn")

  seu@meta.data$Clusters_PA_ViperDist<-as.factor(as.numeric(seu@meta.data$seurat_clusters))

  seu@meta.data$Clusters_PA_ViperDist<-paste0("C", seu@meta.data$Clusters_PA_ViperDist)

  }



  if(distance=="cor"){

    message ("clustering using viper similarity as distance...")

    vpmat<-as.matrix(seu@assays$Activity@data)

    vp_dist <- 1-cor(vpmat)

    seu@graphs <-FindNeighbors(vp_dist, distance.matrix=TRUE)

    #seu@graphs <-FindNeighbors(as.matrix(vp_dist), distance.matrix=TRUE, k.param = 20, prune.SNN = 1/15)

    #seu@graphs <-FindNeighbors(as.matrix(vp_dist), k.param = 50, prune.SNN = 1/15)

    seu <- FindClusters(seu, resolution =resolution, graph.name="snn")

    seu@meta.data$Clusters_PA_CorDist<-as.factor(as.numeric(seu@meta.data$seurat_clusters))

    seu@meta.data$Clusters_PA_CorDist<-paste0("C", seu@meta.data$Clusters_PA_CorDist)

  }


  seu@meta.data$seurat_clusters<- NULL

  return(seu)

}











SeuPAClusteringDefault<-function(seu, dims=50, resolution=0.6, assay=c("Activity","RNA"), slot=c("data", "counts", "scaled.data") ){

  require(Seurat)
  require(ggplot2)

  message ("clustering using seurat default method...")

  seu@active.assay <- 'Activity'

  seu@assays$Activity@scale.data<- as.matrix(seu@assays$Activity@data)

  seu<- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 500)


  #seu <- RunPCA(seu, features = rownames(seu))

  seu <- RunPCA(seu, features = VariableFeatures(seu))

  seu <- FindNeighbors(seu, dims=1:50)

  seu <- FindClusters(seu, resolution = resolution)


  seu <- RunUMAP(seu, dims=1:dims)
  seu <- RunTSNE(seu, dims=1:dims)


  seu@meta.data$PA.Clusters<-as.factor(as.numeric(seu@meta.data$seurat_clusters))

  return(seu)

}








