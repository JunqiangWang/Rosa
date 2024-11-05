#' This function generates an heatmap plot of the differential markers
#' @param seu seurat object
#' @param cluster.by the ID in seu@meata.data
#' @param reverse.clusters whether the order of the cluster id shoud be reversed
#' @param cluster.order the order of the cluster
#' @param groups.use the groups that are used run the presto
#' @param top.n return the top.n dif features
#' @param assay The assay that is used
#' @param slot the slot of the assay that is used
#' @return return the markers and store the markers in a new assay
#' @export
GetSeuMarkersActivity<-function(seu,
                        cluster.by="seurat_clusters",
                        reverse.clusters=FALSE,
                        ident.1=NULL,
                        ident.2=NULL,
                        top.n=25,
                        slot.use="scale.data",
                        only.pos = TRUE,
                       test.use="t",
                      clusters.ordered = NULL,
                      min.pct = 0,
                      myAUC.cutoff=0,
                      logfc.threshold =-Inf,
                      avg_diff.threshold=0,
                      p_val_adj.threshold=1){

  require(dplyr)

  require(Seurat)

  DefaultAssay(object = seu) <- "Activity"

   Idents(seu)<-cluster.by

seu@meta.data$Clusters<-seu@meta.data[,match(cluster.by, colnames(seu@meta.data))]

seu@meta.data$NES.Stouffer<-NA

tmp.mat<-GetAssayData(object = seu[["Activity"]], slot = slot.use)

tmp.cluster.vector<-unname(unlist(subset(seu@meta.data, select=cluster.by)))

clusters<-sort(unique(tmp.cluster.vector), decreasing = FALSE)

if(is.null(clusters.ordered)==FALSE){clusters<-clusters.ordered}


if(is.null(ident.1)==TRUE){

markers <-FindAllMarkers(seu,
                         test.use = "t",
                         only.pos = TRUE,
                         slot="scale.data",
                         assay ="Activity",
                         min.pct = min.pct,
                         logfc.threshold = logfc.threshold)


} else if(is.null(ident.1)==FALSE){

 markers <- FindAllMarkers(seu, ident.1=ident.1, ident.2=ident.2, slot="data", assay ="Activity", test.use = test.use, only.pos = TRUE, min.pct = min.pct.threshold)

  clusters<-c(ident.1, ident.2)

}


print(markers)


if(reverse.clusters==FALSE){

  clusters<-clusters

} else if (reverse.clusters==TRUE){

  clusters<-rev(clusters)

}


if(test.use =="t"){

  top.n.markers<-markers %>%
    arrange(p_val, desc(avg_diff)) %>%
    group_by(cluster) %>%
    slice_max(n = top.n, order_by = desc(p_val)) %>%
    distinct(gene, .keep_all = TRUE)


#
  top.n.markers$cluster %>% unique() ->cluster.unique
lapply(cluster.unique, FUN=function(x, cluster, top.n){
tmp<-grep(x, cluster)
tmp<-tmp[1:min(length(tmp), top.n)]
}, cluster=top.n.markers$cluster, top.n=top.n) %>%
  unlist() -> tmp

top.n.markers<-top.n.markers[tmp,]

print(top.n.markers, n=100)



}else if (test.use =="roc") {

  markers %>%
    group_by(cluster) %>%
      filter(myAUC>myAUC.cutoff) %>%
    filter(avg_diff>avg_diff.threshold) %>%
    top_n(n = top.n, wt=avg_diff) -> top.n.markers

}



all.cluster.markers<-vector()




for (cluster.i in clusters){

  cluster.markers<-top.n.markers%>%
    dplyr::filter(cluster==cluster.i)

  markers.id<-match(cluster.markers$gene, rownames(tmp.mat))

  cells<-as.data.frame(seu@meta.data) %>%
    dplyr::filter(Clusters==cluster.i)

  cells.id<-match(rownames(cells), colnames(tmp.mat))

  mat<-tmp.mat[markers.id, cells.id ] %>% as.matrix()

  NES.Stouffer<-colSums(mat)/sqrt(ncol(mat))

  seu@meta.data$NES.Stouffer[match(names(NES.Stouffer), rownames(seu@meta.data))]<-NES.Stouffer

  all.cluster.markers<-c(all.cluster.markers, cluster.markers$gene)

}




all.markers.id<-match(all.cluster.markers, rownames(tmp.mat))

 mat<-tmp.mat[all.markers.id, ] %>% as.matrix()


tmp<-seu@meta.data%>%
  as.data.frame() %>%
  arrange(Clusters, -NES.Stouffer) %>%
  rownames()


mat<-mat[, match(tmp, colnames(mat))]


seu.markers<-CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)

seu.markers@meta.data<-seu@meta.data[match(colnames(mat), rownames(seu@meta.data)), ]

seu.markers@assays$RNA@misc$diffMarkers<-all.cluster.markers

seu.markers<-AddMarkerAssay(seu=seu.markers, marker.mat=mat)

seu.markers@assays$Marker@misc$diffMarkers<-all.cluster.markers


return(seu.markers)

}






