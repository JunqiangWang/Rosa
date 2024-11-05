#' This function generates a Heatmap of the differential markers
#' @param seu seurat object
#' @param cluster.by the column ID in the seu@meata.data
#' @param reverse.clusters whether the orders of the cluster id should be reversed
#' @param cluster.order the order of the cluster
#' @param groups.use the groups that are used 
#' @param top.n return the top n differential features
#' @param assay The assay that is used
#' @param slot the slot of the assay that is used
#' @export
GetActivityMarkers<-function(seu,
                                cluster.by="seurat_clusters",
                                reverse.clusters=FALSE,
                                ident.1=NULL,
                                ident.2=NULL,
                                top.n=NULL,
                                slot.use="scale.data",
                                only.pos = TRUE,
                                test.use=c("t", "wilcox"),
                                clusters.ordered = NULL,
                                min.pct.threshold = 0,
                                myAUC.cutoff=0,
                                logfc.threshold =-Inf,
                                avg_diff.threshold=-Inf,
                                p_val_adj.threshold=1){

  require(dplyr)

  require(Seurat)

if(length(test.use)==2){ test.use="t"}

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
                             test.use="t",
                             only.pos = TRUE,
                             slot="scale.data",
                             assay ="Activity",
                             min.pct =0,
                             logfc.threshold = -Inf)


  } else if(is.null(ident.1)==FALSE){

    markers <- FindAllMarkers(seu,
                              ident.1=ident.1,
                              ident.2=ident.2,
                              slot="data",
                              assay ="Activity",
                              test.use = test.use,
                              only.pos = TRUE,
                              min.pct = min.pct.threshold)

    clusters<-c(ident.1, ident.2)

  }


  print(markers)


  if(is.null(top.n)){top.n<-nrow(markers)}


  if(reverse.clusters==FALSE){

    clusters<-clusters

  } else if (reverse.clusters==TRUE){

    clusters<-rev(clusters)

  }


  if(test.use =="t"|test.use =="wilcox"){

top.n.markers<- markers %>%
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
      filter(myAUC>myAUC.cutoff) %>%
      filter(avg_diff>avg_diff.threshold) %>%
      arrange(desc(myAUC), desc(avg_diff)) %>%
      group_by(cluster) %>%
      slice_max(n = top.n, order_by = myAUC) %>%
      distinct(gene, .keep_all = TRUE)-> top.n.markers

    #

    top.n.markers$cluster %>% unique() ->cluster.unique
    lapply(cluster.unique, FUN=function(x, cluster, top.n){
      tmp<-grep(x, cluster)
      tmp<-tmp[1:min(length(tmp), top.n)]
    }, cluster=top.n.markers$cluster, top.n=top.n) %>%
      unlist() -> tmp

    top.n.markers<-top.n.markers[tmp,]

    print(top.n.markers, n=100)

  }

  return(top.n.markers)

}












