
#'@description This function generates dot plot of the regulators
#'@param markers The markers to plot
#'@param idents The idents 
#'@param title The title of the plot. 
#'@param top.n The top n markers to plot
#'@param limits The vector of upper and lower limits for the data, for example: limits=c(-2, 2)
#'@export
DotPlotActivity<-function(seu,
                         markers,
                         idents=NULL,
                         title=c("TFs", "coTFs", "Sigs"),
                         limits=NULL
){

  require(tidyverse)
  require(ggdendro)
  require(cowplot)
  require(ggtree)
  require(patchwork)
  require(scales)


if(is.null(idents)==FALSE){
Idents(seu)<-idents
}

  if(is.null(limits)){
    p<-DotPlotActivitySeu(seu,
                       features = markers,
                       # assay="Activity",
                       dot.scale=6,
                       scale=FALSE) +
      scale_colour_gradient2(
        low = "blue",
        mid = "white",
        high = "red"
      )+

      FontSize(7) +
      theme(axis.text.x = element_text(size = 6, angle=45, hjust=1))+
      theme(axis.text.y = element_text(size = 6, angle=0, hjust=1))+
      theme(legend.title  = element_text(size = 10))+
      theme(legend.text  = element_text(size = 10))+
      theme(axis.title.x = element_text(size = 0, face="plain")) +
      theme(axis.title.y = element_text(size = 0, face="plain")) +
      coord_flip() +
      ggtitle(title)

  } else{


    mat.pa<-seu@assays$Activity@data %>% as.matrix()
    mat.pa[mat.pa>limits[2]]<-limits[2]
    mat.pa[mat.pa<limits[1]]<-limits[1]
    seu@assays$Activity@data<-mat.pa
    seu@assays$Activity@scale.data<-mat.pa

    p<-DotPlotActivitySeu(seu,
                       features = markers,
                       # assay="Activity",
                       dot.scale=6,
                       scale=FALSE) +
      scale_colour_gradient2(
        limits=limits,
        low = "blue",
        mid = "white",
        high = "red"
      )+

      FontSize(7) +
      theme(axis.text.x = element_text(size = 6, angle=45, hjust=1))+
      theme(axis.text.y = element_text(size = 6, angle=0, hjust=1))+
      theme(legend.title  = element_text(size = 10))+
      theme(legend.text  = element_text(size = 10))+
      theme(axis.title.x = element_text(size = 0, face="plain")) +
      theme(axis.title.y = element_text(size = 0, face="plain")) +
      coord_flip() +
      ggtitle(title)

  }

  return(p)

}

