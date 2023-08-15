
#' @description  This function do the waterfall plot for a set of genes. The x-axis is the ranked signature and the y-axis 
#' is the Z scores transformed by the ranked signature
#'
#'
#'
#'
#' @export
WaterfallGSEA<-function(
    signature=signature, 
    geneset=geneset,
    show.labels=TRUE,
    label.cutoff.z=1.65
  ){
  
  require(ggplot2)
  require(ggrepel)
  
  gsea<-gsea(signature=signature, geneset=geneset, per=100)
  
  ledge<-gsea$ledge %>% names()
  
  signature<-signature %>% sort(decreasing = TRUE)
  
  rank.signature<-rank(signature)
  
  df<-data.frame(gene=names(rank.signature), Stouffer.signature=signature, rank.signature=rank.signature)
  df$label<-rep(0, lenth=nrow(df))
  df$label[df$gene %in% ledge]<-1
  
  df$label2<-ifelse(df$label==1, df$gene, NA)
  
  df$label<-as.factor(df$label)
  
  df$ledge<-df$label
  
#ggplot visualization
  
 p<- ggplot(df,
            aes(rank.signature,Stouffer.signature),
           fill=label2
            )+
   geom_point(aes(color=ledge)) +
   scale_color_manual(values = c("grey", "red"))+
 theme_bw(base_size = 12) + 
   theme(legend.position = "bottom") +
   geom_text_repel(
     data = subset(df, label==1),
     aes(label = gene),
     size = 3,
     box.padding = unit(0.35, "lines"),
     point.padding = unit(0.3, "lines")
   ) +
   scale_fill_discrete(labels=paste0(df$gene,': ',' Full names here'),
                       name='Significant genes') +
   theme(legend.position = 'right')
 
return(p)
   
}



