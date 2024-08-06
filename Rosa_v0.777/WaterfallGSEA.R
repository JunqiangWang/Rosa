
#' @description  This function do the waterfall plot for a set of genes. The x-axis is the ranked signature and the y-axis
#' is the Z scores transformed by the ranked signature
#'
#' @export
WaterfallGSEA<-function(
    signature=signature,
    geneset=geneset,
    gsea.out=gsea.out,
    show.labels=TRUE,
    max.iter = 1e8,
    max.time = 333,
    label.cutoff.z=0.5
  ){

  require(ggplot2)
  require(ggrepel)

  signature<-signature %>% sort(decreasing = TRUE)
  rank.signature<-rank(signature)

  #rank.signature<-max(rank(signature))-rank(signature)+1

  ledge<-gsea.out$ledge %>% names()

  ledge<-intersect(ledge, names(signature[signature>label.cutoff.z]))

  df<-data.frame(gene=names(rank.signature), Stouffer.signature=signature, rank.signature=rank.signature)
  df$label<-rep(0, lenth=nrow(df))
  df$label[df$gene %in% ledge]<-1

  df$label2<-ifelse(df$label==1, df$gene, NA)

  df$label<-as.factor(df$label)

  df$ledge<-df$label

# ggplot visualization
# reference
# https://ggrepel.slowkow.com/articles/examples

p<- ggplot(df,
           aes(rank.signature,Stouffer.signature),
           fill=label2
)+
  scale_x_continuous(limits = c(min(df$rank.signature), max(df$rank.signature) * 1.3))+
  scale_y_continuous(position = "left") +
  geom_point(aes(color=ledge),
             size=3
             #shape = 21, color = "black", size = 3
             ) +
  scale_color_manual(values = c("0"="grey", "1"="orange"), name = "Ledge")+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank())+
  theme(legend.position = "bottom")+
     geom_text_repel(
     data = subset(df, label==1),
     aes(label = gene),
     force=0.5,
     nudge_x = max(df$rank.signature) * 0.2,
     direction="y",
     hjust=0.5,
     size = 3,
     segment.size=0.2,
     max.iter = max.iter,
     max.time =  max.time
  ) +
  labs(x="Rank Signature", y="Stouffer Signature")+
  scale_fill_discrete(labels=paste0(df$gene,': ',' Full names here'),
                      name='Significant genes') +
  theme(legend.position = 'right')

return(p)

}


