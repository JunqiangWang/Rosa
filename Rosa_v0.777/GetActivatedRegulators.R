
#'@description This function get the activated regulators belongs to the reg types: TFs, coTFs or Sigs
#' @param markers Differentially activated regulators
#' @param cRegulators A list of regulators
#' @param reg.type The regulatory types
#'@export
GetActivatedRegulators<-function(seu,
                         markers,
                         cRegulators,
                         reg.type="TFs",
                         top.n=10
){


  if(reg.type=="TFs"){
    tmp<- markers$gene %in% cRegulators$tfs
  }else if(reg.type=="coTFs"){
    tmp<-markers$gene %in% cRegulators$cotfs
  }else if(reg.type=="Sigs"){
    tmp<-markers$gene %in% cRegulators$sigs
  }

  markers<-markers[which(tmp==TRUE), ]


  top.n.markers<-markers %>%
    arrange(p_val, desc(avg_diff)) %>%
    distinct(gene, .keep_all = TRUE) %>%
    group_by(cluster) %>%
    slice_max(n = top.n, order_by = -p_val, with_ties = FALSE)

  #
  # top.n.markers$cluster %>% unique() ->cluster.unique
  # lapply(cluster.unique, FUN=function(x, cluster, top.n){
  #   tmp<-grep(x, cluster)
  #   tmp<-tmp[1:min(length(tmp), top.n)]
  # }, cluster=top.n.markers$cluster, top.n=top.n) %>%
  #   unlist() -> tmp
  #
  # top.n.markers<-top.n.markers[tmp,]
  #
  # top.n.markers<-top.n.markers %>%
  #   arrange(cluster)

  # Chose the top n markers

  regulators<-top.n.markers$gene %>% unique()

  return(regulators)

}





