#'@description This function selects the activated regulators belongs to the reg types: TFs, coTFs or Sigs
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


  regulators<-top.n.markers$gene %>% unique()

  return(regulators)

}





