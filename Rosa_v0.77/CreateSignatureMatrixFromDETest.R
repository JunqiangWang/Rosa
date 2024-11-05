#' @param X A DEG data frame
#' 
CreateSignatureMatrixFromDETest<-function(X){
  
  require(stringr)
  require(dplyr)
  
  
  split_dfs_0 <- X %>%
    group_by(cluster) %>%
    group_split()
  
  
  split_dfs <- lapply(split_dfs_0, function(x) {
    x<-x %>% 
      select(gene, mean_diff)
  })
  
  
  split_dfs_names <- lapply(split_dfs_0, function(x) {
    cluster_id <- unique(x$cluster) %>% as.character()
  })
  
  names(split_dfs)<-split_dfs_names
  
  merged_df <- reduce(split_dfs, left_join, by = "gene")
  colnames(merged_df)<-c("gene", split_dfs_names)
  
  return(merged_df)
}


