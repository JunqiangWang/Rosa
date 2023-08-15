
#' RSEA - regulon set enrichment analysis algorithm
#' Improvement: no rank transformation 
#' Improvement: use the MAD feature selection measure as weight for MetaRsea
#' 
#' @author  Junqiang Wang
#' @export

Rsea <- function (eset, regulon, wt.features, minsize=20) {
  
  if(length(names(regulon))<100){
    regulon<-regulon[1][[1]]}
  
  
  if(is.null(wt.features)==TRUE){
    
    message("Null wt.features is depreciated!")
    
    wt.features<-rep(1, times=nrow(eset))
    
    names(wt.features)<-rownames(eset)
  }
  
  # Data Preprocessing # Step 1 - Filter the expression data.
  
  # 1.1 Get a list of all targets of every regulator in the network.
  
  tmp <- c(names(regulon), unlist(lapply(regulon, function(x) names(x$tfmode)), use.names = FALSE))
  
  # 1.2 Remove expression data of all the genes which aren't in the list identified above. 
  
  eset <- eset[rownames(eset) %in% unique(tmp),]
  
  # 1.3 Scale the expression matrix on rows [i.e. row = (row - mean(row)) / sd(row) ]. 
  
  # Scale the expression matrix on cols

 # tt<-scale(eset)
  
  tt <- t(scale(t(eset)))
  
  
  # Step 2 - Filter the regulon object (i.e. the regulatory network).
  
  # 2.1 Remove targets in the regulon which are not in the rownames of the expression matrix. 
  
  regulon <- lapply(regulon, function(x, genes) {
  
  filtro <- names(x$tfmode) %in% genes 
  
  x$tfmode <- x$tfmode[filtro] 
  
  if (length(x$likelihood) == length(filtro)) 
    x$likelihood <- x$likelihood[filtro] 
  return(x) }, genes = rownames(eset))

# 2.2 Define minimum regulon size for filtering (default is 20).

# The 'minsize' parameter is specified in the function parameters.

# 2.3 Remove regulators with a regulon size below the 'minsize' parameter. 
  
 # regulon <- regulon[sapply(regulon, function(x) length(x$tfmode)) >= minsize]

# rsea
  
  nes <- Rsea.core(tt, regulon, wt.features=wt.features)

# Return VIPER matrix 
  
  return(nes)

}



# regulon set enrichment analysis
Rsea.core <- function (tt, regulon, wt.features) {
  
  if(length(names(regulon))<100){
    regulon<-regulon[1][[1]]}
    
    # Step 1 - Create the 'Mode of Regulation' and 'Weights' matrices.
    # 1.1 Get a list of all targets of every regulator in the network.
    
    targets <- unique(unlist(lapply(regulon, function(x) names(x$tfmode)), use.names = FALSE))
    
    # 1.2 Create the Mode of Regulation matrix from the regulon object. 
    mor <- sapply(regulon, function(x, genes) {
    
    return(x$tfmode[match(genes, names(x$tfmode))]) }, genes = targets)

# 1.2 Create the Weights matrix from the regulon object. 
    
    wts <- sapply(regulon, function(x, genes) {

tmp <- x$likelihood[match(genes, names(x$tfmode))] 

tmp[is.na(match(genes, names(x$tfmode)))] <- NA 

return(tmp/max(tmp, na.rm = T)) }, genes = targets)

# 1.3 For each regulator, assign values of 0 to genes which are not listed as its targets. 
    
    mor[is.na(mor)] <- 0 
    
    wts[is.na(wts)] <- 0 

# 1.4 Scale the columns of the 'Weights' matrix to the sum of the weights. 
    
    # wtss <- scale(wts, center = FALSE, scale = colSums(wts))

# Step 2 - Calculate the two-tail enrichment scores.

# 2.1 Calculate the 'T2 rank' matrix from the expression dataset. 
    
  #  t2 <- apply(tt, 2, rank)/(nrow(tt) + 1)
    

# This line of code is necessary to match the order of genes # for the subsequent matrix multiplication steps.

pos <- match(targets, rownames(tt))

# 2.2  t2q are the normalized values. 

t2q<-filterRowMatrix(tt, pos)

tmp<-match(targets, names(wt.features))

wt.features<-wt.features[tmp]

wt.features.mat<-matrix(0, nrow=nrow(mor), ncol=ncol(mor))

wt.features[is.na(wt.features)] <- 0 

for(i in 1:ncol(wt.features.mat)){
  
  wt.features.mat[,i]<-wt.features
}
  
# 2.3 Matrix multiplication. 

#total weights 
wtss<-mor*wts*wt.features.mat
#wtss<-mor*wts
#wtss<-apply(wtss, 2, function(x){
#  x<-x*wt.features
#})


message("the integrated weigted is calculated!")

wtss<-t(wtss)

# the Enrichment score
es <- wtss %*% t2q

message("the ES is calculated!")

# The sd
wtss2<-wtss^2

es.var <- apply(wtss2, 1, sum)

es.sd<-sqrt(es.var)


# get the NES
nes<-apply(es, 2, function(x, es.sd){
  
  x.nes<-x/es.sd
  return(x.nes)
  
}, es.sd=es.sd)

message("the NES is calculated!")

return(nes)
}




# metaRSEA

MetaRsea<-function(eset, regulon.list, wt.features=NULL, wt.amplification=1,  minsize=20){
  
  require(parallel)
  
  pa.list<-mclapply(regulon.list, function(x){
    
    x<-Rsea(eset=eset, regulon = x, wt.features=wt.features)
    
  })
  
  message("PAs for each regulome are inferred!")
  
  
# generate the augmented matrices
  
  regulators<-lapply(pa.list, FUN=rownames) %>% unlist() %>% unique() %>% sort()
  
  cells<-colnames(eset)


pa.aug.list<-vector(mode="list", length=length(pa.list))

names(pa.aug.list)<- names(pa.list)

# assign the values by rows


for (i in 1:length(pa.aug.list)){
  
  x<-matrix(0, nrow=length(regulators), ncol=length(cells))
  rownames(x)<-regulators
  colnames(x)<-cells

  tmp<-match(rownames(pa.list[[i]]), rownames(x))
  
  x[tmp,]<-pa.list[[i]]
  
  pa.aug.list[[i]]<-x
  
}


 x<-matrix(0, nrow=length(regulators), ncol=length(cells))
 rownames(x)<-regulators
 colnames(x)<-cells
 
 
# 
 MAD.list<-lapply(pa.aug.list, function(x){
   
   # mean absolute difference as weight
   row.mean<-apply(x, 1, mean)
   x2<-apply(x, 2, FUN=function(x){
     x<-x-row.mean
     x<-abs(x)
     return(x)}
   )
   
   x2<-rowSums(x2)/ncol(x2)
   
   x<-x2
   
 })
 
 
for(i in 1:nrow(x)){
  
  # the wt of each proteins
  wt.i<-lapply(MAD.list, function(x){
    
    x<-x[i]
   
  }) %>% unlist() %>% unname()
  
  wt.i<-wt.i*wt.amplification

  for(j in 1:ncol(x)){

   nes.ij<- lapply(pa.aug.list, function(x){
      y<-x[i,j]
    }) %>% unlist() %>% unname()
  

   # weight by sd
  # nes.ij<-sum(nes.ij*wt.i)/sqrt(sum(wt.i^2))
   
   # weight
   nes.ij<-sum(nes.ij*abs(nes.ij)*wt.i)/sqrt(sum((abs(nes.ij)*wt.i)^2))

   x[i,j]<-nes.ij
  }
}
 
 pa.nes<-x

  message("The NESs are calculated!")
  
  return(pa.nes)
  
}



#' @paramters exprs count matrix
#' @export
WeightFeature<-function(seu){
  #
  
  DefaultAssay(seu)<-"RNA"
  
  seu <- NormalizeData(seu)
  #
  seu <- FindVariableFeatures(seu, 
                              assay="RNA",
                              selection.method = "vst", 
                              nfeatures = dim(seu)[1])
  #
  features.ranked<-VariableFeatures(seu)
  
  features.wt<-1-seq(1:length(features.ranked))/length(features.ranked) 
  
  names(features.wt)<-features.ranked
  
  return(features.wt)
  
  return(features)
}






