
#' RSEA - regulon set enrichment analysis algorithm
#' Improvement: use the sd of the pa as weight in MetaRsea
#' 
#' 
#' @author  Junqiang Wang
#' @export


Rsea <- function (eset, regulon, minsize=20) {
  
  if(length(names(regulon))==1){
    regulon<-regulon[1][[1]]}
  
  # Data Preprocessing # Step 1 - Filter the expression data.
  
  # 1.1 Get a list of all targets of every regulator in the network.
  
  tmp <- c(names(regulon), unlist(lapply(regulon, function(x) names(x$tfmode)), use.names = FALSE))
  
  # 1.2 Remove expression data of all the genes which aren't in the list identified above. 
  
  eset <- eset[rownames(eset) %in% unique(tmp),]
  
  # 1.3 Scale the expression matrix on rows [i.e. row = (row - mean(row)) / sd(row) ]. 
  
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
  
  nes <- Rsea2(tt, regulon)

# Return VIPER matrix 
  
  return(nes)

}






# regulon set enrichment analysis
Rsea2 <- function (tt, regulon) {
    
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
    
    t2 <- apply(tt, 2, rank)/(nrow(tt) + 1)

# This line of code is necessary to match the order of genes # for the subsequent matrix multiplication steps.

pos <- match(targets, rownames(tt))

# 2.2 Transform T2 ranks to Gaussian values. 

t2q <- qnorm(filterRowMatrix(t2, pos))

# 2.3 Matrix multiplication. 

# es follow normal distribution with mean=0

wtss<-t(mor*wts)

es <- t(mor * wts) %*% t2q

wtss2<-wtss^2

es.var <- apply(wtss2, 1, sum)

es.sd<-sqrt(es.var)

nes<-apply(es, 2, function(x, es.sd){
  
  x.nes<-x/es.sd
  return(x.nes)
  
}, es.sd=es.sd)

#t.statistic<-t.test(es, alternative = "less", mu = 0)

#nes<-t.statistic$p.value


#nes<-es/sqrt(es.var)

return(nes)

# Step 3 - Calculate the one-tail score matrix.

# 3.1 Calculate the 'T1 Rank' matrix from the 'T2 Rank' matrix. 

# t1 <- abs(t2 - 0.5) * 2 

# t1 <- t1 + (1 - max(t1))/2

# 3.2 Get qnorm values 

# t1q <- qnorm(filterRowMatrix(t1, pos))

# 3.3 Matrix multiplication.

# sum2 <- t((1 - abs(mor)) * wtss) %*% t1q

# Step 4 - Calculate the Three-tail enrichment score.

# 4.1 Extract the signs of the Two-tail enrichment scores 

#ss <- sign(sum1) ss[ss == 0] <- 1

# 4.2 Combine the Two-tail and One-tail enrichment score matrices. 

#sum3 <- (abs(sum1) + sum2 * (sum2 > 0)) * ss

# Step 5 - Calculate the Normalized Enrichment Scores.

# 5.1 For each regulator, calculate an index proportional to the likelihood value # of all its regulatory interactions.

#lwt <- sqrt(colSums(wts^2))

# 5.2 Adjust 3T enrichment scores proportionally to the weights calculated above. 

#nes <- sum3 * lwt

}


#metaRSEA

MetaRsea<-function(eset, regulon.list, wt.amplification=1,  minsize=20){
  
  require(parallel)
  
  pa.list<-mclapply(regulon.list, function(x){
    
    x<-Rsea(eset=eset, regulon = x)
    
  })
  
  message("PAs for each regulome are inferred!")
  
  # fill na with 0

 # pa.list<-mclapply(pa.list, function(x){
    
  #  if(length(names(x))==1){
   #   x<-x[1][[1]]}
 #   x[x=="NA"]<-0
 # })
  
  # generate the hash
  
  #hash.list<-mclapply(pa.list, function(x){
#names(x)    

#  })
  
 # identical(colnames(pa.list[[1]]), colnames(eset))
# generate the augmented matrices
  
  regulators<-lapply(pa.list, FUN=rownames) %>% unlist() %>% unique() %>% sort()
  
  cells<-colnames(eset)


pa.aug.list<-vector(mode="list", length=length(pa.list))

names(pa.aug.list)<- names(pa.list)

# 
# pa.aug.list<-mclapply(pa.aug.list, function(x){
# x<-matrix(0, nrow=length(regulators), ncol=length(cells))
# rownames(x)<-regulators
# colnames(x)<-cells
# })

# assign the values by rows



for (i in 1:length(pa.aug.list)){
  
  x<-matrix(0, nrow=length(regulators), ncol=length(cells))
  rownames(x)<-regulators
  colnames(x)<-cells

  tmp<-match(rownames(pa.list[[i]]), rownames(x))
  
  x[tmp,]<-pa.list[[i]]
  
  pa.aug.list[[i]]<-x
  
}

#pa.ws.list<-lapply(pa.aug.list, function(x){
#  x<-x*abs(x)
#})

#pa.var.list<-lapply(pa.aug.list, function(x){
#  x<-x*x
#})

 # pa.ws<-Reduce('+', pa.ws.list)
 # pa.var<-Reduce('+', pa.var.list)
  
 # pa.nes<-pa.ws/sqrt(pa.var)
  
  # do stouffer integration
  
 # no weight
 # pa.nes<-pa.sum/sqrt(length(pa.aug.list))

# weight the NES by sd of the proteins

 x<-matrix(0, nrow=length(regulators), ncol=length(cells))
 rownames(x)<-regulators
 colnames(x)<-cells
 
 
# 
for(i in 1:nrow(x)){
  
  # the wt of each proteins
  wt.i<-lapply(pa.aug.list, function(x){
    y<-var(x[i,])
  }) %>% unlist() %>% unname()
  
  wt.i<-wt.i*wt.amplification

  for(j in 1:ncol(x)){

   nes.ij<- lapply(pa.aug.list, function(x){
      y<-x[i,j]
    }) %>% unlist() %>% unname()
  
   nes.ij<-sum(nes.ij*wt.i)/sqrt(sum(wt.i^2))

   x[i,j]<-nes.ij
  }
}
 
 pa.nes<-x

  message("The NESs are calculated!")
  
  return(pa.nes)
  
}


