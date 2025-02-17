
#' Rosa - regulon structure-based analysis of protein activities
#' @param reg_weigth_combiation methods for combining two set of reg weights: multiplicative or additive
#' @param reg_weight The reg_weight used: mi, cor or both
#' @author  Junqiang Wang
#' @export
Rosa <- function (eset,
                  regulon,
                  reg_weight=c("both", "mi", "cor"),
                  reg_weight_combination=c("multiplicative","additive"),
                  wt.features=NULL,
                  DoScaleEset=FALSE,
                  minsize=20) {

  message("the size of regulon is ", length(names(regulon)))

   if(length(names(regulon))<100){
    regulon<-regulon[[1]]}


  if(is.null(wt.features)==TRUE){

    message("Null wt.features is used!")

    wt.features<-rep(1, times=nrow(eset))

    names(wt.features)<-rownames(eset)
  }


  tmp <- c(names(regulon), unlist(lapply(regulon, function(x) names(x$tfmode)), use.names = FALSE))
  eset <- eset[rownames(eset) %in% unique(tmp),]


   if(isTRUE(DoScaleEset)){
   eset <- t(scale(t(eset)))
   }


  regulon <- lapply(regulon, function(x, genes) {
  filtro <- names(x$tfmode) %in% genes
  x$tfmode <- x$tfmode[filtro]
  if (length(x$likelihood) == length(filtro)){
    x$likelihood <- x$likelihood[filtro] }
   return(x) }, genes = rownames(eset))


  regulon <- regulon[sapply(regulon, function(x) length(x$tfmode)) >= minsize]

  
  nes <- RosaModelStat(eset,
                        regulon,
                        reg_weight=reg_weight,
                        reg_weight_combination=reg_weight_combination,
                        wt.features=wt.features)

  return(nes)

}


RosaModelStat <- function (eset,
                           regulon,
                           reg_weight=c("both", "mi", "cor"),
                           reg_weight_combination=c("multiplicative","additive"),
                           wt.features) {

  if(is.null(wt.features)==TRUE){
    message("Null wt.features is used!")
    wt.features<-rep(1, times=nrow(eset))
    names(wt.features)<-rownames(eset)
  }

    reg_targets <- unique(unlist(lapply(regulon, function(x) names(x$tfmode)), use.names = FALSE))
    reg_targets<-sort(reg_targets)


    
    reg_weight_cor <- sapply(regulon, function(x, genes) {
                  tmp<-x$tfmode[match(genes, names(x$tfmode))]
                  return(tmp)
                    }, genes = reg_targets)

    rownames(reg_weight_cor)<-reg_targets
    reg_weight_cor[is.na(reg_weight_cor)]<-0


    reg_mode<-sign(reg_weight_cor)
    reg_mode[is.na(reg_mode)]<-0


    reg_weight_cor<-abs(reg_weight_cor)

    
    reg_weight_mi <- sapply(regulon, function(x, genes) {
                  tmp <- x$likelihood[match(genes, names(x$tfmode))]
                  return(tmp)
                 }, genes = reg_targets)

    rownames(reg_weight_mi)<-reg_targets
    reg_weight_mi<-reg_weight_mi/max(reg_weight_mi, na.rm=T)

    
   reg_weight_mi[is.na(reg_weight_mi)]<-0

   
  z_exp<-matrix(0, nrow=length(reg_targets), ncol=ncol(eset))
  rownames(z_exp)<-reg_targets
  colnames(z_exp)<-colnames(eset)


common_genes <-intersect(rownames(z_exp), rownames(eset))
common_genes_index_z_exp<-match(common_genes, rownames(z_exp))
common_genes_index_eset<-match(common_genes, rownames(eset))
z_exp[common_genes_index_z_exp,]<-eset[common_genes_index_eset,]


wt.features<-wt.features/max(wt.features)
tmp<-match(reg_targets, names(wt.features))
feature_weight<-wt.features[tmp]
names(feature_weight)<-reg_targets

feature_weight[is.na(feature_weight)]<-0
feature_weight<-round(feature_weight, digits=3)



feature_weight<-matrix(rep(feature_weight, times=ncol(reg_mode)), ncol=ncol(reg_mode))
rownames(feature_weight)<-reg_targets
colnames(feature_weight)<-colnames(reg_mode)



message("calculate the total weights ...")

if(reg_weight=="both" && reg_weight_combination=="multiplicative"){

wtss<-reg_mode*reg_weight_cor*reg_weight_mi*feature_weight
}


if(reg_weight=="both" && reg_weight_combination=="additive"){

  wtss<-reg_mode*0.5*(reg_weight_cor+reg_weight_mi)*feature_weight
}



if(reg_weight=="mi"){

  wtss<-reg_mode*reg_weight_mi*feature_weight
}


if(reg_weight=="cor"){

  wtss<-reg_mode*reg_weight_cor*feature_weight
}



wtss<-t(wtss)

message("the total weights is calculated!")



message("calculate the ES ...")

es <- wtss %*% z_exp

message("the ES is calculated!")


message("calculate the NES ...")

wtss2<-wtss^2
es.var <- apply(wtss2, 1, sum)
es.sd<-sqrt(es.var)


nes<-apply(es, 2, function(x, es.sd){

  x.nes<-x/es.sd
  return(x.nes)

}, es.sd=es.sd)

message("the NES is calculated!")

return(nes)
}




MetaRosa<-function(eset,
                   regulon.list,
                   wt.features,
                   reg_weight=c("both", "mi", "cor"),
                   reg_weight_combination=c("multiplicative","additive"),
                   DoScaleEset=FALSE,
                   minsize=20,
                   integration.method=c("average", "maxAbs")
                   ){

  require(parallel)

  pa.list<-mclapply(regulon.list, function(x){

    y<-Rosa(eset=eset,
            regulon = x,
            reg_weight=reg_weight,
            reg_weight_combination=reg_weight_combination,
            wt.features=wt.features,
            DoScaleEset=DoScaleEset
            )
  })

  message("PAs for each regulome are inferred!")


  message("aggregate pa matrices...")

  regulators<-lapply(pa.list, FUN=rownames) %>% unlist() %>% unique() %>% sort()
  cells<-colnames(eset)

pa.aug.list<-vector(mode="list", length=length(pa.list))
names(pa.aug.list)<- names(pa.list)


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


 if(integration.method=="average"){
 x<-Reduce('+', pa.aug.list)
 pa.nes<-x/length(pa.aug.list)
 }


 if(integration.method=="maxAbs"){
   x_max <- do.call(pmax, pa.aug.list)
   x_min<-do.call(pmin, pa.aug.list)
   x_sign<-sign(x_max+x_min)

   pa.aug.list.abs <- lapply(pa.aug.list, abs)
   x_max_abs<-do.call(pmax, pa.aug.list.abs)
   pa.nes<-x_sign*x_max_abs

 }

 #
  message("The NESs are calculated!")
  return(pa.nes)

}



#' @paramters exprs count matrix
#' @export
WeightFeature<-function(seu,
                      signature.wt.method="vst"
                        ){
  
  DefaultAssay(seu)<-"RNA"

  seu <- FindVariableFeatures(seu,
                              assay="RNA",
                              selection.method = "vst",
                              nfeatures = length(rownames(seu)))

  features.vst<-seu@assays$RNA@meta.features$vst.variance.standardized

  names(features.vst)<-rownames(seu@assays$RNA@meta.features)


  if( signature.wt.method=="rank"){

    features.wt<- rank(features.vst)/length(features.vst)

  }

  if( signature.wt.method=="sd"){

    features.wt<-features.vst/max(features.vst)
  }

  return(features.wt)

}



