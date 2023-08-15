#'@description This function do GSEA based on Zscores.
#'@description The Zscores are generated from the normal distribution transformed from the ranked signature
#'
#'
#'
#'@export
GSEAz<-function(geneset,
                signature,
                symmertrize.signature=TRUE
                ){

  if(symmertrize.signature==TRUE){
  rank.signature<-rank(signature)

  q.rank<-rank.signature/(length(rank.signature)+1)

  z<-qnorm(q.rank)

  }else{

    z<-signature/sd(signature)

  }

    tmp<- intersect(names(z), geneset)

    z.geneset<-z[match(tmp, names(signature))]

    z.stouffer<-sum(z.geneset)/sqrt(length(z.geneset))

    p<-pnorm(z.stouffer, lower.tail=FALSE)

    out<-list(NES=z.stouffer, p=p)

    return(out)

}






