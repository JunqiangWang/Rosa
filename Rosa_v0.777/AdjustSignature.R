

#' @description This function adjusts the signatures of a matrix.
#' @param  x A matrix with column signature.
#' @param  n_bins The bins used. can be 100,1000,10000
AdjustSignatureMatrix<-function(x,
                                n_bins=100){

x<-pnorm(x)
x<-trunc(x*n_bins)/n_bins
x[x==0]<-1/n_bins
x[x==1]<-1-1/n_bins

x<-qnorm(x)

  return(x)

}

# Run example
# z<-rnorm(100, mean=0, sd=1)
# z<-matrix(z, nrow=10)
# z<-AdjustSignatureMatrix(z)
# z<-AdjustSignatureMatrix(z, n_bins=100)
#
