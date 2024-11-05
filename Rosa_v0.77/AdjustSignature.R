

#' @description This function adjusts the signatures in a matrix.
#' @param  x A matrix with column signatures.
#' @param  n_bins The bins used. The bins can be set to values like 100, 1,000, 10,000, and so on.
#' @author Junqiang Wang
#' @export
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
