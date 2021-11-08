#' @title  FracMatrix
#'
#' @description  Fractional distance matrix
#'
#' @param    H Hurst parameter (a real in interval [0,1))
#' @param    K number of observations
#' @param    X Coordinates
#'
#' @return   Returns a fractional distance matrix, which depends on the Hurst parameter.
#'
#' @examples
#' # Load FracKrigingR library
#' library(FracKrigingR)
#' #Fractional Brownian vector field
#'     K = 10; H = 0.5; p = 2
#' #Generate coordinates
#'     X<-matrix(0,ncol=p, nrow=K)
#'     for(j in 1:p){
#'         for(i in 1:K){
#'             X[i,j] = rnorm(1, 0, 1)
#'         }
#'     }
#'     FracMatrix(H, K, X)
#' @export




FracMatrix<-function(H, K, X){
  a<-matrix(0,nrow=K,ncol=K)
  for (i in 1:K){
    for (j in 1:K){
      a[i,j]<- ((t(X[i,]-X[j,]))%*%((X[i,])-X[j,]))^H
    }
  }
  a
}
