#' @title FracField
#'
#' @description Generates fractional Brownian vector field data
#'
#' @param  K number of observations
#' @param  m number of criteria
#' @param  H Hurst parameter (a real in interval [0,1))
#' @param  X Coordinates

#' @return Returns a fractional Brownian vector field matrix.
#'
#' @examples
#' # Load FracKrigingR library
#' library(FracKrigingR)
#' # generate Coordinates
#'    p=2; K=10;
#'    X<-matrix(0,ncol=p, nrow=K)
#'    for(j in 1:p){
#'      for(i in 1:K){
#'        X[i,j] = rnorm(1, 0, 1)
#'      }
#'    }
#'    # generate fractional Brownian vector field
#'    H = 0.5; m = 3
#'    FracField(K,m,H,X)
#'
#' @importFrom graphics mtext
#'
#' @importFrom points stats rnorm
#' @export




FracField<- function(K,m,H,X){
  suppressWarnings({
   positive_definite_matrix <- (clusterGeneration::genPositiveDefMat(dim = m)$Sigma)
   cholecky_matrix<- t(chol(positive_definite_matrix))

   E<-as.matrix( c(rep(c(1), times = K)))

    gB<-function(H){
      gg<-((E%*%t(E))-((FracMatrix(H,K,X)*c(t(E)%*%solve(FracMatrix(H,K,X))%*%E))/2))
      gg

    }

    vector<- function(m,K){
    ksi<-matrix(0,ncol =m,nrow=K)
    for (l in seq(from =1,to =m, by=1)){ksi[,l]<-rnorm(K,0,1)};
    ksi
    }

    y<-t(cholecky_matrix%*%t(vector(m,K))%*%((chol(gB(H)))))

    return(y)
})
}


