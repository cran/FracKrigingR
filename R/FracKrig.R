#' @title FracKrig
#'
#' @description Performs extrapolation for spatial multivariate data
#'
#' @param   X Coordinates
#' @param   Z observations
#' @param   Xnew Coordinates of points where the prognosis should be made
#' @param   H Hurst parameter (a real in interval [0,1))
#'
#' @return   Returns a matrix of fractional kriging prognosis.
#'
#'
#' @examples
#'
#' library(sp)
#' library(gstat)
#'  data(meuse)
#'  xy<-cbind(meuse$x,meuse$y)

#'  X<-xy[1:50,]
#'  min_max_norm <- function(x) {
#'      (x - min(x)) / (max(x) - min(x))
#'  }
#'  normalize <- function(x) {
#'  return ((x - min(x)) / (max(x) - min(x)))
#'  }


#'  dat<-cbind(meuse[3],meuse[4],meuse[5])
#'  data<-dat[51:100,]
#'  zz1 <- as.data.frame(lapply(dat, normalize))
#'  data1=as.data.frame(lapply(as.data.frame(data), normalize))
#'  Z<-as.matrix(zz1[1:50,])

# Load FracKrigingR library
#' library(FracKrigingR)
#'  K<-50

#' #Hurst parameter estimation
#'  H<-0.2
#'  Xnew<-xy[51:100,]
#'  results<- FracKrig(X,Z,Xnew,H)
#'  denormalize <- function(x, bottom, top){
#'     (top - bottom) * x + bottom
#'  }

#'z1 = denormalize(
#'  results[,1], top = max(data[,1]), bottom = min(data[,1])
#')
#'z2 = denormalize(
#' results[,2], top = max(data[,2]), bottom = min(data[,2])
#')
#'z3 = denormalize(
#'  results[,3], top = max(data[,3]), bottom = min(data[,3])
#')


#'RMSE<-function(z,prognosis){
#'  rmse<-sqrt(((1/(length(z))))*sum((z-prognosis)^2))
#'  rmse
#'}

#'Cd<-RMSE(data[,1],z1)
#'Cu<-RMSE(data[,2],z2)
#'Pb<-RMSE(data[,3],z3)


#'Cd
#'Cu
#'Pb
#'
#' @export
#'


FracKrig<- function(X,Z,Xnew,H){
  X<-as.matrix(X)
  Z<-as.matrix(Z)
  Xnew<-as.matrix(Xnew)
  E<-c(rep(1,K))
  E<-as.matrix(E)
  K<-nrow(X)
  p<-ncol(X)
  P<-nrow(Xnew)
  m<-ncol(Z)

  prognosis<-function(xpr,H){

    b<-c()
    for (i in seq(from =1, to=K, by=1)){
      b<-c(b,(((xpr - t(X[i,]))%*%t(xpr-t(X[i,])))^H))
    }

    prog<-(as.matrix(t(Z)%*%solve(FracMatrix(H,K,X))))%*%((as.matrix(b))+(E%*%(((1-c(t(E)%*%solve(FracMatrix(H,K,X))%*%(as.matrix(b))))/c(t(E)%*%solve(FracMatrix(H,K,X))%*%E)))))
    prog

  }


  final_result = matrix(c(0), nrow= P, ncol =m)
  for (i in seq(from =1, to=P, by=1)){

    final_result[i,]<-t(prognosis(Xnew[i,],H))
  }
  final_result
}
