#' @title MaxLikelihood
#'
#' @description  Maximum likelihood method for Hurst parameter estimation of multivariate data
#'
#' @param X Coordinates
#' @param Z Observations
#'
#' @return Returns the estimate of the Hurst parameter (a real in [0,1))
#'  and a graph indicating the minimized maximum likelihood function with the Hurst parameter.
#' @examples
#' # Load FracKrigingR library
#' library(FracKrigingR)
#' # generate Coordinates
#'    p<-2; K<-20;
#'    X<-matrix(0,ncol=p, nrow=K)
#'    for(j in 1:p){
#'      for(i in 1:K){
#'        X[i,j] = rnorm(1, 0, 1)
#'      }
#'    }
#'    # generate fractional Brownian vector field
#'    H <- 0.8; m <- 3
#'    Z<-FracField(K,m,H,X)
#'   # Hurst parameter estimation
#'    MaxLikelihood(X,Z)
#'
#' @export

MaxLikelihood<- function(X,Z){
  Z<-as.matrix(Z)
  K = nrow(X)
  m = ncol(Z)
  E<-c(rep(1,K))
  E<-as.matrix(E)


  golden.section.search = function(f, lower.bound, upper.bound, tolerance)
  {
    golden.ratio = 2/(sqrt(5) + 1)

    x1 = upper.bound - golden.ratio*(upper.bound - lower.bound)
    x2 = lower.bound + golden.ratio*(upper.bound - lower.bound)

    f1 = f(x1)
    f2 = f(x2)

    ite = 0
    while (abs(upper.bound - lower.bound) > tolerance)
    {
      ite = ite + 1
      if (f2 > f1)
      {
        upper.bound = x2
        x2 = x1
        f2 = f1
        x1 = upper.bound - golden.ratio*(upper.bound - lower.bound)
        f1 = f(x1)
      }
      else
      {
        lower.bound = x1
        x1 = x2
        f1 = f2
        x2 = lower.bound + golden.ratio*(upper.bound - lower.bound)
        f2 = f(x2)
      }
    }
    minn = (lower.bound + upper.bound)/2
    minn
  }




  A<-function(H){

    a<-matrix(0,nrow=K,ncol=K)
    for (i in 1:K){
      for (j in 1:K){
        a[i,j]<- ((t(X[i,]-X[j,]))%*%((X[i,])-X[j,]))^H
      }
    }
    a
  }

  bc<-function(H){
    bc<-((1/K)*(((((t(Z)%*%solve(A(H))%*%E)%*%t((t(Z)%*%solve(A(H))%*%E)))/c(t(E)%*%solve(A(H))%*%E)))-(t(Z)%*%solve(A(H))%*%Z)))
    bc
  }

  Ff<-function(H){
    f<-(log(det(bc(H)))/m)+log(abs(det(A(H))))/K
    f
  }


  t<-seq(from=0.01,to=0.9,by=0.01)
  vect<-c()
  for (i in t){
    vect<-c(vect,Ff(i))
  }
  gsc<-golden.section.search(Ff, 0.01, 0.9, 1e-5)
  ff<-Ff(gsc)

  plot(t,vect,type="l", col="darkblue", lwd=1, xlab="", ylab="",cex.axis=1.3)
  points(gsc,ff,pch=8,col="darkgreen",cex=2)
  mtext(text = "Hurst Parameter", side = 1,line = 2.4,font=1,cex=1.4)
  mtext(text = "Maximum Likelihood", side = 2, line = 2.4,font=1,cex=1.4)
  gsc
}
