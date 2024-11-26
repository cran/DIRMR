#' @param data is a highly correlated data set with missing
#' @param df1 is a highly correlated data set
#' @param maxiter is the maximum number of iterations
#'
#' @return Y01,Yhat
#' @export
#'
#' @examples
#' ECME(data,df1,maxiter)
ECME<-function(data,df1,maxiter){
  X0=as.matrix(data[,c(2,3,4)])
  Z0=as.matrix(data[,c(5,6)])
  Y0=df1[,1]
  n <- dim(X0)[1]
  p <- dim(X0)[2]
  q <- dim(Z0)[2]
  nna=sum(is.na(data))
  nob=n-nna
  data[is.na(data)]=0
  Y=data[,1]
  delta=Y/Y0
  mr=which(delta==0)
  r=matrix(rbind(1:(n-nob),mr),nrow=n-nob,byrow=T)
  R=matrix(rep(0,(n-nob)*n),ncol=n)
  R[r]=1
  nr=which(delta==1)
  v=matrix(rbind(1:nob,nr),nrow=nob,byrow=T)
  V=matrix(rep(0,nob*n),ncol=n)
  V[v]=1
  J=rbind(R,V)
  Y01=J%*%Y0
  Y1=J%*%Y
  X1=J%*%X0
  Z1=J%*%Z0
  Yobs=V%*%Y
  Ymis=R%*%Y
  Xobs=as.matrix(X1[(n-nob+1):n,])
  Zobs=V%*%Z0
  sigma2_a <- sigma2_e <- NULL
  sigma2_e[1]=sigma2_e.start=1
  sigma2_a[1]=sigma2_a.start=1
  Kappa <- list()
  Kappa.change <- LL <- NULL
  Kappa[[1]] <- c(sigma2_a[1], sigma2_e[1])
  nn <- dim(Xobs)[1]
  pp <- dim(Xobs)[2]
  qq <- dim(Zobs)[2]
  K <- diag(nn) - Xobs%*%ginv(t(Xobs)%*%Xobs)%*%t(Xobs)
  for(w in 1:maxiter)
  {
    G <- sigma2_a[w]*diag(qq)
    Ginv=ginv(G)
    R=sigma2_e[w]*diag(nn)
    Rinv=ginv(R)
    U=ginv(Ginv+t(Zobs)%*%Rinv%*%Zobs)
    H <- Zobs%*%G%*%t(Zobs) + sigma2_e[w]*diag(nn)
    Hinv <- ginv(H)
    P <- Hinv - Hinv%*%Xobs%*%ginv(t(Xobs)%*%Hinv%*%Xobs)%*%t(Xobs)%*%Hinv
    beta <- ginv(t(Xobs)%*%Hinv%*%Xobs)%*%t(Xobs)%*%Hinv%*%Yobs
    sigma2_e[w+1]=(1/nn)*(t(Yobs-Xobs%*%beta)%*%K%*%(Yobs-Xobs%*%beta))
    b=U%*%t(Zobs)%*%P%*%(Yobs-Xobs%*%beta)
    sigma2_a[w+1] <- (1/qq)*(sigma2_a[w]*t(b)%*%b+tr(U))
    Kappa[[w+1]] <- c(sigma2_a[w+1], sigma2_e[w+1])
    Kappa.change[w] <- sqrt( (t(Kappa[[w+1]] - Kappa[[w]])%*%(Kappa[[w+1]] - Kappa[[w]])) / (t(Kappa[[w]])%*%Kappa[[w]]) )
    if (Kappa.change[w] < 10^-8) {
      message("Converged in ", w, " iterations")
      break
    }
  }
  for (i in 1:(n-nob)){
    betahat=beta
    bhat=b
    Y1[i]=X1[i,]%*%betahat+Z1[i,]%*%bhat}
  Yhat=Y1

  result <- list(Y01=Y01, Yhat = Yhat)
  return(result)
}
