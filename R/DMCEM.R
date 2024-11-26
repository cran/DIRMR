#' @param data is a highly correlated data set with missing
#' @param df1 is a highly correlated data set
#' @param M is the number of Blocks
#' @param maxiter is the maximum number of iterations
#'
#' @return Y011,Yhat,Ymean,Yhatmean
#' @export
#'
#' @examples
#' DMCEM(data,df1,M,maxiter)
DMCEM<-function(data,df1,M,maxiter){
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
  Kappa <- list()
  Kappa.change <- LL <- NULL
  nm = n/M
  Wm=matrix(rep(0,nm*M),ncol=M)
  wr=matrix(rep(0,M*nm),ncol=nm)
  Yhat=matrix(rep(0,nm,M),nrow=nm,ncol=M)
  Y011=matrix(rep(0,nm,M),nrow=nm,ncol=M)
  Ymean=matrix(rep(0,nm,M),nrow=nm,ncol=M)
  Yhatmean=matrix(rep(0,nm,M),nrow=nm,ncol=M)
  beta=matrix(rep(0,p,M),nrow=p,ncol=M)
  b.tilde=matrix(rep(0,q,M),nrow=q,ncol=M)
  for (m in 1:M) {
    wr[m,] = sample(1:n,nm,replace=TRUE)
    w=matrix(c(1:nm,wr[m,]),ncol=nm,byrow=T)
    Wm[,m]=w[2,]
    W=matrix(rep(0,nm*n),ncol=n)
    W[t(w)]=1
    Y111=W%*%Y0
    Y100=W%*%Y
    X100=W%*%X0
    Z100=W%*%Z0
    delta=Y100/Y111
    mr=which(delta==0)
    l1=length(mr)
    r=matrix(rbind(1:l1,mr),nrow=l1,byrow=T)
    R=matrix(rep(0,l1*nm),ncol=nm)
    R[r]=1
    nr=which(delta==1)
    l2=length(nr)
    v=matrix(rbind(1:l2,nr),nrow=l2,byrow=T)
    V=matrix(rep(0,l2*nm),ncol=nm)
    V[v]=1
    J=rbind(R,V)
    Y011[,m]=J%*%Y111
    Y001=J%*%Y100
    X001=J%*%X100
    Z001=J%*%Z100
    Y101=V%*%Y100
    X101=V%*%X100
    Z101=V%*%Z100
    C=cbind(X101,Z101)
    thetahat0=solve(t(C)%*%C)%*%t(C)%*%Y101
    betahat0=thetahat0[1:p]
    bhat0=thetahat0[(p+1):(p+q)]
    Xbe=c(rep(0,l1))
    betahat=betahat0
    bhat=bhat0
    Kappa[[1]] <- c(betahat,bhat)
    for(w in 1:maxiter)
    {
      for (i in 1:(l1)){
        Y001[i]=X001[i,]%*%betahat+Z001[i,]%*%bhat
        Xbe=Y001-Z001%*%bhat
        D=cbind(X001,Z001)
        IP=diag(c(rep(1,p)))
        RP=matrix(0,p,q)
        M=cbind(IP,RP)
        betahat=M%*%(solve(t(D)%*%D)%*%t(D)%*%Y001)
        bhat=solve(t(Z001)%*%Z001)%*%t(Z001)%*%(Y001-X001%*%betahat)}
      Kappa[[w+1]] <- c(betahat,bhat)
      Kappa.change[w] <- sqrt( (t(Kappa[[w+1]] - Kappa[[w]])%*%(Kappa[[w+1]] - Kappa[[w]])) / (t(Kappa[[w]])%*%Kappa[[w]]) )
      if (Kappa.change[w] < 10^-8) {
        message("Converged in ", w, " iterations")
        break
      }
    }
    Yhat[,m]=Y001
    Ymean[,m]=matrix(mean(Y011[,m]),nrow=nm,ncol=1)
    Yhatmean[,m]=matrix(mean(Yhat[,m]),nrow=nm,ncol=1)
  }
  result <- list(Y011 = Y011, Yhat = Yhat, Ymean = Ymean, Yhatmean = Yhatmean)
  return(result)
}
