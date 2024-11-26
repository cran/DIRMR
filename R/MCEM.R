  #' @param data is a highly correlated data set with missing
  #' @param df1 is a highly correlated data set
  #' @param maxiter is the maximum number of iterations
  #'
  #' @return Y01,Yhat
  #' @export
  #'
  #' @examples
  #' MCEM(data,df1,maxiter)
  MCEM<-function(data,df1,maxiter){
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
    Xobs=X1[(n-nob+1):n,]
    Xmis=X1[1:(n-nob),]
    Zobs=V%*%Z0
    Zmis=R%*%Z0
    C=cbind(Xobs,Zobs)
    thetahat0=ginv(t(C)%*%C)%*%t(C)%*%Yobs
    betahat0=thetahat0[1:p]
    bhat0=thetahat0[(p+1):(p+q)]
    Xbe=c(rep(0,n-nob))
    betahat=betahat0
    bhat=bhat0
    Kappa <- list()
    Kappa.change <- LL <- NULL
    Kappa[[1]] <- c(betahat,bhat)
    for(w in 1:maxiter)
    {
      for (i in 1:(n-nob)){
        Y1[i]=X1[i,]%*%betahat+Z1[i,]%*%bhat
        Xbe=Y1-Z1%*%bhat
        D=cbind(X1,Z1)
        IP=diag(c(rep(1,p)))
        RP=matrix(0,p,q)
        M=cbind(IP,RP)
        betahat=M%*%(solve(t(D)%*%D)%*%t(D)%*%Y1)
        bhat=ginv(t(Z1)%*%Z1)%*%t(Z1)%*%(Y1-X1%*%betahat)}
      Kappa[[w+1]] <- c(betahat,bhat)
      Kappa.change[w] <- sqrt( (t(Kappa[[w+1]] - Kappa[[w]])%*%(Kappa[[w+1]] - Kappa[[w]])) / (t(Kappa[[w]])%*%Kappa[[w]]) )
      if (Kappa.change[w] < 10^-8) {
        message("Converged in ", w, " iterations")
        break
      }}
    Yhat=Y1
    return(list(Y01=Y01,Yhat=Yhat))}
