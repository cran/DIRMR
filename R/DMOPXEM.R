#' @param data is a highly correlated data set with missing
#' @param df1 is a highly correlated data set
#' @param M is the number of Blocks
#' @param omega is A variable of this method
#' @param maxiter is the maximum number of iterations
#'
#' @return Y011,Yhat,Ymean,Yhatmean
#' @export
#'
#' @examples
#' DMOPXEM(data,df1,M,omega,maxiter)
DMOPXEM<-function(data,df1,M,omega,maxiter){
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
  dd<-lambda<-sigma2_a<-sigma2_e<-MOEMsigma2_e<-MOEMsigma2_a<-oldMOEMsigma2_e<-oldMOEMsigma2_a<-NULL
  dd[1]=dd.start=1
  oldMOEMsigma2_e=MOEMsigma2_e=sigma2_e[1]=sigma2_e.start=1
  lambda[1] <- 2
  oldMOEMsigma2_a=MOEMsigma2_a=sigma2_a[1]=lambda[1]^2*dd[1]
  Kappa <- list()
  Kappa.change <- LL <- NULL
  Kappa[[1]] <- c(MOEMsigma2_a[1], MOEMsigma2_e[1])
  nm = n/M
  Wm=matrix(rep(0,nm*M),ncol=M)
  wr=matrix(rep(0,M*nm),ncol=nm)
  Yhat=matrix(rep(0,nm,M),nrow=nm,ncol=M)
  Y011=matrix(rep(0,nm,M),nrow=nm,ncol=M)
  Ymean=matrix(rep(0,nm,M),nrow=nm,ncol=M)
  Yhatmean=matrix(rep(0,nm,M),nrow=nm,ncol=M)
  beta=matrix(rep(0,p,M),nrow=p,ncol=M)
  b.tilde=matrix(rep(0,q,M),nrow=q,ncol=M)
  MOEMbeta=matrix(rep(0,p,M),nrow=p,ncol=M)
  MOEMb.tilde=matrix(rep(0,q,M),nrow=q,ncol=M)
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
    K <- diag(nm-l1) - X101%*%solve(t(X101)%*%X101)%*%t(X101)
    W <- cbind(X101,Z101)
    for(w in 1:maxiter)
    {
      G <- sigma2_a[w]*diag(q)
      Ginv <- (1/sigma2_a[w])*diag(q)
      H <- Z101%*%G%*%t(Z101) + sigma2_e[w]*diag(nm-l1)
      Hinv <- solve(H)
      beta[,m]=solve(t(X101)%*%Hinv%*%X101)%*%t(X101)%*%Hinv%*%Y101
      P <- Hinv - Hinv%*%X101%*%solve(t(X101)%*%Hinv%*%X101)%*%t(X101)%*%Hinv
      S <- (1/sigma2_e[w])*K
      Rinv <- (1/sigma2_e[w])*diag(nm-l1)
      C.ZZ <- solve(t(Z101)%*%S%*%Z101+ Ginv)
      C.XZ <- -solve(t(X101)%*%Hinv%*%X101)%*%t(X101)%*%Hinv%*%Z101%*%G
      C_ZZ <- t(Z101)%*%Rinv%*%Z101+ Ginv
      C_XX <- t(X101)%*%Rinv%*%X101
      C_XZ <- t(X101)%*%Rinv%*%Z101
      C <- rbind(cbind(C_XX,C_XZ), cbind(t(C_XZ),C_ZZ))
      Cinv <- solve(C)
      b.tilde[,m]<- G%*%t(Z101)%*%P%*%Y101
      e.tilde<- sigma2_e[w]*P%*%Y101
      LL[w] <- -0.5 * ( determinant(t(X101)%*%Hinv%*%X101, logarithm=TRUE)$modulus + determinant(H, logarithm=TRUE)$modulus + t(Y101)%*%P%*%Y101)
      M1=t(Y101-Z101%*%b.tilde[,m])%*%K%*%(Y101-Z101%*%b.tilde[,m])+sum(diag(t(Z101)%*%K%*%Z101%*%C.ZZ))
      M2=t(Y101)%*%K%*%Z101%*%b.tilde[,m]/(t(b.tilde[,m])%*%t(Z101)%*%K%*%Z101%*%b.tilde[,m]+ sum(diag(t(Z101)%*%K%*%Z101%*%C.ZZ)) )
      sigma2_e[w+1]=1/(nm-l1-p) *M1
      dd[w+1] <- (1/q)*(t(b.tilde[,m])%*%b.tilde[,m] + sum(diag(C.ZZ)))
      lambda[w+1]<-M2
      sigma2_a[w+1] <- lambda[w+1]^2*dd[w+1]
      Kappa[[w+1]] <- c(sigma2_a[w+1],sigma2_e[w+1])
      Kappa.change[w] <- sqrt( (t(Kappa[[w+1]] - Kappa[[w]])%*%(Kappa[[w+1]] - Kappa[[w]])) / (t(Kappa[[w]])%*%Kappa[[w]]) )
      if (Kappa.change[w] < 10^-8) {
        message("Converged in ", w, " iterations")
        break
      }
    }
    rrr=sigma2_a[w+1]/oldMOEMsigma2_a
    omegat=omega*(min(rrr))/(1+omega-omega*(min(rrr)))
    MOEMsigma2_a[w+1]=(1+omegat)*sigma2_a[w+1]-omegat*oldMOEMsigma2_a
    MOEMsigma2_e[w+1]=(1+omega)*sigma2_e[w+1]-omega*oldMOEMsigma2_e
    G0 <- MOEMsigma2_a[w+1]*diag(q)
    Ginv0 <- (1/MOEMsigma2_a[w+1])*diag(q)
    H0 <- Z101%*%G0%*%t(Z101) + MOEMsigma2_e[w+1]*diag(nm-l1)
    Hinv0 <- solve(H0)
    MOEMbeta[,m]=solve(t(X101)%*%Hinv0%*%X101)%*%t(X101)%*%Hinv0%*%Y101
    P0<-Hinv0-Hinv0%*%X101%*%solve(t(X101)%*%Hinv0%*%X101)%*%t(X101)%*%Hinv0
    MOEMb.tilde[,m]<- G0%*%t(Z101)%*%P0%*%Y101
    for (i in 1:l1){
      Y001[i]=X001[i,]%*%MOEMbeta[,m]+Z001[i,]%*%MOEMb.tilde[,m]}
    Yhat[,m]=Y001
    Ymean[,m]=matrix(mean(Y011[,m]),nrow=nm,ncol=1)
    Yhatmean[,m]=matrix(mean(Yhat[,m]),nrow=nm,ncol=1)
  }
  result <- list(Y011 = Y011, Yhat = Yhat, Ymean = Ymean, Yhatmean = Yhatmean)
  return(result)
}
