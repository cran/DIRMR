#' @param data is a highly correlated data set with missing
#' @param df1 is a highly correlated data set
#' @param M is the number of Blocks
#' @param maxiter is the maximum number of iterations
#'
#' @return Y011,Yhat,Ymean,Yhatmean
#' @export
#'
#' @examples
#'set.seed(99)
#'library(MASS)
#'library(mvtnorm)
#'n=600;p=3;q=2;M=5;omega=0.15;ratio=0.1;maxiter=1500;nob=round(n-(n*ratio))
#'dd.start=1;sigma2_e.start=1
#'X0=matrix(runif(n*p,0,2),ncol=p)
#'beta=matrix(rnorm(p*1,0,3),nrow=p)
#'Z0=matrix(runif(n*q,2,3),ncol=q)
#'e=matrix(rnorm(n*1,0,sigma2_e.start),n,1)
#'b=matrix(rnorm(q*1,0,1),q,1)
#'Y0=X0%*%beta+Z0%*%b+e
#'df1=data.frame(Y=Y0,X=X0,Z=Z0)
#'misra=function(data,ratio){
#'  nob=round(n-(n*ratio))
#'  data[sample(n,n-nob),1]=NA
#'  return(data)}
#'data=misra(data=df1,ratio=0.1)
#'DEM(data,df1,M,maxiter)
DEM<-function(data,df1,M,maxiter){
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
dd <- lambda <- sigma2_a <- sigma2_e <- NULL
dd[1]=dd.start=1
sigma2_e[1]=sigma2_e.start=1
lambda[1] <- 2
sigma2_a[1] <- lambda[1]^2*dd[1]
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
Kappa[[1]] <- c(betahat0,bhat0)
for(w in 1:maxiter)
{
I=diag(c(rep(1,nm-l1)))
G <- sigma2_a[w]*diag(q)
Ginv <- (1/sigma2_a[w])*diag(q)
H <- Z101%*%G%*%t(Z101) + sigma2_e[w]*diag(nm-l1)
Hinv <- solve(H)
beta[,m] <- solve(t(X101)%*%Hinv%*%X101)%*%t(X101)%*%Hinv%*%Y101
P <- Hinv%*%X101%*%solve(t(X101)%*%Hinv%*%X101)%*%t(X101)%*%Hinv
b.tilde[,m] <- G%*%t(Z101)%*%Hinv%*%(Y101-X101%*%beta[,m])
e.tilde=(I-Z101%*%G%*%t(Z101)%*%Hinv)%*%(I-X101%*%solve(t(X101)%*%Hinv%*%X101)%*%t(X101)%*%Hinv)%*%Y101
sigma2_e[w+1]=sigma2_e[w]+t(e.tilde)%*%e.tilde-(sigma2_e[w])^2*tr(Hinv-P)
sigma2_a[w+1]=sigma2_a[w]+t(b.tilde[,m])%*%b.tilde[,m]-tr(sigma2_a[w]*t(Z101)%*%(Hinv-P)%*%Z101*sigma2_a[w])
Kappa[[w+1]] <- c(beta[,m],b.tilde[,m])
Kappa.change[w] <- sqrt( (t(Kappa[[w+1]] - Kappa[[w]])%*%(Kappa[[w+1]] - Kappa[[w]])) / (t(Kappa[[w]])%*%Kappa[[w]]) )
if (Kappa.change[w] < 10^-8) {
  message("Converged in ", w, " iterations")
  break
}
}
for (i in 1:l1){
Y001[i]=X001[i,]%*%beta[,m]+Z001[i,]%*%b.tilde[,m]}
Yhat[,m]=Y001
Ymean[,m]=matrix(mean(Y011[,m]),nrow=nm,ncol=1)
Yhatmean[,m]=matrix(mean(Yhat[,m]),nrow=nm,ncol=1)
}

result <- list(Y011 = Y011, Yhat = Yhat, Ymean = Ymean, Yhatmean = Yhatmean)
return(result)
}
