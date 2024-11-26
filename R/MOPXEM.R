#' @param data is a highly correlated data set with missing
#' @param df1 is a highly correlated data set
#' @param omega is a variable of this method
#' @param maxiter is the maximum number of iterations
#'
#' @return Y01,Yhat
#' @export
#'
#' @examples
#' MOPXEM(data,df1,omega,maxiter)
MOPXEM<-function(data,df1,omega,maxiter){
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
Zobs=V%*%Z0
dd<-lambda<-sigma2_a<-sigma2_e<-MOEMsigma2_e<-MOEMsigma2_a<-oldMOEMsigma2_e<-oldMOEMsigma2_a<-NULL
dd[1]=dd.start=1
MOEMsigma2_e=sigma2_e[1]=sigma2_e.start=1
lambda[1] <- 4
MOEMsigma2_a=sigma2_a[1]= lambda[1]^2*dd[1]
Kappa <- list()
Kappa.change <- LL <- NULL
Kappa[[1]] <- c(sigma2_a[1], sigma2_e[1])
nn <- dim(Xobs)[1]
p <- dim(Xobs)[2]
q <- dim(Zobs)[2]
K <- diag(nn) - Xobs%*%solve(t(Xobs)%*%Xobs)%*%t(Xobs)
for(w in 1:maxiter)
{
G <- sigma2_a[w]*diag(q)
Ginv <- (1/sigma2_a[w])*diag(q)
H <- Zobs%*%G%*%t(Zobs) + sigma2_e[w]*diag(nn)
Hinv <- solve(H)
beta=solve(t(Xobs)%*%Hinv%*%Xobs)%*%t(Xobs)%*%Hinv%*%Yobs
P <- Hinv - Hinv%*%Xobs%*%solve(t(Xobs)%*%Hinv%*%Xobs)%*%t(Xobs)%*%Hinv
S <- (1/sigma2_e[w])*K
Rinv <- (1/sigma2_e[w])*diag(nn)
C.ZZ <- solve(t(Zobs)%*%S%*%Zobs+ Ginv)
C.XZ <- -solve(t(Xobs)%*%Hinv%*%Xobs)%*%t(Xobs)%*%Hinv%*%Zobs%*%G
C_ZZ <- t(Zobs)%*%Rinv%*%Zobs+ Ginv
C_XX <- t(Xobs)%*%Rinv%*%Xobs
C_XZ <- t(Xobs)%*%Rinv%*%Zobs
C <- rbind(cbind(C_XX,C_XZ), cbind(t(C_XZ),C_ZZ))
Cinv <- solve(C)
b.tilde<- G%*%t(Zobs)%*%P%*%Yobs
e.tilde<- sigma2_e[w]*P%*%Yobs
LL[w] <- -0.5 * ( determinant(t(Xobs)%*%Hinv%*%Xobs,logarithm=TRUE)$modulus + determinant(H, logarithm=TRUE)$modulus + t(Yobs)%*%P%*%Yobs)
M1=t(Yobs-Zobs%*%b.tilde)%*%K%*%(Yobs-Zobs%*%b.tilde)+sum(diag(t(Zobs)%*%K%*%Zobs%*%C.ZZ))
M2=t(Yobs)%*%K%*%Zobs%*%b.tilde/(t(b.tilde)%*%t(Zobs)%*%K%*%Zobs%*%b.tilde+ sum(diag(t(Zobs)%*%K%*%Zobs%*%C.ZZ)) )
sigma2_e[w+1]=1/(nn-p) *M1
dd[w+1] <- (1/q)*(t(b.tilde)%*%b.tilde+sum(diag(C.ZZ)))
lambda[w+1]<-M2
sigma2_a[w+1] <- lambda[w+1]^2*dd[w+1]
Kappa[[w+1]] <- c(sigma2_a[w+1],sigma2_e[w+1])
Kappa.change[w] <- sqrt( (t(Kappa[[w+1]] - Kappa[[w]])%*%(Kappa[[w+1]] - Kappa[[w]])) / (t(Kappa[[w]])%*%Kappa[[w]]) )
if (Kappa.change[w] < 10^-8) {
  message("Converged in ", w, " iterations")
  break
}}
r=sigma2_a[w+1]/MOEMsigma2_a
omegat=omega*(min(r))/(1+omega-omega*(min(r)))
MOEMsigma2_a[w+1]=(1+omegat)*sigma2_a[w+1]-omegat*MOEMsigma2_a
MOEMsigma2_e[w+1]=(1+omega)*sigma2_e[w+1]-omega*MOEMsigma2_e
G0 <- MOEMsigma2_a[w+1]*diag(q)
Ginv0 <- (1/MOEMsigma2_a[w+1])*diag(q)
H0 <- Zobs%*%G0%*%t(Zobs) + MOEMsigma2_e[w+1]*diag(nn)
Hinv0 <- solve(H0)
MOEMbeta=solve(t(Xobs)%*%Hinv0%*%Xobs)%*%t(Xobs)%*%Hinv0%*%Yobs
P0<-Hinv0-Hinv0%*%Xobs%*%solve(t(Xobs)%*%Hinv0%*%Xobs)%*%t(Xobs)%*%Hinv0
MOEMb.tilde<- G0%*%t(Zobs)%*%P0%*%Yobs
for (i in 1:(n-nob)){
Y1[i]=X1[i,]%*%MOEMbeta+Z1[i,]%*%MOEMb.tilde}
Yhat=Y1
result <- list(Y01=Y01, Yhat = Yhat)
return(result)}
