#' @param data is a highly correlated data set with missing
#' @param df1 is a highly correlated data set
#' @param maxiter is the maximum number of iterations
#'
#' @return Y01,Yhat
#' @export
#'
#' @examples
#' PXEM(data,df1,maxiter)
PXEM<-function(data,df1,maxiter){
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
dd <- lambda <- sigma2_a <- sigma2_e <- NULL
dd[1]=dd.start=1
sigma2_e[1]=sigma2_e.start=1
lambda[1] <- 1
sigma2_a[1] <- lambda[1]^2*dd[1]
Kappa <- list()
Kappa.change <- LL <- NULL
Kappa[[1]] <- c(sigma2_a[1], sigma2_e[1])
nn <- dim(Xobs)[1]
pp <- dim(Xobs)[2]
qq <- dim(Zobs)[2]
K <- diag(nn) - Xobs%*%ginv(t(Xobs)%*%Xobs)%*%t(Xobs)
W <- cbind(Xobs,Zobs)
for(w in 1:maxiter)
{
G <- sigma2_a[w]*diag(qq)
Ginv <- (1/sigma2_a[w])*diag(qq)
H <- Zobs%*%G%*%t(Zobs) + sigma2_e[w]*diag(nn)
Hinv <- ginv(H)
beta <- ginv(t(Xobs)%*%Hinv%*%Xobs)%*%t(Xobs)%*%Hinv%*%Yobs
P <- Hinv - Hinv%*%Xobs%*%ginv(t(Xobs)%*%Hinv%*%Xobs)%*%t(Xobs)%*%Hinv
S <- (1/sigma2_e[w])*K
Rinv <- (1/sigma2_e[w])*diag(nn)
C.ZZ <- ginv(t(Zobs)%*%S%*%Zobs+ Ginv)
C.XZ <- -ginv(t(Xobs)%*%Hinv%*%Xobs)%*%t(Xobs)%*%Hinv%*%Zobs%*%G
C_ZZ <- t(Zobs)%*%Rinv%*%Zobs+ Ginv
C_XX <- t(Xobs)%*%Rinv%*%Xobs
C_XZ <- t(Xobs)%*%Rinv%*%Zobs
C <- rbind(cbind(C_XX,C_XZ), cbind(t(C_XZ),C_ZZ))
Cinv <- ginv(C)
b.tilde <- G%*%t(Zobs)%*%P%*%Yobs
e.tilde <- sigma2_e[w]*P%*%Yobs
LL[w] <- -0.5 * ( determinant(t(Xobs)%*%Hinv%*%Xobs, logarithm=TRUE)$modulus + determinant(H, logarithm=TRUE)$modulus + t(Yobs)%*%P%*%Yobs)
sigma2_e[w+1] <- 1/(nn-pp) * ( t(Yobs-Zobs%*%b.tilde)%*%K%*%(Yobs-Zobs%*%b.tilde) + sum(diag(t(Zobs)%*%K%*%Zobs%*%C.ZZ)) )
dd[w+1] <- (1/qq)*(t(b.tilde)%*%b.tilde + sum(diag(C.ZZ)))
lambda[w+1]<-t(Yobs)%*%K%*%Zobs%*%b.tilde/(t(b.tilde)%*%t(Zobs)%*%K%*%Zobs%*%b.tilde + sum(diag(t(Zobs)%*%K%*%Zobs%*%C.ZZ)) )
sigma2_a[w+1] <- lambda[w+1]^2*dd[w+1]
Kappa[[w+1]] <- c(sigma2_a[w+1], sigma2_e[w+1])
Kappa.change[w] <- sqrt( (t(Kappa[[w+1]] - Kappa[[w]])%*%(Kappa[[w+1]] - Kappa[[w]])) / (t(Kappa[[w]])%*%Kappa[[w]]) )
if (Kappa.change[w] < 10^-8) {
  message("Converged in ", w, " iterations")
  break
}}
for (i in 1:(n-nob)){
betahat=beta
bhat=b.tilde
Y1[i]=X1[i,]%*%betahat+Z1[i,]%*%bhat}
Yhat=Y1
result <- list(Y01=Y01, Yhat = Yhat)
return(result)}
