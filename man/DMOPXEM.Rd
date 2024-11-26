\name{DMOPXEM}
\alias{DMOPXEM}
\title{
DMOPXEM
}
\description{
In DMOPXEM method, data is allocated to different computing nodes for parallel processing. Each node independently executes the EM algorithm and updates the local model parameters. Then, each node passes the local model parameters to other nodes for the merging and updating of global model parameters.}
\usage{
DMOPXEM(data,df1,M,omega,maxiter)
}
\arguments{
\item{data}{The real data sets with missing data used in the method}
 \item{df1}{The real data sets used in the method}
 \item{M}{The number of Blocks}
  \item{omega}{A variable of this method}
  \item{maxiter}{The maximum number of iterations}
}
\value{\item{Y011}{The response variable value after projection for each block}
\item{Yhat}{The estimated response variable value after projection for each block}
\item{Ymean}{The mean of response variable value after projection for each block}
\item{Yhatmean}{The mean of response variable value after projection for each block}
}
\author{Guangbao Guo,Yu Li}
\examples{
set.seed(99)
library(MASS)
library(mvtnorm)
n=50;p=6;q=5;M=2;omega=0.15;ratio=0.1;maxiter=15;nob=round(n-(n*ratio))
dd.start=1;sigma2_e.start=1
X0=matrix(runif(n*p,0,2),ncol=p)
beta=matrix(rnorm(p*1,0,3),nrow=p)
Z0=matrix(runif(n*q,2,3),ncol=q)
e=matrix(rnorm(n*1,0,sigma2_e.start),n,1)
b=matrix(rnorm(q*1,0,1),q,1)
Y0=X0%*%beta+Z0%*%b+e
df1=data.frame(Y=Y0,X=X0,Z=Z0)
misra=function(data,ratio){
  nob=round(n-(n*ratio))
  data[sample(n,n-nob),1]=NA
  return(data)}
data=misra(data=df1,ratio=0.1)
DMOPXEM(data,df1,M=2,omega=0.15,maxiter=15)
}
