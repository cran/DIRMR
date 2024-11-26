\name{MOPXEM}
\alias{MOPXEM}
\title{
MOPXEM
}
\description{
The MOPXEM method is an improved EM algorithm that combines the monotonic super-relaxation strategy with the PXEM strategy. The main idea of the MOPXEM method is to accelerate the EM algorithm using the ULS strategy, while simultaneously expanding and optimizing the model parameters using the PX-EM strategy.
}
\usage{
MOPXEM(data,df1,omega,maxiter)
}
\arguments{
\item{data}{The real data sets with missing data used in the method}
  \item{df1}{The real data sets used in the method}
  \item{omega}{A variable of this method}
  \item{maxiter}{The maximum number of iterations}
}
\value{
\item{Y01}{The response variable value after projection}
\item{Yhat}{The estimated response variable value after projection}
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
MOPXEM(data,df1,omega=0.15,maxiter=15)
}
