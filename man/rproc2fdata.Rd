\name{rproc2fdata}
\Rdversion{1.1}
\alias{rproc2fdata}
\title{ Generate random process of fdata class.}
\description{
Generate Functional data from: Ornstein-Uhlenbeck process, Brownian process, Gaussian process or exponential variogram process.
}
\usage{
rproc2fdata(n,t,mu=rep(0,length(t)),sigma=1,
par.list=list("scale"=1,"theta"=1/3),norm=FALSE,...) 
}
\arguments{
	\item{n}{ Number of functional data to generate.}       
  \item{t}{ Discretization points.} 
  \item{mu}{ a vector giving the means of the discretization poins.}
  \item{sigma}{ if is a character: create a covariance matrix of the variables indicated in the argument,  if is a    positive-definite symmetric matrix specifying the covariance matrix of the variables. If is a integer, create the sigma as: sigma=diag(ncol(fdataobj))*sigma, by default sigma=1.} 
  \item{par.list}{ List of parameter to process.}  
  \item{norm}{=TRUE the norm of random projection is 1.}  
 \item{\dots}{ Further arguments passed to or from other methods.}	    
}
\value{
  Return a \code{fdata} class object.
  }  
\references{
cuales
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\examples{

par(mfrow=c(3,2))
plot(rproc2fdata(200,seq(0,1,len=30),sigma="OU",par.list=list("sca"=1)))
plot(rproc2fdata(200,seq(0,1,len=30),sigma="vexponential"))
plot(rproc2fdata(200,seq(0,1,len=30),sigma=1))
plot(rproc2fdata(200,seq(0,1,len=30),sigma=diag(30)))
plot(rproc2fdata(200,seq(0,1,len=30),sigma="brownian"))
plot(rproc2fdata(200,seq(0,1,len=30),sigma="wiener"))
#plot(rproc2fdata(200,seq(0,1,len=30),sigma="oo")) # this is a error 


}
\keyword{datagen}