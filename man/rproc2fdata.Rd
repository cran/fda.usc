\name{rproc2fdata}
\Rdversion{1.1}
\alias{rproc2fdata}
\title{ Generate random process of fdata class.}
\description{
Generate Functional data from: Ornstein Uhlenbeck process, Brownian process, Gaussian process or Exponential variogram process.
}
\usage{
rproc2fdata(n,t=NULL,mu=rep(0,length(t)),sigma=1,
par.list=list("scale"=1,"theta"=.2*diff(rtt)),norm=FALSE,verbose=FALSE,...)
}
\arguments{
	\item{n}{ Number of functional data to generate.}       
  \item{t}{ Discretization points.} 
\item{mu}{ \code{vector} which specifies the trend values at the discretization points, by default \code{mu}=\eqn{\mu(t)=0}.
 If \code{mu} is a \code{fdata} class object, \code{t}\eqn{=}\code{argvals(mu)}.}  
\item{sigma}{ 
 A positive-definite symmetric matrix, \eqn{\Sigma_{s,t}}{\Sigma}, specifying the covariance matrix of the variables. 
 If \code{sigma} is a \code{scalar}, creates a random Gaussian process with  \eqn{\Sigma_{s,t}=}{\Sigma=}\code{sigmaI} (by default \code{sigma=1}).\cr
 If  \code{sigma} is a \code{vector}, creates a random Gaussian process with \eqn{\Sigma_{s,t}=}{\Sigma=}\code{diag(sigma)}.\cr
If \code{sigma} is a character: create a random process using the covariance matrix \eqn{\Sigma_{s,t}}{\Sigma} indicated in the argument,
\itemize{              
 \item \code{"OU"} or \code{"OrnsteinUhlenbeck"}, creates an random Ornstein Uhlenbeck process with
  \eqn{\Sigma_{s,t}=\frac{\sigma^2}{2\theta}e^{-\theta\left(s+t\right)} \left(e^{2\theta\left(s+t\right)}-1\right)}{\Sigma=\sigma^2/\theta exp(-\theta(s+t))(exp(2\theta(s+t)-1))},
   by default \eqn{\theta=1/(3range(t))}{\theta=1/(3range(t))}, \eqn{\sigma^2={1}}{\sigma^2=1}.
 \item \code{"wiener"}, creates an random Wiener process with \eqn{\Sigma_{s,t}=\sigma^2 min(s,t)}{\Sigma=\sigma^2 min(s,t)}, by default \eqn{\sigma^2=1}{\sigma^2=1}.
 \item \code{"vexponential"}, creates an random Exponential variogram process with
  \eqn{\Sigma_{s,t}=\sigma^2 e^{\left(-\frac{\left|s-t\right|}{\theta}\right)}}{\Sigma=\sigma^2 exp(-|s-t|/\theta )}, 
  by default \eqn{\theta={0.2 range(t)}}{\theta=.2*range(t)}, \eqn{\sigma^2={1}}{\sigma^2=1}.
}}
 \item{par.list}{ List of parameter to process, by default \code{"scale"} \eqn{\sigma^2=1}{\sigma^2=1} and \code{"theta"} \eqn{\theta=0.2 range(t)}{\theta=0.2 range(t)}.}  
 \item{norm}{ If \code{TRUE} the norm of random projection is 1. Default is \code{FALSE}}  
 \item{verbose}{ If \code{TRUE}, information about procedure is printed. Default is \code{FALSE}.} 
 \item{\dots}{ Further arguments passed to or from other methods.}	    
}
\value{
  Return the functional random process of \code{fdata} class object.
  }  
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\examples{
par(mfrow=c(3,2))
lent<-30
tt<-seq(0,1,len=lent)
mu<-fdata(rep(0,lent),tt)
plot(rproc2fdata(200,t=tt,sigma="OU",par.list=list("scale"=1)))
plot(rproc2fdata(200,mu=mu,sigma="OU",par.list=list("scale"=1)))
plot(rproc2fdata(200,t=tt,sigma="vexponential"))
plot(rproc2fdata(200,t=tt,sigma=1:lent))
plot(rproc2fdata(200,t=tt,sigma="brownian"))
plot(rproc2fdata(200,t=tt,sigma="wiener"))
#plot(rproc2fdata(200,seq(0,1,len=30),sigma="oo")) # this is a error 
}
\keyword{datagen}