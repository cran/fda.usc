\name{norm.fdata}
\Rdversion{1.1}
\alias{norm.fdata}
\alias{norm.fd}
\title{
Approximates Lp-norm  for functional data. 
}
\description{	
Approximates Lp-norm for functional data (fdata) object  using metric or semimetric functions. Norm for functional data  using by default Lp-metric.
}
\usage{
norm.fdata(fdataobj,metric=metric.lp,...)
norm.fd(fdobj)
}
\arguments{  
   \item{fdataobj}{ \code{\link{fdata}} class object.}  
  \item{fdobj}{Functional data  or curves of  \code{\link{fd}} class. }
   \item{metric}{ Metric function, by default \code{\link{metric.lp}}.}
   \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{
By default it computes the L2-norm with \code{p = 2} and weights \code{w} with length=\code{(m-1)}. 
\deqn{Let \ \  f(x)= fdataobj(x)\ }{}
\deqn{\left\|f\right\|_p=\left (  \frac{1}{\int_{a}^{b}w(x)dx}  \int_{a}^{b} \left|f(x)\right|^{p}w(x)dx \right)^{1/p}}{}

The observed points on each curve are equally spaced (by default) or not. 
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{
See also  \code{\link{metric.lp}} and \code{\link{norm}}\cr
Alternative method:  \code{\link{inprod}} of fda-package
} 
\examples{


x<-seq(0,2*pi,length=1001)
fx1<-sin(x)/sqrt(pi)
fx2<-cos(x)/sqrt(pi)
argv<-seq(0,2*pi,len=1001)
fdat0<-fdata(rep(0,len=1001),argv,range(argv))
fdat1<-fdata(fx1,x,range(x))
metric.lp(fdat1)
metric.lp(fdat1,fdat0)
norm.fdata(fdat1)
# The same
integrate(function(x){(abs(sin(x)/sqrt(pi))^2)},0,2*pi)
integrate(function(x){(abs(cos(x)/sqrt(pi))^2)},0,2*pi)

bspl1<- create.bspline.basis(c(0,2*pi),21)
fd.bspl1 <- fd(basisobj=bspl1)
fd.bspl2<-fdata2fd(fdat1,nbasis=21)
norm.fd(fd.bspl1)
norm.fd(fd.bspl2)

}
\keyword{ math }