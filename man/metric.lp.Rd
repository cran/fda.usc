\name{metric.lp}
\Rdversion{1.1}
\alias{metric.lp}
\title{
Aproximates Lp-metric distances for functional data. 
}
\description{	
Measures the  proximity between the functional data and curves approximating Lp-metric.
If \code{w = 1} approximates the Lp-metric by Simpson's rule.	By default it uses \code{p = 2} and weights \code{w = 1}.
}
\usage{
metric.lp(fdata1,fdata2=NULL,p=2,w=1,...)
}
\arguments{  
  \item{fdata1}{Functional data 1 or curve 1. If \code{fdata} class, the dimension of \code{fdata1$data} object is (\code{n1} x \code{m}), where \code{n1} is the number of curves and \code{m} are the points observed in each curve.}
  \item{fdata2}{Functional data 2 or curve 2. If \code{fdata} class, the dimension of \code{fdata2$data} object is (\code{n2} x \code{m}), where \code{n2} is the number of curves and \code{m} are the points observed in each curve.}
  \item{p}{ Lp norm, by default it uses \code{p = 2}}
  \item{w}{ Vector of weights with length \code{m}, If \code{w = 1} approximates the metric Lp by Simpson's rule. By default it uses \code{w = 1} }
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{
By default it uses the L2-norm with \code{p = 2}. 
\deqn{Let \ \  f(x)= fdata1(x)-fdata2(x)} {}
\deqn{\left\|f\right\|_p=\left (  \frac{1}{\int_{a}^{b}w(x)dx} \int_{a}^{b} \left|f(x)\right|^{p}w(x)dx \right)^{1/p}}{}

The observed points on each curve are equally spaced (by default) or not. 
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See also  \code{\link{semimetric.basis}} and \code{\link{semimetric.NPFDA}}
} 
\examples{

#	INFERENCE PHONDAT
data(phoneme)
mlearn<-phoneme$learn
mtest<-phoneme$test
glearn<-phoneme$classlearn
gtest<-phoneme$classtest
# Matrix of distances of curves of DATA1
mdist1<-metric.lp(mlearn)

# Matrix of distances beetween curves of DATA1 and curves of DATA2
mdist2<-metric.lp(mlearn,mtest,p=2)
# mdist with L1 norm and weigth=v
v=dnorm(seq(-3,3,len=dim(mlearn)[2]))
mdist3<-metric.lp(mlearn,mtest,p=1,w=v)
plot(1:250,mdist2[1,],type="l",ylim=c(1,max(mdist3[1,])))
lines(mdist3[1,],type="l",col="2")

# mdist with mlearn with diferents discretization points.
#mlearn2=mlearn
#mlearn2[["argvals"]]=seq(0,1,len=150)
#mdist5<-metric.lp(mlearn,mlearn2)
#mdist6<-metric.lp(mlearn2,mlearn) 
#sum(mdist5-mdist6)
#sum(mdist1-mdist6)

x<-seq(0,2*pi,length=1001)
fx<-fdata(sin(x)/sqrt(pi),x)
fx0<-fdata(rep(0,length(x)),x)
metric.lp(fx,fx0)
# The same
integrate(function(x){(abs(sin(x)/sqrt(pi))^2)},0,2*pi)
}
\keyword{ cluster }