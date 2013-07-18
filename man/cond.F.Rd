\name{cond.F}
\Rdversion{1.1}
\alias{cond.F}
\title{ Conditional Distribution Function }
\description{
 Calculate the conditional distribution function of a scalar response with functional data.
}
\usage{
cond.F(fdata0,y0,fdataobj,y,h=0.15,g=0.15,metric=metric.lp,
Ker=list(AKer=AKer.epa,IKer=IKer.epa),...)
}
\arguments{
   \item{fdata0}{ Conditional explanatory functional data of \code{\link{fdata}} class.}  
	\item{y0}{ Vector of conditional response with length \code{n}.}
	\item{fdataobj}{ \code{\link{fdata}} class object.}  
   \item{y}{ Vector of scalar response with length \code{nn}.}	
	\item{h}{ Smoothing parameter or bandwidth of response \code{y}.}
  \item{g}{ Smoothing parameter or bandwidth of explanatory functional data \code{fdataobj}.}
  \item{metric}{ Metric function, by default \code{\link{metric.lp}}.}
  \item{Ker}{ List of 2 arguments. The fist argument is a character string that determines the type of asymetric kernel (see \code{\link{Kernel.asymmetric}}). Asymmetric Epanechnikov kernel is selected by default. The second argumentis a string that determines the type of integrated kernel(see \code{\link{Kernel.integrate}}). Integrate Epanechnikov kernel is selected by default.\cr.}
    \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{ If x.dist=NULL the distance matrix between \code{fdata} objects is calculated by function passed in \code{metric} argument.}
\value{
  \item{Fc}{ Conditional distribution function.}
  \item{y0}{ Vector of conditional response.}
  \item{g}{ Smoothing parameter or bandwidth of explanatory functional data (\code{fdataobj}).}
  \item{h}{ Smoothing parameter or bandwidth of respone, \code{y}.}
  \item{x.dist}{ Distance matrix between curves of \code{fdataobj} object.}
 \item{xy.dist}{ Distance matrix between cuves of \code{fdataobj} and \code{fdata0} objects.}
}
\references{
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional data analysis.}
Springer Series in Statistics, New York.
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{
See Also as:  \code{\link{cond.mode}} and \code{\link{cond.quantile}}. 
}
\examples{

# Read data
n= 500
t= seq(0,1,len=101)
beta = t*sin(2*pi*t)^2
x = matrix(NA, ncol=101, nrow=n)
y=numeric(n)
x0<-rproc2fdata(n,seq(0,1,len=101),sigma="wiener")
x1<-rproc2fdata(n,seq(0,1,len=101),sigma=0.1)
x<-x0*3+x1
fbeta = fdata(beta,t)
y<-inprod.fdata(x,fbeta)+rnorm(n,sd=0.1)

prx=x[1:100];pry=y[1:100]
ind=101;ind2=102:110
pr0=x[ind];pr10=x[ind2,]
ndist=61
gridy=seq(-1.598069,1.598069, len=ndist)
# Conditional Function
res1 = cond.F(pr10, gridy, prx, pry,p=1)
# res2 = cond.F(pr10, gridy, prx, pry,h=0.3)
# res3 = cond.F(pr10, gridy, prx, pry,g=0.25,h=0.3)

# plot(res1$Fc[,1],type="l",ylim=c(0,1))
# lines(res2$Fc[,1],type="l",col=2)
# lines(res3$Fc[,1],type="l",col=3)




}
\keyword{ distribution }
