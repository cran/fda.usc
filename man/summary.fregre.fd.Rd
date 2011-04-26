\name{summary.fregre.fd}
\Rdversion{1.1}
\alias{summary.fregre.fd}
\alias{print.fregre.fd}

\title{ Summarizes information from fregre.fd objects. }
\description{
	Summary function for \code{\link{fregre.pc}}, \code{\link{fregre.basis}}, \code{\link{fregre.pls}}, \code{\link{fregre.np}} and \code{\link{fregre.plm}} functions.
}
\usage{
\method{summary}{fregre.fd}(object,times.influ=3,times.sigma=3,draw=TRUE,\dots)
\method{print}{fregre.fd}(x, digits = max(3, getOption("digits") - 3),\dots)
}
\arguments{
  \item{object,x}{ Estimated by functional regression, \code{fregre.fd} object.}
	\item{times.influ}{ Limit for detect possible infuence curves.}
	\item{times.sigma}{ Limit for detect possible oultiers or atypical curves.}
	\item{draw}{ =TRUE draw   estimation and residuals graphics.}
	\item{digits}{a non-null value for digits specifies the minimum number of significant digits to be printed in values. The default, NULL, uses \code{\link{getOption}}(digits).}	

  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{
Shows:\cr
\tabular{ll}{
\tab -Call.\cr
\tab -R squared.\cr
\tab -Residual variance.\cr
\tab -Index of possible atypical curves or possible outliers.\cr
\tab -Index of possible influence curves.\cr
}

If the \code{fregre.fd} object comes from the \code{\link{fregre.pc}} then shows:
\tabular{ll}{
\tab -Variability of explicative variables  explaines by Principal Components.\cr
\tab -Variability for each  principal components -PC-.\cr
}
     
If draw=TRUE plot: \cr
\tabular{ll}{
\tab -y vs y fitted values.\cr
\tab -Residuals vs fitted values.\cr
\tab -Standarized residuals vs fitted values.\cr
\tab -Levarage.\cr
\tab -Residual boxplot.\cr
\tab -Quantile-Quantile Plot (qqnorm).\cr
}
If \code{ask}=FALSE draw graphs in one window, by defalult. If \code{ask}=TRUE, draw each graph in a window, waiting to confirm.
}
\value{
\item{Influence}{ Vector of influence measures.}
\item{i.influence}{ Index of possible influence curves.}
\item{i.atypical}{ Index of possible atypical curves or possible outliers.}
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
Summary function for  \code{\link{fregre.pc}}, \code{\link{fregre.basis}}, \code{\link{fregre.pls}},  \code{\link{fregre.np}} and \code{\link{fregre.plm}}.
}
\examples{

# Ex 1. Simulated data
library(e1071)
n=500
t=seq(0,1,len=101)
beta=t*sin(2*pi*t)^2 # beta con 2 colinas
x=matrix(NA,ncol=101,nrow=n)
y=numeric(n)
for(i in 1:n){
	x[i, ]=rwiener(1,101)*3+rnorm(101,sd=0.1)
	y[i]=mean(x[i,]*beta)+ rnorm(1,sd=0.1)
}
x.fdata<-fdata(x,t)

# Functional regression
res=fregre.pc(x.fdata,y,l=c(1:5))
summary.fregre.fd(res)
summary.fregre.fd(res,3)
summary.fregre.fd(res,3,draw=FALSE)
summary.fregre.fd(res,3,ask=TRUE)

res2=fregre.pls(x.fdata,y,l=c(1:5))
summary.fregre.fd(res2)
}
\keyword{print}