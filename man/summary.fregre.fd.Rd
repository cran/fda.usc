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
\tab -Variability of explicative variables  explained by Principal Components.\cr
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
If \code{ask}=FALSE draw graphs in one window, by default. If \code{ask}=TRUE, draw each graph in a window, waiting to confirm.
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
n= 200;tt= seq(0,1,len=101)
x0<-rproc2fdata(n,tt,sigma="wiener")
x1<-rproc2fdata(n,tt,sigma=0.1)
x<-x0*3+x1
beta = tt*sin(2*pi*tt)^2
fbeta = fdata(beta,tt)
y<-inprod.fdata(x,fbeta)+rnorm(n,sd=0.1)

# Functional regression
res=fregre.pc(x,y,l=c(1:5))
summary(res,3,ask=TRUE)

# res2=fregre.pls(x,y,l=c(1:4))
# summary(res2)

# res3=fregre.pls(x,y)
# summary(res3)

}
\keyword{print}