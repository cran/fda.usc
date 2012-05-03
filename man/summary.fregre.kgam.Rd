\name{summary.fregre.kgam}
\Rdversion{1.1}
\alias{summary.fregre.kgam}
\alias{print.fregre.kgam}

\title{ Summarizes information from fregre.kgam objects. }
\description{
	Summary function for \code{\link{fregre.kgam}} function.
}
\usage{
\method{summary}{fregre.kgam}(object,draw=TRUE,selec = NULL, times.influ = 3,\dots)
\method{print}{fregre.kgam}(x, digits = max(3, getOption("digits") - 3),\dots)
}
\arguments{
  \item{object,x}{ Estimated by functional regression, \code{fregre.fd} object.}
	\item{draw}{ =TRUE draw   estimation and residuals graphics.}
       \item{selec}{ Allows the plot for a single model term to be selected for printing. e.g. if you just want the plot for the second smooth term set select=2. 
.}
	\item{times.influ}{ Limit for detect possible infuence curves.}

	\item{digits}{a non-null value for digits specifies the minimum number of significant digits to be printed in values. The default, NULL, uses \code{\link{getOption}}(digits).}	
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{
\tabular{ll}{
\tab -Family used.\cr
\tab -Number or iteration of algorithm and if it has converged. \cr
\tab -Residual and null deviance.\cr
\tab -Number of data.\cr
}
Produces a list of summary information for a fitted fregre.np object for each functional covariate. 
\tabular{ll}{
\tab -Call.\cr
\tab -R squared.\cr
\tab -Residual variance.\cr
\tab -Index of possible atypical curves or possible outliers.\cr
\tab -Index of possible influence curves.\cr
}     
If draw=TRUE plot: \cr
\tabular{ll}{
\tab -y vs y fitted values.\cr
\tab -Residuals vs fitted values.\cr
\tab -Residual boxplot.\cr
\tab -Quantile-Quantile Plot (qqnorm).\cr
\tab -Plot for a each single model term.\cr 
}
If \code{ask}=FALSE draw graphs in one window, by default. If \code{ask}=TRUE, draw each graph in a window, waiting to confirm.
}
\value{
\item{object}{ Object.}
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
Summary function for  \code{\link{fregre.kgam}}.
}
\examples{
data(tecator)
ab=tecator$absorp.fdata
ab2=fdata.deriv(ab,2)
yfat=as.integer(cut(tecator$y[,"Fat"],c(0,15,100)))-1
# Ajuste simple
xlist=list("df"=data.frame(yfat),"ab2"=ab2,"ab"=ab)
f<-yfat~ab+ab2
res=fregre.kgam(f,data=xlist,family=binomial("logit"),control=list(maxit=2))
summary(res)
res
}
\keyword{print}