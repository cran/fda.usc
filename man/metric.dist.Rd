\name{metric.dist}
\Rdversion{1.1}
\alias{metric.dist}
\title{
Distance Matrix Computation
}
\description{	
This function computes the distances between the rows of a data matrix  by using the specified distance measure.
}
\usage{
metric.dist(x, y = NULL, method = "euclidean", p = 2, ...) 
}
\arguments{  
  \item{x}{Data frame 1. The dimension is (\code{n1} x \code{m}).}
  \item{y}{Data frame 2. The dimension is (\code{n2} x \code{m}).}
  \item{method}{ The distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".}
  \item{p}{ The power of the Minkowski distance.}
  \item{\dots}{ Further arguments passed to \code{\link{dist}} function.}
}
\details{
This function returns a distance matrix by using \code{\link{dist}} function. \cr
The matrix dimension is (\code{n1} x \code{n1}) if \code{y=NULL},  (\code{n1} x \code{n2})  otherwise.
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{
See also  \code{\link{dist}} for multivariate date case and \code{\link{metric.lp} for functional data case}
} 
\examples{
# data(iris)
# d<-metric.dist(iris[,1:4])
# matplot(d,type="l",col=as.numeric(iris[,5]))
}
\keyword{ cluster }