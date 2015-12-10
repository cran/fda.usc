\name{order.fdata}
\Rdversion{1.1}
\alias{order.fdata}
\title{ A wrapper for the  \code{\link{order}} function}
\description{  The funcional data is ordered w.r.t the sample order of the values of vector.}
\usage{
\method{order}{fdata}(y,fdataobj,na.last = TRUE, decreasing = FALSE)
}
\arguments{
  \item{y}{ a sequence of numeric, complex, character or logical vectors, all of the same length, or a classed R object. }
  \item{fdataobj}{\code{fdata} class object.}
  \item{na.last}{ for controlling the treatment of NAs. If TRUE, missing values in the data are put last; if FALSE, they are put first; if NA, they are removed; if "keep" they are kept with rank NA.}   
  \item{decreasing}{logical. Should the sort order be increasing or decreasing?.}
}
\value{ 
 \code{\link{order.fdata}} returns the functional data \code{fdataobj} w.r.t. a permutation which rearranges its first argument into ascending or descending order. 
}
%\references{
%Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
%Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
%}
\keyword{math}  