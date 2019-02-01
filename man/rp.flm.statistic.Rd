\name{rp.flm.statistic}
\alias{rp.flm.statistic}
%\alias{Adot}

\title{ Statistic for testing the FLM using random projections}
\description{
  Random projections statistic (PCvM) for the Functional Linear Model with scalar response (FLM): \eqn{Y=\big<X,\beta\big>+\varepsilon}{Y=<X,\beta>+\epsilon}.
 }

\usage{
%PCvM.statistic (X, residuals, p, Adot.vec)
rp.flm.statistic(proj.X,residuals,proj.X.ord=NULL,verbose=FALSE,F.code=TRUE)
}

\arguments{
  \item{proj.X}{Functional projections covariate for the FLM. The object must be either in the class \code{\link{fdata}} or in the class \code{\link{fd}}.
  It is used to compute the matrix of inner products.}
  \item{residuals}{Residuals of the estimated FLM.}
  \item{proj.X.ord}{Order projections.}
  \item{verbose}{Either to show or not information about computing progress.}
  \item{F.code}{ \code{logical}, Check if the fortran dll is loaded}  
}

%\details{
%In order to optimize the computation of the statistic, the critical parts of these two functions are programmed in FORTRAN. The hardest part corresponds to the function \code{Adot}, which involves the computation of a symmetric matrix of dimension \eqn{n\times n}{n x n} where each entry is a sum of \eqn{n}{n} elements.
%As this matrix is symmetric, the order of the method can be reduced from \eqn{O(n^3)}{O(n^3)} to \eqn{O\big(\frac{n^3-n^2}{2}\big)}{O((n^3-n^2)/2)}. The memory requirement can also be reduced to \eqn{O\big(\frac{n^2-n+2}{2}\big)}{O((n^2-n+2)/2)}. The value of \code{Adot} is a vector of length \eqn{\frac{n^2-n+2}{2}}{(n^2-n+2)/2} 
%where the first element is the common diagonal element and the rest are the lower triangle entries of the matrix, sorted by rows (see Examples).
%}
\value{
 Return a list containing the following components:  
 \item{rp.flm.statistic}{ The value of the statistic.}
 \item{proj.X.ord}{ Order projections.} 
}
\note{No NA's are allowed in the functional covariate.}

\author{Eduardo Garcia-Portugues. Please, report bugs and suggestions to \if{latex}{\cr}\email{eduardo.garcia@usc.es}}

\seealso{\code{\link{rp.flm.test}}}

\keyword{htest}
