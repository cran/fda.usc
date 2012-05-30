\name{PCvM.statistic}
\alias{PCvM.statistic}
\alias{Aijr}

\title{PCvM statistic for the Functional Linear Model with scalar response}
\description{
  Projected Cramer-von Mises statistic (PCvM) for the Functional Linear Model with scalar response (FLM), \eqn{Y=\big<X,\beta\big>+\varepsilon,}{Y=<X,\beta>+\epsilon.}
 }

\usage{
PCvM.statistic (X, residuals, p, Aijr0)
Aijr (X, inpr)
}

\arguments{
  \item{X}{Functional covariate for the FLM. The object must be either in the class \code{\link{fdata}} or in the class \code{\link{fd}}.
  It is used to compute the matrix of inner products.}
  \item{residuals}{Residuals of the estimated FLM.}
  \item{p}{Number of elements of the functional basis where the functional covariate is represented.}
  \item{Aijr0}{Output from the \code{Aijr} function (see Details). Computed if not given.}
  \item{inpr}{Matrix of inner products of \code{X}. Computed if not given.}
 }

\details{
In order to optimize the computation of the statistic, the critical parts of these two functions are programmed in FORTRAN. The hardest part corresponds to the function \code{Aijr}, which involves the computation of an \code{\link{array}} of dimension \code{c(n,n,n)}.
As this array is symmetric for the first pair of dimensions, the order of the method can be reduced from \eqn{O(n^3)}{O(n^3)} to \eqn{O\big(\frac{n^3+n^2}{2}\big)}{O((n^3+n^2)/2)}. The value of \code{Aijr} is a vector of length \eqn{\frac{n^3+n^2}{2}}{(n^3+n^2)/2} 
containing the unique values of the original three dimensional array.
}

\value{
 For \code{PCvM.statistic}, the value of the statistic. For \code{Aijr}, a suitable output to be used in the argument \code{Aijr0}.
}

\references{

Escanciano, J. C. (2006). A consistent diagnostic test for regression models using projections. Econometric
Theory, 22, 1030-1051.

Garcia-Portugues, E., Gonzalez-Manteiga, W. and Febrero-Bande, M. (2012). A goodness--of--fit test for the functional linear model with scalar
response. \url{http://arxiv.org/abs/1205.6167}

}

\note{No NA's are allowed in the functional covariate.}

\author{Eduardo Garcia-Portugues. Please, report bugs and suggestions to \email{eduardo.garcia@usc.es}}

\seealso{\code{\link{flm.test}}}

\examples{

# Functional process
X=rproc2fdata(n=10,t=seq(0,1,l=101))

# Aijr
Aijr(X)

# Statistic
PCvM.statistic(X,residuals=rnorm(10),p=5)

}

\keyword{htest}

\keyword{flm}
\keyword{test}
\keyword{functional}
\keyword{linear}
\keyword{model}

