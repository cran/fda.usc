% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/int.simpson.R
\name{int.simpson}
\alias{int.simpson}
\alias{int.simpson2}
\title{Simpson integration}
\usage{
int.simpson(fdataobj, method = NULL)

int.simpson2(x, y, equi = TRUE, method = NULL)
}
\arguments{
\item{fdataobj}{fdata objtect.}

\item{method}{Method for numerical integration, see details.}

\item{x}{Sorted vector of x-axis values: \code{argvals}.}

\item{y}{Vector of y-axis values.}

\item{equi}{=TRUE, the observed points on each curve are equally spaced (by
default).}
}
\description{
Computes the integral of \code{fdataobj$data} with respect to
\code{fdataobj$argvals} using simpson or trapezoid rule integration.
}
\details{
Posible values for \code{method} are: \itemize{ \item \code{"TRAPZ"}:
Trapezoid rule integration. \item \code{"CSR"}: Composite Simpson's rule
integration.  \item \code{"ESR"}: Extended Simpson's rule integration. } If
\code{method=NULL} (default), the value of \code{par.fda.usc$int.method} is
used.
}
\examples{
\dontrun{
x<-seq(0,2*pi,length=1001)
fx<-fdata(sin(x)/sqrt(pi),x)
fx0<-fdata(rep(0,length(x)),x)
int.simpson(fx0)
int.simpson(fx)
}
}
\seealso{
See also \code{\link{integrate}}.
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente
\email{manuel.oviedo@udc.es}
}
\keyword{cluster}
