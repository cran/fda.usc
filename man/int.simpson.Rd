\name{int.simpson}
\Rdversion{1.1}
\alias{int.simpson}
\alias{int.simpson2}
\title{ Simpson integration }
\description{	
Computes the integral of y with respect to x using simpson or trapezoid rule integration.
}
\usage{
int.simpson(fdataobj,equi=TRUE,method="TRAPZ")
int.simpson2(x, y, equi = TRUE, method = "TRAPZ") 
}
\arguments{  
 \item{fdataobj}{ fdata objtect.}
  \item{x}{Sorted vector of x-axis values: \code{argvals}.}
  \item{y}{Vector of y-axis values.}
  \item{equi}{ =TRUE, the observed points on each curve are equally spaced (by default). }
  \item{method}{ Method for numerical integration, see details.}   
}
\details{
If method="CSR", composite Simpson's rule integration is used. \cr
If method="ESR", extended Simpson's rule integration is used.\cr
If method="TRAPZ", Trapezoid rule integration is used.
If \code{x} is a fdata class object it computes the integral of \code{x$data} with Respect to \code{x$argvals}.
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See also  \code{\link{integrate}}.
} 
\examples{
x<-seq(0,2*pi,length=1001)
fx<-fdata(sin(x)/sqrt(pi),x)
fx0<-fdata(rep(0,length(x)),x)
int.simpson(fx0$data,fx$data)
int.simpson(fx)
int.simpson2(fx$argvals,fx$data^2)
int.simpson2(fx$argvals,fx$data^2,method="TRAPZ")
int.simpson2(fx$argvals,fx$data^2,method="ESR")

}
\keyword{ cluster }