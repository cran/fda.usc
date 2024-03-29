% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fdata.deriv.R
\name{fdata.deriv}
\alias{fdata.deriv}
\title{Computes the derivative of functional data object.}
\usage{
fdata.deriv(
  fdataobj,
  nderiv = 1,
  method = "bspline",
  class.out = "fdata",
  nbasis = NULL,
  ...
)
}
\arguments{
\item{fdataobj}{\code{\link{fdata}} class object.}

\item{nderiv}{Order of derivation, by defalult \code{nderiv}=1.}

\item{method}{Type of derivative method, for more information see
\bold{details}.}

\item{class.out}{Class of functional data returned: \code{fdata} or
\code{fd} class.}

\item{nbasis}{Number of Basis for \code{fdatataobj\$DATA}. It is only used
if method =\emph{"bspline"}, \emph{"exponential"}, \emph{"fourier"},
\emph{"monomial"} or \emph{"polynomial"}}

\item{\dots}{Further arguments passed to or from other methods.}
}
\value{
Returns the derivative of functional data of \code{fd} class if
\code{class.out}="\emph{fd}" or  \code{fdata} class if \code{class.out}="\emph{fdata}".
}
\description{
Computes the derivative of functional data.
\itemize{ 
\item If method =\emph{"bspline"}, \emph{"exponential"}, \emph{"fourier"}, 
\emph{"monomial"} or \emph{"polynomial"}.
 \code{fdata.deriv} function creates a basis to represent the functional
data. 
The functional data are converted to class \code{fd} using the
\code{\link{Data2fd}} function and the basis indicated in the \code{method}.
Finally, the function calculates the derivative of order \code{nderiv} of
curves using \code{\link{deriv.fd}} function.\cr 
\item If \code{method}=\emph{"fmm"}, \emph{"periodic"}, \emph{"natural"} or
\emph{"monoH.FC"} is used \code{\link{splinefun}} function.
\item If \code{method}=\emph{"diff"}, raw derivation is applied.  Not recommended to
use this method when the values are not equally spaced.\cr 
}
}
\examples{
data(tecator)
absorp=tecator$absorp.fdata
tecator.fd1=fdata2fd(absorp)
tecator.fd2=fdata2fd(absorp,"fourier",9)
tecator.fd3=fdata2fd(absorp,"fourier",nbasis=9,nderiv=1)
#tecator.fd1;tecator.fd2;tecator.fd3
tecator.fdata1=fdata(tecator.fd1)
tecator.fdata2=fdata(tecator.fd2)
tecator.fdata3=fdata(tecator.fd3)
tecator.fdata4=fdata.deriv(absorp,nderiv=1,method="bspline",
class.out='fdata',nbasis=9)
tecator.fd4=fdata.deriv(tecator.fd3,nderiv=0,class.out='fd',nbasis=9)
plot(tecator.fdata4)
plot(fdata.deriv(absorp,nderiv=1,method="bspline",class.out='fd',nbasis=11))

}
\seealso{
See also \code{\link{deriv.fd}}, \code{\link{splinefun}} and
\code{\link{fdata}}
}
\keyword{manip}
