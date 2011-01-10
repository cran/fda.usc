\name{semimetric.basis}
\Rdversion{1.1}
\alias{semimetric.basis}
\title{Proximities between  functional data}
\description{
Aproximates semi-metric distances for functional data of class \code{fdata} or \code{fd}.
}
\usage{
semimetric.basis(DATA1, DATA2 = DATA1,nderiv=0,type.basis1=NULL,
nbasis1=NULL,type.basis2=NULL,nbasis2=NULL,...)
}
\arguments{  
  \item{DATA1}{Functional data 1 or curve 1.}
  \item{DATA2}{Functional data 2 or curve 2.}
  \item{nderiv}{  Order of derivation, used in \code{deriv.fd}}
  \item{type.basis1}{ Type of Basis for \code{DATA1}.}  
  \item{nbasis1}{ Number of Basis for \code{DATA1}.}  
  \item{type.basis2}{ Type of Basis for \code{DATA2}.}  
  \item{nbasis2}{ Number of Basis for \code{DATA2.}}    
  \item{\dots}{ Further arguments passed to or from other methods.}    
}
\details{
Aproximates semi-metric distances for functional data of two \code{fd} class objects. 
If functional data are not functional \code{fd} class, the \code{semimetric.basis} function creates a basis to represent the functional data, by default is used \code{\link{create.bspline.basis}} and the \code{fdata} class object  is converted to \code{fd} class  using the \code{\link{Data2fd}}.\cr 
The function  calculates distances between the derivative of order \code{nderiv} of curves using \code{\link{deriv.fd}} function.
}
\value{
 Returns a proximities matrix between functional data.
}
\references{
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional data analysis.}
Springer Series in Statistics, New York.
}
\seealso{
See also  \code{\link{metric.lp}}, \code{\link{semimetric.NPFDA}} and \code{\link{deriv.fd}}  
}
\examples{
data(phoneme)
DATA1<-phoneme$learn
DATA2<-phoneme$test
a1=semimetric.basis(DATA1,DATA2)
a2=semimetric.basis(DATA1,DATA2,type.basis1="fourier",
nbasis1=11,type.basis2="fourier",nbasis2=11)
fd1 <- fdata2fd(DATA1)
fd2 <- fdata2fd(DATA2)
a3=semimetric.basis(fd1,fd2)
}
\keyword{ cluster }
