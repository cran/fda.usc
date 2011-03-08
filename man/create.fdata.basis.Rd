\name{create.fdata.basis}
\Rdversion{1.1}
\alias{create.fdata.basis}
\alias{create.pc.basis}
\title{  Create Basis Set for Functional Data of fdata class}
\description{Compute  basis for functional data.}
\usage{
create.fdata.basis(fdataobj,l=1:5,maxl=max(l),type.basis="bspline",
rangeval=fdataobj$rangeval,class.out="fd")
create.pc.basis(fdataobj,l=1:5)
}
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
  \item{l}{ Vector of basis index.}
  \item{maxl}{maximum number of basis}
  \item{type.basis}{ Type of basis (see create.basis function).}   
  \item{rangeval}{ A vector of length 2 giving the lower and upper limits of the range of permissible values for the function argument. }
  \item{class.out}{=="fd" basisfd class, =="fdata" fdata class.}
}
%\details{ }
\value{
    \item{basis}{ basis}
    \item{x}{ Is true the value of the rotated data (the centred data multiplied by the rotation matrix) is returned}     
    \item{mean}{ functional mean of \code{fdataobj}}
    \item{type}{ type of basis}  
}
\references{
Ramsay, James O. and Silverman, Bernard W. (2006), \emph{Functional Data Analysis}, 2nd ed., Springer, New York. \cr
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as  \link[fda]{create.basis} and   \code{\link{pc.fdata}}.
}
\examples{
data(tecator)
basis.pc<-create.pc.basis(tecator$absorp.fdata,c(1,4,5))
plot(basis.pc$basis)

basis.fd<-create.fdata.basis(tecator$absorp.fdata,c(1,4,5),
type.basis="fourier")
plot(basis.pc$basis)
basis.fdata<-create.fdata.basis(tecator$absorp.fdata,c(1,4,5),
type.basis="fourier",class.out="fdata")
plot(basis.fd,col=2,lty=1)
lines(basis.fdata,col=3,lty=1)
}
\keyword{multivariate} 