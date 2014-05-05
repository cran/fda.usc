\name{create.fdata.basis}
\Rdversion{1.1}
\alias{create.fdata.basis}
\alias{create.pc.basis}
\alias{create.pls.basis}
\alias{create.raw.fdata}
\title{  Create Basis Set for Functional Data of fdata class}
\description{Compute  basis for functional data.}
\usage{
create.fdata.basis(fdataobj,l=1:5,maxl=max(l),type.basis="bspline",
rangeval=fdataobj$rangeval,class.out="fd")
create.pc.basis(fdataobj,l=1:5,norm=TRUE,basis=NULL,lambda=0,P=c(0,0,1),...)
create.pls.basis(fdataobj, y, l = 1:5,norm=TRUE,lambda=0,P=c(0,0,1),...)
create.raw.fdata(fdataobj, l = 1:ncol(fdataobj))
}
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.} 
  \item{y}{ Vector of response (scalar).} 
  \item{l}{ Vector of basis index.}
  \item{maxl}{maximum number of basis}
  \item{type.basis}{ Type of basis (see create.basis function).}   
  \item{rangeval}{ A vector of length 2 giving the lower and upper limits of the range of permissible values for the function argument.}
  \item{norm}{=TRUE the norm of eigenvectors \code{basis} is 1.} 
  \item{class.out}{=="fd" basisfd class, =="fdata" fdata class.}
	\item{basis}{ "fd" basis object.}    
  \item{lambda}{ Amount of penalization. Default value is 0, i.e. no penalization is used.}
  \item{P}{ If P is a vector:  coefficients to define the penalty matrix object. By default P=c(0,0,1) penalize the second derivative (curvature) or acceleration.
  If P is a matrix:  the penalty matrix object.} 
   
 \item{\dots}{ Further arguments passed to or from other methods.}
}
%\details{ }
\value{
    \item{basis}{ basis}
    \item{x}{ Is true the value of the rotated data (the centred data multiplied by the rotation matrix) is returned}     
    \item{mean}{ functional mean of \code{fdataobj}}
    \item{df}{degree of freedom}    
    \item{type}{ type of basis}  
}
\references{
Ramsay, James O. and Silverman, Bernard W. (2006), \emph{Functional Data Analysis}, 2nd ed., Springer, New York. 

N. Kraemer, A.-L. Boulsteix, and G. Tutz (2008). Penalized Partial Least Squares with Applications to B-Spline Transformations and Functional Data. Chemometrics and Intelligent Laboratory
Systems, 94, 60 - 69. \url{http://dx.doi.org/10.1016/j.chemolab.2008.06.009}

}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{
See Also as  \link[fda]{create.basis} and   \code{\link{fdata2pc}}.
}
\examples{
data(tecator)
basis.pc<-create.pc.basis(tecator$absorp.fdata,c(1,4,5))
plot(basis.pc$basis,col=1)
basis.pls<-create.pls.basis(tecator$absorp.fdata,y=tecator$y[,1],c(1,4,5))
lines(basis.pls$basis,col=2)

basis.fd<-create.fdata.basis(tecator$absorp.fdata,c(1,4,5),
type.basis="fourier")
plot(basis.pc$basis)
basis.fdata<-create.fdata.basis(tecator$absorp.fdata,c(1,4,5),
type.basis="fourier",class.out="fdata")
plot(basis.fd,col=2,lty=1)
lines(basis.fdata,col=3,lty=1)
}
\keyword{multivariate} 