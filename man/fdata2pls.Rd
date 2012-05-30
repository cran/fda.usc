\name{fdata2pls}
\Rdversion{1.1}
\alias{fdata2pls}
\title{ Partial least squares components for functional data.}
\description{ Compute partial least squares (PLS) components for functional data. }
\usage{ 
fdata2pls(fdataobj,y,ncomp=2,...)
 }
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
	\item{y}{ Scalar response with length \code{n}.}
  \item{ncomp}{The number of components to include in the model.}
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{
 %SIMPLS algorithm (de Jong, 1993)
The partial least squares are calculated by \code{NIPALS} algorithm.
}
\value{
\code{fdata2pls} function return:
    \item{rotation}{ \code{\link{fdata}} class object.} 
    \item{x}{ Is true the value of the rotated data (the centred data multiplied by the rotation matrix) is returned.}      
%	  \item{res.pls}{ Fitted regression object. See \code{\link{plsr}} for details.}     
    \item{fdataobj.cen}{ The centered \code{fdataobj} object.}  
    \item{mean}{ mean of \code{fdataobj}.}         
 	  \item{l}{ Vector of index of principal components.}    
     \item{C}{ The matched call.}     
}
\references{ 
Kraemer, N., Sugiyama M. (2011). \emph{The Degrees of Freedom of Partial Least Squares Regression}. Journal of the American Statistical Association. Volume 106, 697-705.
}
\author{ Febrero-Bande, M. and Oviedo de la Fuente, M. }
\seealso{
%See Also as:  \code{\link{fdata2pc}}. \cr
Used in: \code{\link{fregre.pls}}, \code{\link{fregre.pls.cv}}.\cr
Alternative method: \code{\link{fdata2pc}}.
}
\examples{
n= 500;tt= seq(0,1,len=101)
x0<-rproc2fdata(n,tt,sigma="wiener")
x1<-rproc2fdata(n,tt,sigma=0.1)
x<-x0*3+x1
beta = tt*sin(2*pi*tt)^2
fbeta = fdata(beta,tt)
y<-inprod.fdata(x,fbeta)+rnorm(n,sd=0.1)
pls1=fdata2pls(x,y)
norm.fdata(pls1$rotation)
}
\keyword{multivariate} 