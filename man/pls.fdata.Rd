\name{pls.fdata}
\Rdversion{1.1}
\alias{pls.fdata}
\title{ Partial least squares components for functional data.}
\description{ Compute partial least squares (PLS) components for functional data. }
\usage{ pls.fdata(fdataobj,y,ncomp=2,norm=TRUE,...) }
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
	\item{y}{ Scalar response with length \code{n}.}
  \item{ncomp}{The number of components to include in the model.}
  \item{norm}{=TRUE the norm of eigenvectors \code{(rotation)}  is 1.} 
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{
The partial least squares  are calculated by \code{\link{plsr}} function.
}
\value{
\code{pls.fdata} function return:
    \item{rotation}{ \code{\link{fdata}} class object.} 
    \item{x}{ Is true the value of the rotated data (the centred data multiplied by the rotation matrix) is returned.}      
	  \item{res.pls}{ Fitted regression object. See \code{\link{plsr}} for details.}     
    \item{fdataobj.cen}{ The centered \code{fdataobj} object.}  
    \item{mean}{ mean of \code{fdataobj}.}      
}
\references{ 
Preda C. and Saporta G. \emph{PLS regression on a stochastic process}. Comput. Statist. Data Anal. 48 (2005): 149{-}158.
}
\author{ Febrero-Bande, M. and Oviedo de la Fuente, M. }
\seealso{
See Also as:  \code{\link{plsr}}. \cr
Used in: \code{\link{fregre.pls}}, \code{\link{fregre.pls.cv}}.\cr
Alternative method: \code{\link{pc.svd.fdata}}.
}
\examples{
data(tecator)
x<-tecator$absorp.fdata
y<-tecator$y$Fat
res1=pls.fdata(x,y,5)
norm.fdata(res1$rotation)
res2=fregre.pls(x,y,c(1:8))
summary(res2)
}
\keyword{multivariate} 