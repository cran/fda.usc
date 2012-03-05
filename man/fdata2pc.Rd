\name{fdata2pc}
\Rdversion{1.1}
\alias{fdata2pc}
\title{ Principal components for functional data}
\description{
Compute  principal components for functional data.
}
\usage{
fdata2pc(fdataobj,norm=TRUE,ncomp=2)
%pc.fdata(out,l=NULL,loadings=FALSE,draw=TRUE,...)
}
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
    \item{norm}{=TRUE the norm of eigenvectors \code{(rotation)}  is 1.} 
    \item{ncomp}{ Number of  principal comoponents.}
  %\item{loadings}{ =TRUE "rotate"' loading matrices by \link[stats]{varimax} function.}
  %\item{draw}{ =TRUE draw correlations.}   
  %\item{\dots}{ Further arguments passed to or from other methods.}
}
\value{
    \item{lambda}{ Eigenvalues of the principal components.} 
    \item{rotation}{ \code{\link{fdata}} class object.} 
    \item{x}{ Is true the value of the rotated data (the centred data multiplied by the rotation matrix) is returned}      
    \item{fdataobj.cen}{ The centered \code{fdataobj} object.}  
    \item{mean}{ The functional mean of \code{fdataobj} object.}      
 	  \item{l}{ Vector of index of principal comoponents.}    
     \item{C}{ The matched call.}   	  
}
\references{
Venables, W. N. and B. D. Ripley (2002). \emph{Modern Applied Statistics with S}. Springer-Verlag. 
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as  \link[base]{svd} and  \link[stats]{varimax}.
}
\examples{

n= 500;tt= seq(0,1,len=101)
x0<-rproc2fdata(n,tt,sigma="wiener")
x1<-rproc2fdata(n,tt,sigma=0.1)
x<-x0*3+x1
pc=fdata2pc(x)

}
\keyword{multivariate} 