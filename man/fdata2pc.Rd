\name{fdata2pc}
\Rdversion{1.1}
\alias{fdata2pc}
\title{ Principal components for functional data}
\description{
Compute  principal components for functional data.
}
\usage{
fdata2pc(fdataobj,  ncomp = 2,norm = TRUE,...)
}
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
  \item{ncomp}{ Number of  principal comoponents.}
  \item{norm}{=TRUE the norm of eigenvectors \code{(rotation)}  is 1.} 
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\value{
    \item{d}{  The standard deviations of the functional principal components.} 
    \item{rotation}{ are also known as loadings.  A \code{fdata} class object whose rows contain the eigenvectors).} 
    \item{x}{  are also known as scores. The value of the rotated functional data is returned. }      
    \item{fdataobj.cen}{ The centered \code{fdataobj} object.}  
    \item{mean}{ The functional mean of \code{fdataobj} object.}      
    \item{l}{ Vector of index of principal components.}    
     \item{C}{ The matched call.}   	  
}
\references{
Venables, W. N. and B. D. Ripley (2002). \emph{Modern Applied Statistics with S}. Springer-Verlag. 

Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{
See Also as  \link[base]{svd} and  \link[stats]{varimax}.
}
\examples{

n= 200;tt= seq(0,1,len=101)
x0<-rproc2fdata(n,tt,sigma="wiener")
x1<-rproc2fdata(n,tt,sigma=0.1)
x<-x0*3+x1
pc=fdata2pc(x)

}
\keyword{multivariate} 