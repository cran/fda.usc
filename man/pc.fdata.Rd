\name{pc.fdata}
\Rdversion{1.1}
\alias{pc.fdata}
\alias{pc.svd.fdata}
\title{ Principal components for functional data}
\description{
Compute  principal components for functional data.
}
\usage{
pc.svd.fdata(fdataobj,norm=TRUE)
pc.fdata(out,l=NULL,loadings=FALSE,draw=TRUE,...)
}
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
     \item{norm}{=TRUE the norm of eigenvectors \code{(rotation)}  is 1.} 
  \item{out}{If  \code{out}  is \code{fregre.fd} class see \code{\link{fregre.pc}} else if \code{out}  is a \code{fdata} class object the principal components  are calculated by \code{pc.svd.fdata} function.}
  \item{l}{ Vector of index of principal comoponents.}
  \item{loadings}{ =TRUE "rotate"' loading matrices by \link[stats]{varimax} function.}
  \item{draw}{ =TRUE draw correlations.}   
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\value{
    \item{lambda}{ Eigenvalues of the principal components.} 
    \item{rotation}{ \code{\link{fdata}} class object.} 
    \item{x}{ Is true the value of the rotated data (the centred data multiplied by the rotation matrix) is returned}      
    \item{fdataobj.cen}{ The centered \code{fdataobj} object.}  
    \item{mean}{ The functional mean of \code{fdataobj} object.}      
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

library(e1071)
n=500
tt=seq(0,1,len=101)
beta=tt*sin(2*pi*tt)^2
x=matrix(NA,ncol=101,nrow=n)
y=numeric(n)
for(i in 1:n){
	x[i, ]=rwiener(1,101)*3+rnorm(101,sd=0.1)
	y[i]=mean(x[i,]*beta)+ rnorm(1,sd=0.1)
}
x.fdata<-fdata(x,tt)
# by svd() 
pc1=pc.svd.fdata(x.fdata)

# by svd in functional regression
res=fregre.pc(x.fdata,y,l=c(1,4,5))
pc2=pc.fdata(res)
pc3=pc.fdata(res,ask=TRUE)
}
\keyword{multivariate} 