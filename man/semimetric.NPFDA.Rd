\name{semimetric.NPFDA}
\Rdversion{1.1}
\alias{semimetric.NPFDA}
\alias{semimetric.deriv}
\alias{semimetric.fourier}
\alias{semimetric.hshift}
\alias{semimetric.mplsr}
\alias{semimetric.pca}
\alias{mplsr}
\alias{hshift}
\title{Proximities between  functional data (semi-metrics)}
\description{
Aproximates semi-metric distances for functional data  by  Ferraty F and Vieu, P. (2006) aproach.
\code{semimetric.deriv}: computes a semi-metric between curves based on their derivatives.\cr
\code{semimetric.fourier}: computes a semi-metric between curves based on their Fourier expansion.\cr
\code{semimetric.hshift}: computes between curves a semi-metric taking into account an horizontal shift effect.\cr
\code{semimetric.mplsr}: computes between curves a semi-metric based on the partial least  squares method.\cr 
\code{semimetric.pca}: computes between curves a pca-type semi-metric based on the  functional principal components analysis method.
}
\usage{
semimetric.hshift(fdataobj1, fdataobj2, t=1:ncol(DATA1),...)
semimetric.mplsr(fdataobj1, fdataobj2, q=2, class1,...)
semimetric.pca(fdataobj1, fdataobj2, q=1,...)
semimetric.deriv(fdataobj1, fdataobj2, nderiv=1,
nknot=ifelse(floor(ncol(DATA1)/3)>floor((ncol(DATA1)-nderiv-4)/2),
floor((ncol(DATA1)-nderiv-4)/2),floor(ncol(DATA1)/3)),
range.t=c(0,1),...)
semimetric.fourier(fdataobj1, fdataobj2, nderiv=0,
nbasis=ifelse(floor(ncol(DATA1)/3)>floor((ncol(DATA1)-nderiv-4)/2),
floor((ncol(DATA1) - nderiv - 4)/2), floor(ncol(DATA1)/3)),
range.t=c(0,1), period=NULL,...)
}
\arguments{  
  \item{fdataobj1}{Functional data 1 or curve 1. \code{DATA1} with dimension (\code{n1} x \code{m}), where \code{n1} is the number of curves and \code{m} are the points observed in each curve.}
  \item{fdataobj2}{Functional data 2 or curve 2. \code{DATA1} with dimension (\code{n2} x \code{m}), where \code{n2} is the number of curves and \code{m} are the points observed in each curve.}
  \item{q}{  If \code{semimetric.pca}: the retained number of principal components.\cr
	If \code{semimetric.mplsr}: the retained number of factors.}
  \item{nknot}{  semimetric.deriv argument: number of interior knots (needed for defining the B-spline basis).}
  \item{range.t}{ semimetric.deriv argument: vector of length 2 containing the range of the t at which the curve are evaluated (i.e. range of the discretization).}
 \item{nderiv}{  Order of derivation, used in \code{semimetric.deriv} and \code{semimetric.fourier}}
\item{nbasis}{ \code{semimetric.fourier}: size of the basis.}
\item{period}{ \code{semimetric.fourier}:allows to select the period for the fourier expansion.} 
\item{t}{ \code{semimetric.hshift}: vector which defines  \code{t} (one can choose \code{1,2,...,nbt} where \code{nbt} is the number of points of the discretization)}
\item{class1}{ \code{semimetric.mplsr}: vector containing a categorical response which  corresponds to class number for units stored in \code{DATA1}.} 
\item{\dots}{ Further arguments passed to or from other methods.}

}
\details{
The following instructions have been used file: \url{http://www.math.univ-toulouse.fr/staph/npfda/}\cr
}
\value{
 Returns a \code{semimetric.NPFDA} matrix containing the semimetric computed  between the curves lying to the first sample and the curves lying
 to the second one.
}
\source{
\url{http://www.math.univ-toulouse.fr/staph/npfda/}
}
\references{
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional data analysis.}
Springer Series in Statistics, New York.\cr
Ferraty, F. and Vieu, P. (2006). \emph{NPFDA in practice}. \cr
Free access on line at \url{http://www.lsp.ups-tlse.fr/staph/npfda/}
}
\seealso{
See also  \code{\link{metric.lp}} and  \code{\link{semimetric.basis}}
}
\examples{

#	INFERENCE PHONDAT
data(phoneme)
#ind=c(1:50,201:250)
ind<-1:250
mlearn<-phoneme$learn[ind,]
mtest<-phoneme$test[ind,]
n=nrow(mlearn[["data"]])
np=ncol(mlearn[["data"]])
mdist1=semimetric.pca(mlearn,mtest)
mdist2=semimetric.pca(mlearn,mtest,q=2)
mdist3=semimetric.deriv(mlearn,mtest,nderiv=0)
mdist4=semimetric.fourier(mlearn,mtest,nderiv=2,nbasis=21,c(1,np))
#uses hshift function
#mdist5=semimetric.hshift(mlearn,mtest) #takes a lot
#mdist6=semimetric.hshift(mlearn,mtest,1:np) #takes a lot
glearn<-phoneme$classlearn[ind]
#uses mplsr function
mdist7=semimetric.mplsr(mlearn,mtest,5,glearn)
mdist0=metric.lp(mlearn,mtest)

b=as.dist(mdist7)
c2=hclust(b)
plot(c2)
memb <- cutree(c2, k = 2)
table(memb,phoneme$classlearn[ind])
   

}
\keyword{ cluster }