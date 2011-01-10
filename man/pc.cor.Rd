\name{pc.cor}
\Rdversion{1.1}
\alias{pc.cor}
\title{ Correlation for functional data by Principal Component Analysis}
\description{
Compute  correlation principal components of functional data and scalar response \code{y}.
}
\usage{
pc.cor(out,y=NULL,l=NULL,draw=TRUE,...)
}

\arguments{
  \item{out}{If \code{out} is  \code{fregre.pc} object is calculated  correlations between principal components and scalar response \code{y}.  If  \code{out} are functional data  the principal components  are calculated by \code{pc.svd.fdata} function and after are calculated the correlations between the principal  components and scalar response.}
  \item{y}{ Scalar Response }  
  \item{l}{ Vector of index of principal comoponents.}
	\item{draw}{ =TRUE draw correlations.} 	
	\item{\dots}{ Further arguments passed to or from other methods.}
}
\value{
    If \code{draw}=TRUE are drawn the correlations between the principal components and response \code{y}.\cr
If \code{ask}=TRUE, draw each graph in a window, waiting to confirm the change of page with a click of the mouse or pressing ENTER. 
If \code{ask}=FALSE draw graphs in one window.	 
    \item{cor}{ Return correlation between each principal component selected by \code{l} and scalar response \code{y}.} 
}
\references{
Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988). \emph{The New S Language}. Wadsworth \& Brooks/Cole. \cr

Venables, W. N. and B. D. Ripley (2002). \emph{Modern Applied Statistics with S}. Springer-Verlag. 
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as \code{\link{pc.fdata}}  and  \link[stats]{cor}
}
\examples{

library(e1071)
n=500
t=seq(0,1,len=101)
beta=t*sin(2*pi*t)^2 # beta con 2 colinas
x=matrix(NA,ncol=101,nrow=n)
y=numeric(n)
for(i in 1:n){
	x[i, ]=rwiener(1,101)*3+rnorm(101,sd=0.1)
	y[i]=mean(x[i,]*beta)+ rnorm(1,sd=0.1)
}
x<-fdata(x)
# Correlation
pc.cor(x,y,c(1,4,5.8),ask=TRUE)
}
\keyword{multivariate} 