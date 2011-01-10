\name{fregre.pc}
\Rdversion{1.1}
\alias{fregre.pc}
\title{ Functional Regression with scalar response using Principal Components Analysis.}
\description{
 	Computes functional regression between functional explanatory variables and scalar response using Principal Components Analysis.
}
\usage{fregre.pc(fdataobj,y,l=1:3,norm=TRUE)
}
\arguments{
   \item{fdataobj}{ \code{\link{fdata}} class object.}  
	\item{y}{ Scalar response with length \code{n}.}
	\item{l}{ Vector of principal comoponents.}
	\item{norm}{=TRUE the norm of eigenvectors \code{(rotation)}  is 1.} 
  %\item{\dots}{ Further arguments passed to or from other methods.}
}
\details{
The principal components are calculated by \code{\link{svd}} decomposition in the \code{\link{pc.svd.fdata}} function.
}
\value{
Return:\cr
\item{call}{ The matched call of \code{\link{fregre.pc}} function.}
\item{beta.est}{ beta coefficient estimated of class \code{fdata}}
%\item{a.est}{ Intercept parameter estimated}
\item{coefficients}{ A named vector of coefficients.}
\item{fitted.values}{  Estimated scalar response.}
\item{lambda}{ Eigenvalues of the covariance/correlation matrix.}
\item{vs}{ A matrix whose columns contain the eigenvectors selected by \code{l} index.}
%\item{vs}{ Eigenvectors for each  principal components -PC-.}
%\item{S}{ Smoothing matrix.}
\item{Z}{ Matrix of centered \code{fdata} multiplied by the eigenvalues \code{vs}.}
\item{H}{ Hat matrix.}
\item{df}{ The residual degrees of freedom.}
\item{residuals}{ \code{y}-\code{fitted values}.}
\item{r2}{ Coefficient of determination.}
\item{sr2}{ Residual variance.}
\item{fdataobj}{ Functional explanatory data.}
\item{y}{ Scalar response.}
\item{pf}{ Formula used in \code{lm} function}
%\item{kn}{ Number of principal components.}
\item{pc}{ \code{pc.svd.fdata} object.}
\item{l}{ Index of principal components selected.}
\item{lm}{ \code{lm} object.}
}



\references{

Cai TT, Hall P. 2006. \emph{Prediction in functional linear regression}. Annals of Statistics 34: 2159{-}2179.\cr 

Cardot H, Ferraty F, Sarda P. 1999. \emph{Functional linear model}. Statistics and Probability Letters 45: 11{-}22.\cr 

Hall P, Hosseini{-}Nasab M. 2006. \emph{On properties of functional principal components analysis}. Journal of the Royal Statistical
Society B 68: 109{-}126.
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:   \code{\link{fregre.pc.cv}}, \code{\link{summary.fregre.fd}} and \code{\link{predict.fregre.fd}}.\cr
Alternative method: \code{\link{fregre.basis}} and   \code{\link{fregre.np}}.
}
\examples{

data(tecator)
absorp=tecator$absorp.fdata
ind=1:129
x=absorp[ind,]
y=tecator$y[ind,2]

res=fregre.pc(x,y)
summary(res)
pc.cor(x,y,1:8)
pc.fdata(x,1:8,draw=FALSE)
res2=fregre.pc(x,y,l=c(1,3,4))
summary(res2)
x.d=fdata.deriv(x)
pc.cor(x.d,y,1:8)
pc.fdata(x.d,1:8,draw=FALSE)
res3=fregre.pc(x.d,y,l=c(1,3,4))
summary(res3)


}
\keyword{regression}