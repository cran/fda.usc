\name{fregre.pc}
\Rdversion{1.1}
\alias{fregre.pc}
\title{ Functional (Ridge) Regression with scalar response using Principal Components Analysis.}
\description{
 	Computes functional (ridge) regression between functional explanatory variables and scalar response using Principal Components Analysis.
}
\usage{
fregre.pc(fdataobj, y, l = NULL,rn = 0, ...) 
}
\arguments{
   \item{fdataobj}{ \code{\link{fdata}} class object or \code{fdata.comp} class object created by \code{\link{create.pc.basis}} function.}  
	\item{y}{ Scalar response with length \code{n}.}
	\item{l}{ Vector of principal components. If is null \code{l} (by default), \code{l=1:3}.}
	\item{rn}{Ridge parameter.  Default value is rn=0, i.e. no penalization is used.} 
	%If rn = FALSE (by default), rn=0.  If =TRUE, rn is a percentage of a larger singular values of centered functional data.} 
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{
The principal components are calculated by \code{\link{svd}} decomposition in the \code{\link{fdata2pc}} or \code{\link{fdata2pc}}functions.  \code{rn} controls the amount of shrinkage.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%se puede pasar una base en lugar del dato fucnional
}
\value{
Return:\cr
\item{call}{ The matched call of \code{\link{fregre.pc}} function.}
\item{beta.est}{ beta coefficient estimated of class \code{fdata}}
%\item{a.est}{ Intercept parameter estimated}
\item{coefficients}{ A named vector of coefficients.}
\item{fitted.values}{  Estimated scalar response.}
\item{residuals}{ \code{y}-\code{fitted values}.}
\item{H}{ Hat matrix.}
\item{df}{ The residual degrees of freedom. In ridge regression, \code{df(rn)} is the effective degrees of freedom.}
\item{r2}{ Coefficient of determination.}
\item{GCV}{ GCV criterion.}
\item{sr2}{ Residual variance.}
\item{l}{ Index of principal components selected.}
\item{rn}{ Amount of shrinkage.}
\item{fdata.comp}{ Fitted object in \code{\link{fdata2pc}} function.}
\item{lm}{ \code{lm} object.}
\item{fdataobj}{ Functional explanatory data.}
\item{y}{ Scalar response.}
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
y=tecator$y$Fat[ind]

res=fregre.pc(x,y)
summary(res)
res2=fregre.pc(x,y,l=c(1,3,4))
summary(res2)
# Functional Ridge Regression
res3=fregre.pc(x,y,l=c(1,3,4),rn=1)
summary(res3)
res4=fregre.pc(x,y,l=c(1,3,4),rn=10)
summary(res4)
betas<-c(res$beta.est,res2$beta.est,res3$beta.est,res4$beta.est)
plot(betas)

}
\keyword{regression}