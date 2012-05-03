\name{fregre.pls}
\Rdversion{1.1}
\alias{fregre.pls}
\title{ Functional PLS regression with scalar response}
\description{ Compute partial least squares regression (PLSR) for functional data. }
\usage{fregre.pls(fdataobj, y=NULL, l = NULL,...)}

\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
	\item{y}{ Scalar response with length \code{n}.}
  \item{l}{ Index of components to include in the model.}  
  \item{\dots}{ Further arguments passed to or from other methods.} 
}
\details{
The principal components are calculated by \code{NIPALS} algorithm in the \code{\link{fdata2pls}} function.  
}
\value{
Return:
	\item{call}{ The matched call of \code{\link{fregre.pls}} function.}
	\item{beta.est}{ Beta coefficient estimated of class \code{fdata}.}
	\item{coefficients}{ A named vector of coefficients.}
	\item{fitted.values}{  Estimated scalar response.}
	\item{residuals}{ \code{y}-\code{fitted values}.}	
	\item{H}{ Hat matrix.}	
	\item{df}{ The residual degrees of freedom.}	
	\item{r2}{ Coefficient of determination.}
	\item{GCV}{ GCV criterion.}
	\item{sr2}{ Residual variance.}	
  \item{l}{ Index of components to include in the model.}  	
	\item{rn}{ Amount of shrinkage.}
	\item{fdata.comp}{ Fitted object in \code{\link{fdata2pls}} function.}
	\item{lm}{ Fitted object in \code{\link{lm}} function}	
	\item{fdataobj}{ Functional explanatory data.}
	\item{y}{ Scalar response.}
} 
\references{
Preda C. and Saporta G. \emph{PLS regression on a stochastic process}. Comput. Statist. Data Anal. 48 (2005): 149{-}158.\cr
Kraemer, N., Sugiyama M. (2011). \emph{The Degrees of Freedom of Partial Least Squares Regression}. Journal of the American Statistical Association. Volume 106, 697-705.
}
\author{ Febrero-Bande, M. and Oviedo de la Fuente, M. }
\seealso{
See Also as:   \code{\link{fregre.pls.cv}}, \code{\link{summary.fregre.fd}} and \code{\link{predict.fregre.fd}}.\cr
Alternative method: \code{\link{fregre.pc}}, \code{\link{fregre.basis}} and   \code{\link{fregre.np}}.
}
\examples{

data(tecator)
x<-tecator$absorp.fdata
y<-tecator$y$Fat
res=fregre.pls(x,y,c(1:8))
summary(res)
}

\keyword{multivariate} 