\name{fregre.pls}
\Rdversion{1.1}
\alias{fregre.pls}
\title{ Functional PLS regression with scalar response}
\description{ Compute partial least squares regression (PLSR) for functional data. }
\usage{fregre.pls(fdataobj,y,l=1:3,...)}
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
  \item{l}{ Index of components to include in the model.}  
	\item{y}{ Scalar response with length \code{n}.}
  \item{\dots}{ Further arguments passed to or from other methods.} 
}
\details{
The partial least squares  are calculated by \code{\link{plsr}} function.
}
\value{
Return:
	\item{call}{ The matched call of \code{\link{fregre.pls}} function.}
	\item{beta.est}{ Beta coefficient estimated of class \code{fdata}.}
	\item{coefficients}{ A named vector of coefficients.}
	\item{fitted.values}{  Estimated scalar response.}
	%\item{pf}{ Formula used.}
	\item{coefficents}{ A named vector of coefficients.}	
	\item{residuals}{ \code{y}-\code{fitted values}.}	
	\item{df}{ The residual degrees of freedom.}	
	\item{r2}{ Coefficient of determination.}
	\item{sr2}{ Residual variance.}
	\item{H}{ Hat matrix.}	
	\item{fdataobj}{ Functional explanatory data.}
	\item{y}{ Scalar response.}
  \item{l}{ Index of components to include in the model.}  	
	%\item{Z}{ Matrix of centered \code{fdata} multiplied by the factors .}
	\item{fdata2pls}{ Fitted object in \code{\link{fdata2pls}} function.}
	\item{lm}{ Fitted object in \code{\link{lm}} function}	
	%\item{Z}{ Design matrix}   
}
\references{
Preda C. and Saporta G. \emph{PLS regression on a stochastic process}. Comput. Statist. Data Anal. 48 (2005): 149{-}158.
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