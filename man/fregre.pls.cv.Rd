\name{fregre.pls.cv}
\Rdversion{1.1}
\alias{fregre.pls.cv}
\title{ Functional PLS regression with scalar response using selection of number of  PLS components}
\description{																						
Functional Regression with scalar response using selection of number of partial leas squares (PLS) comoponents through cross-validation. The algoritm selects the  PLS components with best estimates the response. The	select is performed by cross-validation (CV) or  Model Selecttion Criteria (MSC).  After is computing  functional regression using the best selection of PLS components.
	}
\usage{fregre.pls.cv(fdataobj, y, kmax=8, criteria = "CV",...)}
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
  \item{kmax}{The number of PLS components to include in the model.}
  \item{y}{ Scalar response with length \code{n}.}
 \item{criteria}{ Type of cross-validation (CV) or  Model Selecttion Criteria (MSC) applied. Possible values are \emph{"CV"}, \emph{"AIC"}, \emph{"AICc"}, \emph{"SIC"}.} 	
 \item{\dots}{ Further arguments passed to \code{\link{fregre.pls}}}  
}
\details{
The partial least squares  are calculated by \code{\link{plsr}} function.\cr
The method selects the PLS components with minimum MSC criteria by stepwise regression.
}    
\value{
Return:\cr
\item{pls.opt}{ Index of PLS components selected.}
\item{MSC.min}{ Minimum Model Selecttion Criteria (MSC) value for the code{pls.opt} components.}
\item{MSC}{ Minimum  Model Selecttion Criteria (MSC) value for \code{kmax} components.}
\item{fregre.pls}{ Fitted regression object by the best \code{pls.opt} components.}
}
\author{ Febrero-Bande, M. and Oviedo de la Fuente, M. }
\seealso{
See Also as:   \code{\link{fregre.pls}}, \code{\link{summary.fregre.fd}} and \code{\link{predict.fregre.fd}}.\cr
Alternative method: \code{\link{fregre.pc}}, \code{\link{fregre.basis}} and   \code{\link{fregre.np}}.
}
\references{
Preda C. and Saporta G. \emph{PLS regression on a stochastic process}. Comput. Statist. Data Anal. 48 (2005): 149{-}158.
}
\examples{
data(tecator)
x<-tecator$absorp.fdata
y<-tecator$y$Fat
res1=fdata2pls(x,y,5)
norm.fdata(res1$rotation)
res2=fregre.pls(x,y,c(1:8))
summary(res2)
}
\keyword{multivariate} 