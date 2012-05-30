\name{fregre.pls}
\Rdversion{1.1}
\alias{fregre.pls}
\title{ Functional PLS regression with scalar response}
\description{ 
Computes functional linear regression between functional explanatory variable \eqn{X(t)} and scalar response \eqn{Y} using partial least squares (PLS). 
\deqn{Y=\big<X,\beta\big>+\epsilon=\int_{T}{X(t)\beta(t)dt+\epsilon}}{Y=<X,\beta>+\epsilon}

where \eqn{ \big<  \cdot , \cdot \big>}{<.,.>} denotes the inner product on \eqn{L_2} and \eqn{\epsilon} are random errors with mean zero , finite variance \eqn{\sigma^2} and \eqn{E[X(t)\epsilon]=0}{E[X(t)\epsilon]=0}.\cr

The functional partial least squared (FPLS) are calculated in \code{\link{fdata2pls}}  by alternative formulation of the NIPALS algorithm proposed by Kraemer and Sugiyama (2011).  
}
\usage{fregre.pls(fdataobj, y=NULL, l = NULL,...)}
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
	\item{y}{ Scalar response with length \code{n}.}
  \item{l}{ Index of components to include in the model.}  
  \item{\dots}{ Further arguments passed to or from other methods.} 
}
\details{
Functional partial least squares (FPLS) algorithm maximizes the covariance between \eqn{X(t)} and the scalar response \eqn{Y} via the partial least squares (PLS) components. Let \eqn{\left\{\tilde{\nu}_k\right\}_{k=1}^{\infty}}{{\nu_k}_k=1:\infty} the functional PLS components and \eqn{X_i(t)=\sum_{k=1}^{\infty}\tilde{\gamma}_{ik}\tilde{\nu}_k}{X_i(t)=\sum{k=1:\infty} \gamma_{ik} \nu_k} and  \eqn{\beta(t)=\sum_{k=1}^{\infty}\tilde{\beta}_k\tilde{\nu}_k}{\beta(t)=\sum{k=1:\infty} \beta_k \nu_k}. The functional linear model is estimated by:
\deqn{ \hat{y}=\big<  X,\hat{\beta} \big> \approx \sum_{k=1}^{k_n}\tilde{\gamma}_{k}\tilde{\beta}_k  }{ y.est=<  X,\beta.est > \approx \sum{k=1:k_n} \gamma_k \beta_k }
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
\keyword{regression}
