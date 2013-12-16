\name{fregre.pls}
\Rdversion{1.1}
\alias{fregre.pls}
\title{ Functional Penalized  PLS regression with scalar response}
\description{ 
Computes functional linear regression between functional explanatory variable \eqn{X(t)} and scalar response \eqn{Y} using penalized Partial Least Squares (PLS)%,  where \eqn{\tilde{X(t)}=MX(t)} with  \eqn{M=(I+\lambda P)^{-1}}. 
\deqn{Y=\big<\tilde{X},\beta\big>+\epsilon=\int_{T}{\tilde{X}(t)\beta(t)dt+\epsilon}}{Y=<\tilde{X},\beta>+\epsilon}

where \eqn{ \big<  \cdot , \cdot \big>}{<.,.>} denotes the inner product on \eqn{L_2} and \eqn{\epsilon} are random errors with mean zero , finite variance \eqn{\sigma^2} and \eqn{E[\tilde{X}(t)\epsilon]=0}{E[X(t)\epsilon]=0}.\cr
}
\usage{
fregre.pls(fdataobj, y=NULL, l = NULL,lambda=0,P=c(0,0,1),...)
}
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
	\item{y}{ Scalar response with length \code{n}.}
  \item{l}{ Index of components to include in the model.}  
   \item{lambda}{ Amount of penalization. Default value is 0, i.e. no penalization is used.}
 \item{P}{ If \code{P} is a vector:  \code{P} are coefficients to define the penalty matrix object. By default \code{P=c(0,0,1)} penalize the second derivative (curvature) or acceleration.
  If \code{P} is a matrix: P is the penalty matrix object.} 
  \item{\dots}{ Further arguments passed to or from other methods.} 
} 
\details{

%The function computes the \eqn{\left\{\nu_k\right\}_{k=1}^{\infty}}{\nu_1,...,\nu_\infty} orthonormal basis of PLS to represent the functional data as \eqn{X_i(t)=\sum_{k=1}^{\infty}\gamma_{ik}\nu_k}{X(t)=\sum_(k=1:\infty) \gamma_k \nu_k}.%, where \eqn{\tilde{X}=MX} with  \eqn{M=(I+\lambda P)^{-1}},\eqn{\gamma_{ik}=\Big< \tilde{X}_i(t),\nu_k\Big>}{\gamma_k= < \tilde{X}_i(t), \nu_k >} . \cr

Functional (FPLS) algorithm maximizes the covariance between \eqn{X(t)} and the scalar response \eqn{Y} via the partial least squares (PLS) components. The functional penalized  PLS are calculated in  \code{\link{fdata2pls}}  by alternative formulation of the NIPALS algorithm proposed by Kraemer and Sugiyama (2011).   

Let \eqn{\left\{\tilde{\nu}_k\right\}_{k=1}^{\infty}}{{\nu_k}_k=1:\infty} the functional PLS components and \eqn{\tilde{X}_i(t)=\sum_{k=1}^{\infty}\tilde{\gamma}_{ik}\tilde{\nu}_k}{X_i(t)=\sum{k=1:\infty} \gamma_{ik} \nu_k} and  \eqn{\beta(t)=\sum_{k=1}^{\infty}\tilde{\beta}_k\tilde{\nu}_k}{\beta(t)=\sum{k=1:\infty} \beta_k \nu_k}. The functional linear model is estimated by:
\deqn{ \hat{y}=\big< X,\hat{\beta} \big> \approx \sum_{k=1}^{k_n}\tilde{\gamma}_{k}\tilde{\beta}_k  }{ y.est=<  X,\beta.est > \approx \sum{k=1:k_n} \gamma_k \beta_k }

The response can be fitted by:
\itemize{
\item \eqn{\lambda=0}{\lambda=0}, no penalization, \deqn{\hat{y}=\nu_k^{\top}(\nu_k^{\top}\nu_k)^{-1}\nu_k^{\top}y}{y.est= \nu'(\nu'\nu)^{-1}\nu'y}


\item Penalized regression, \eqn{\lambda>0}{\lambda>0} and \eqn{P\neq0}{P!=0}. For example, \eqn{P=c(0,0,1)}{P=c(0,0,1)} penalizes the second derivative (curvature) by \code{P=P.penalty(fdataobj["argvals"],P)}, \deqn{\hat{y}=\nu_k^{\top}(\nu_k\top \nu_k+\lambda \nu_k^{\top} \textbf{P}\nu_k)^{-1}\nu_k^{\top}y}{y.est=\nu'(\nu'\nu+\lambda v'Pv)^{-1}\nu'y}
}
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
	\item{lambda}{ Amount of shrinkage.}
	\item{fdata.comp}{ Fitted object in \code{\link{fdata2pls}} function.}
	\item{lm}{ Fitted object in \code{\link{lm}} function}	
	\item{fdataobj}{ Functional explanatory data.}
	\item{y}{ Scalar response.}
} 
\references{
Preda C. and Saporta G. \emph{PLS regression on a stochastic process}. Comput. Statist. Data Anal. 48 (2005): 149{-}158.

N. Kraemer, A.-L. Boulsteix, and G. Tutz (2008). \emph{Penalized Partial Least Squares with Applications
to B-Spline Transformations and Functional Data}. Chemometrics and Intelligent Laboratory Systems, 94, 60 - 69. \url{http://dx.doi.org/10.1016/j.chemolab.2008.06.009}

Martens, H., Naes, T. (1989) \emph{Multivariate calibration.} Chichester: Wiley.

Kraemer, N., Sugiyama M. (2011). \emph{The Degrees of Freedom of Partial Least Squares Regression}. Journal of the American Statistical Association. Volume 106, 697-705.

Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
}
\author{ Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es} }
\seealso{
See Also as:   \code{\link{P.penalty}} and \code{\link{fregre.pls.cv}}.\cr
Alternative method: \code{\link{fregre.pc}}.
}
\examples{

\dontrun{
data(tecator)
x<-tecator$absorp.fdata
y<-tecator$y$Fat
res=fregre.pc(x,y,c(1:8))
summary(res)
res2=fregre.pls(x,y,c(1:8),lambda=10)
summary(res2)
}
}
\keyword{regression}
