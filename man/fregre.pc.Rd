\name{fregre.pc}
\Rdversion{1.1}
\alias{fregre.pc}
\title{ Functional (Ridge) Regression with scalar response using Principal Components Analysis.}
\description{
Computes functional (ridge) regression between functional explanatory variable \eqn{X(t)} and scalar response \eqn{Y} using Principal Components Analysis. 	

\deqn{Y=\big<X,\beta\big>+\epsilon=\int_{T}{X(t)\beta(t)dt+\epsilon}}{Y=<X,\beta>+\epsilon}

where \eqn{ \big<  \cdot , \cdot \big>}{<.,.>} denotes the inner product on \eqn{L_2} and \eqn{\epsilon} are random errors with mean zero , finite variance \eqn{\sigma^2} and \eqn{E[X(t)\epsilon]=0}{E[X(t)\epsilon]=0}.
 
}
\usage{
fregre.pc(fdataobj, y, l = NULL,rn = 0,weights= rep(1,len=n),...) 
}
\arguments{
   \item{fdataobj}{ \code{\link{fdata}} class object or \code{fdata.comp} class object created by \code{\link{create.pc.basis}} function.}  
	\item{y}{ Scalar response with length \code{n}.}
	\item{l}{ Vector of principal components. If is null \code{l} (by default), \code{l=1:3}.}
	\item{rn}{Ridge parameter.  Default value is rn=0, i.e. no penalization is used.} 
	\item{weights}{ weights}
	%If rn = FALSE (by default), rn=0.  If =TRUE, rn is a percentage of a larger singular values of centered functional data.} 
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{
The function computes the \eqn{\left\{\nu_k\right\}_{k=1}^{\infty}}{\nu_1,...,\nu_\infty} orthonormal basis of functional principal components to represent the functional data as \eqn{X_i(t)=\sum_{k=1}^{\infty}\gamma_{ik}\nu_k}{X(t)=\sum_(k=1:\infty) \gamma_k \nu_k} and the functional parameter as \eqn{\beta(t)=\sum_{k=1}^{\infty}\beta_k\nu_k}{\beta(t)=\sum_(k=1:\infty) \beta_k \nu_k}, where \eqn{\gamma_{ik}=\Big< X_i(t),\nu_k\Big>}{\gamma_k= < X(t), \nu_k >} and \eqn{\beta_{k}=\Big<\beta,\nu_k\Big>}{\beta_k=<\beta,\nu_k>}. \cr

The estimation of \eqn{\beta(t)} can be made by a few principal components (PC) under the assumption that \eqn{\beta_{k}=0}, for \eqn{k>k_n}  and the integral can be approximated by:\cr
    
\deqn{\hat{y}=\Big<X,\hat{\beta} \Big>=\sum_{k=1}^{k_{n}}{\gamma_{ik}\hat{\beta}_{k}}}{y.est= < X ,\beta.est > = \sum_{k=1:k_n} \gamma_{ik}\beta.est_k}
where, \deqn{\hat{\beta}_{(1:k_{n})}=\sum_{k=1}^{k_n}\frac{\gamma_{.k}^{T}y}{n\lambda_{k}}}{\beta.est_(1:k_n)=\sum_{k=1:k_n} \gamma_{.k}' y / (n\lambda_k)}
and  \eqn{\gamma_{(1:k_{n})}}{\gamma_(1:k_n)} is the (\eqn{n},\eqn{k_{n}}) matrix with \eqn{k_n} principal components estimation of \eqn{\beta} scores and \eqn{\lambda_i} the eigenvalues of the PC.\cr

\deqn{\hat{y}=\Big<X,\hat{\beta^{Ridge}} \Big>=\sum_{k=1}^{k_{n}}{\gamma_{ik}\hat{\beta}_{k}}}{y.est=<X,\beta.est^{Ridge}>=\sum_{k=1:k_n} \gamma_{ik} \beta.est_k}
where, \deqn{\hat{\beta}_{(1:k_{n})}=\sum_{k=1}^{k_n}\frac{\gamma_{.k}^{T}y}{n\lambda_{k}+rn}}{\beta.est_(1:k_{n})=\sum_{k=1:k_n} \gamma_{.k}'y / (n \lambda_k+rn)}

The response can also be fitted by: \deqn{\hat{y}=\sum_{k=1}^{k_n} \nu_k^{T}(\nu_k^{T}\nu_k)^{-1}\nu_k^{T}y}{y.est=\sum_{k=1:k_n} \nu_k'(\nu_k'\nu_k)^{-1}\nu_k'y}\cr
and for ridge regression by: \deqn{\hat{y}=\nu_k^{T}(\nu_k^{T}\nu_k+rnI)^{-1}\nu_k^{T}y}{y.est=\nu_k'(\nu_k'\nu_k+rnI)^{-1}\nu_k'y}
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

Cai TT, Hall P. 2006. \emph{Prediction in functional linear regression}. Annals of Statistics 34: 2159{-}2179.

Cardot H, Ferraty F, Sarda P. 1999. \emph{Functional linear model}. Statistics and Probability Letters 45: 11{-}22.

Hall P, Hosseini{-}Nasab M. 2006. \emph{On properties of functional principal components analysis}. Journal of the Royal Statistical Society B 68: 109{-}126.

Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}

}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
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