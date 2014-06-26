\name{fregre.gkam}
\Rdversion{1.1}
\alias{fregre.gkam}
\alias{kgam.H}
\title{ Fitting Functional Generalized Kernel Additive Models.}
\description{
	Computes functional regression between functional explanatory variables \eqn{(X^{1}(t_1),...,X^{q}(t_q))}{(X(t_1),...,X(t_q))} and scalar response \eqn{Y} using backfitting algorithm.
	
}
\usage{
fregre.gkam(formula,family = gaussian(),data, weights= rep(1,nobs),
par.metric = NULL,par.np=NULL,offset=NULL,
control = list(maxit = 100,epsilon = 0.001, trace = FALSE, inverse="solve"),...)
kgam.H(object,inverse="svd")
       }
\arguments{
\item{formula}{ an object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under \code{Details}.}
  \item{data}{ List that containing the variables in the model. }
 \item{family}{ a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions). }
  \item{weights}{ weights}
  \item{par.metric}{ List of arguments by covariate to pass to the \code{metric} function by covariate.}
  \item{par.np}{  List of arguments to pass to the \code{fregre.np.cv} function }
  \item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting.}	     
  \item{control}{ a list of parameters for controlling the fitting process, by default: \code{maxit}, \code{epsilon}, \code{trace} and \code{inverse}}
  \item{object}{ List that containing the Hat matrix for each variable in the model.}
  \item{inverse}{ ="svd" (by default) or ="solve" method.}
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{  

The smooth  functions  \eqn{f(.)} are estimated nonparametrically  using a iterative local scoring algorithm by applying  Nadaraya-Watson weighted kernel smoothers using \code{\link{fregre.np.cv}} in each step, see Febrero-Bande and Gonzalez-Manteiga (2011) for more details.  

Consider the fitted response \eqn{\hat{Y}=g^{-1}(H_{Q}y)}{g(Y.est)=Hy}, where \eqn{H_{Q}}{H} is the weighted hat matrix.\cr
Opsomer and Ruppert (1997)  solves a system  of equations for fit the unknowns \eqn{f(\cdot)}{f(.)} computing the additive smoother matrix \eqn{H_k}{H_k} such that \eqn{\hat{f}_k (X^k)=H_{k}Y}{f.est_k(X_k)=H_k Y} and \eqn{H_Q=H_1+,\cdots,+H_q}{H= H_1+,...,+H_q}. The additive model is fitted as follows:
     \deqn{\hat{Y}=g^{-1}\Big(\sum_i^q \hat{f_i}(X_i)\Big)}{g(y.est)=\sum(i:q) f.est_i(X_i)}

}
\value{
\item{result}{ List of non-parametric estimation by covariate.}
\item{fitted.values}{  Estimated scalar response.}
\item{residuals}{ \code{y} minus \code{fitted values}.}
\item{effects}{ The residual degrees of freedom.}
\item{alpha}{  Hat matrix.}
\item{family}{ Coefficient of determination.}
\item{linear.predictors}{ Residual variance.}
\item{deviance}{ Scalar response.}
\item{aic}{ Functional explanatory data.}
\item{null.deviance}{ Non functional explanatory data.}
\item{iter}{ Distance matrix between curves.}
\item{w}{ beta coefficient estimated}
\item{eqrank}{  List that containing the variables in the model.}
\item{prior.weights}{ Asymmetric kernel used.}
\item{y}{ Scalar response.}
\item{H}{ Hat matrix, see Opsomer and Ruppert(1997) for more details.}
\item{converged}{ conv.}
}
\references{
Febrero-Bande M. and Gonzalez-Manteiga W. (2012). \emph{Generalized Additive Models for Functional Data}. TEST. Springer-Velag. 
\url{http://dx.doi.org/10.1007/s11749-012-0308-0}

Opsomer J.D. and Ruppert D.(1997). \emph{Fitting a bivariate additive model by local polynomial regression}.Annals of Statistics, \code{25}, 186-211.
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:   \code{\link{fregre.gsam}}, \code{\link{fregre.glm}} and  \code{\link{fregre.np.cv}}\cr
}
\examples{
\dontrun{
data(tecator)
ab=tecator$absorp.fdata[1:100]
ab2=fdata.deriv(ab,2)
yfat=tecator$y[1:100,"Fat"]

# Example 1: # Changing the argument par.np and family
yfat.cat=ifelse(yfat<15,0,1)
xlist=list("df"=data.frame(yfat.cat),"ab"=ab,"ab2"=ab2)
f2<-yfat.cat~ab+ab2

par.NP<-list("ab"=list(Ker=AKer.norm,type.S="S.NW"),
"ab2"=list(Ker=AKer.norm,type.S="S.NW"))
res2=fregre.gkam(f2,family=binomial(),data=xlist,
par.np=par.NP)
res2

# Example 2: Changing the argument par.metric and family link
par.metric=list("ab"=list(metric=semimetric.deriv,nderiv=2,nbasis=15),
"ab2"=list("metric"=semimetric.basis))
res3=fregre.gkam(f2,family=binomial("probit"),data=xlist,
par.metric=par.metric,control=list(maxit=2,trace=FALSE))
summary(res3)

# Example 3: Gaussian family (by default)
# Only 1 iteration (by default maxit=100)
xlist=list("df"=data.frame(yfat),"ab"=ab,"ab2"=ab2)
f<-yfat~ab+ab2
res=fregre.gkam(f,data=xlist,control=list(maxit=1,trace=FALSE))
res
}
}
\keyword{regression}