\name{fregre.kgam}
\Rdversion{1.1}
\alias{fregre.kgam}
\alias{kgam.H}
\title{ Fitting Functional Generalized Additive Models.}
\description{
	Computes functional regression between functional  explanatory variables and scalar response using backfitting algorithm.
}
\usage{
fregre.kgam(formula,family = gaussian(),data, weights= rep(1,nobs),
     par.metric = NULL,par.np=NULL,offset=NULL,
     control = list(maxit = 100,epsilon = 0.001, trace = FALSE, inverse="solve"),...)
kgam.H(object,inverse="solve")
       }
\arguments{
\item{formula}{ an object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under \code{Details}.}
  \item{data}{ List that containing the variables in the model. }
 \item{family}{ a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions). }
  \item{weights}{ weights}
  \item{par.metric}{ List of arguments by covariable to pass to the \code{metric} function by covariable.}
  \item{par.np}{  List of arguments to pass to the \code{fregre.np.cv} function }
  \item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting.}	     
  \item{control}{ a list of parameters for controlling the fitting process, by default: \code{maxit}, \code{epsilon}, \code{trace} and \code{inverse}}
  \item{object}{ List that containing the Hat matrix for each variable in the model.}
  \item{inverse}{ ="solve" (by default) or ="svd" method.}
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{  
The \code{kgam.H} function returns the hat matrix. Joint estimation of the hat matrix is performed combined  the hat matrix calculated in the \code{\link{fregre.np.cv}} function.
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
\item{H}{ Hat matrix.}
\item{converged}{ conv.}
}
\references{
Febrero-Bande M. and Gonzalez-Maneiga W. (2011). \emph{Generalized Additive Models for Functional Data}. Recent Advances in functional Data Analyisis and Related Topics. Physica-Velag. \cr
McCullagh and Nelder (1989), \emph{Generalized Linear Models} 2nd ed. Chapman and Hall. \cr
Opsomer J.D. and Ruppert D.(1997). \emph{Fitting a bivariate additive model by local polynomial regression}.Annals of Statistics, \code{25}, 186-211.
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:  \code{\link{fregre.glm}} and  \code{\link{fregre.np.cv}}\cr
}
\examples{
data(tecator)
ab=tecator$absorp.fdata
ab1=fdata.deriv(ab,1)
ab2=fdata.deriv(ab,2)
yfat=as.integer(cut(tecator$y[,"Fat"],c(0,15,100)))-1
# Ajuste simple
xlist=list("df"=data.frame(yfat),"ab2"=ab2,"ab"=ab)
f<-yfat~ab+ab2
res=fregre.kgam(f,data=xlist,family=binomial("logit"),control=list(maxit=2))
}
\keyword{regression}