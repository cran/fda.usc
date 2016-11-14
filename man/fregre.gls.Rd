\name{fregre.gls}
\title{Fit Functional Linear Model Using Generalized Least Squares}
\alias{fregre.gls}
\usage{
fregre.gls(formula, data, correlation = NULL, basis.x = NULL, 
    basis.b = NULL, rn, lambda, weights = NULL, subset, method = c("REML", 
        "ML"), control = list(), verbose = FALSE, criteria="GCCV1",...) 
}
\arguments{
  \item{formula}{a two-sided linear formula object describing the
    model, with the response on the left of a \code{~} operator and the
    terms, separated by \code{+} operators, on the right.}
 \item{data}{an optional data frame containing the variables named in
   \code{model}, \code{correlation}, \code{weights}, and
   \code{subset}. By default the variables are taken from the
   environment from which \code{gls} is called.}
 \item{correlation}{an optional \code{\link{corStruct}} object describing the
   within-group correlation structure. See the documentation of
   \code{\link{corClasses}} for a description of the available \code{corStruct}
   classes. If a grouping variable is to be used, it must be specified in
   the \code{form} argument to the \code{corStruct}
   constructor. Defaults to \code{NULL}, corresponding to uncorrelated 
   errors.}  
  \item{basis.x}{ List of basis for functional explanatory data estimation.}
	\item{basis.b}{ List of basis for \eqn{\beta(t)} parameter estimation.}
 	\item{rn}{List of Ridge parameter.} 
  \item{lambda}{List of Roughness penalty parameter.} 
  \item{weights}{an optional \code{\link{varFunc}} object or one-sided formula
   describing the within-group heteroscedasticity structure. If given as
   a formula, it is used as the argument to \code{\link{varFixed}},
   corresponding to fixed variance weights. See the documentation on
   \code{\link{varClasses}} for a description of the available \code{\link{varFunc}}
   classes. Defaults to \code{NULL}, corresponding to homoscedastic
   errors.} 
 \item{subset}{an optional expression indicating which subset of the rows of
   \code{data} should  be  used in the fit. This can be a logical
   vector, or a numeric vector indicating which observation numbers are
   to be included, or a  character  vector of the row names to be
   included.  All observations are included by default.}
 \item{method}{a character string.  If \code{"REML"} the model is fit by
   maximizing the restricted log-likelihood.  If \code{"ML"} the
   log-likelihood is maximized.  Defaults to \code{"REML"}.}
 %\item{na.action}{a function that indicates what should happen when the data contain \code{NA}s.  The default action (\code{\link{na.fail}}) causes\code{gls} to print an error message and terminate if there are any incomplete observations.}
 \item{control}{a list of control values for the estimation algorithm to
   replace the default values returned by the function \code{\link{glsControl}}.
   Defaults to an empty list.}
 \item{verbose}{an optional logical value. If \code{TRUE} information on
   the evolution of the iterative algorithm is printed. Default is
   \code{FALSE}.}
 \item{criteria}{ GCCV criteria, see \code{\link{GCCV.S}}.}    
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this methodl.} 
 %\item{evaluate}{If \code{TRUE} evaluate the new call else return the call.}
}
\description{
  This function fits a functional linear model using generalized least
  squares. The errors are allowed to be correlated and/or have unequal
  variances.  
}
\value{
  an object of class \code{"gls"} representing the functional linear model
  fit. Generic functions such as \code{print}, \code{plot}, and 
  \code{summary} have methods to show the results of the fit. 
  
  See \code{\link{glsObject}} for the components of the fit. The functions
  \code{\link{resid}}, \code{\link{coef}} and \code{\link{fitted}},
  can be used to extract some of its components. 

Beside, the class(z) is  "gls", "lm" and "fregre.lm" with the following objects:
\item{sr2}{ Residual variance.}
\item{Vp}{ Estimated covariance matrix for the parameters.}
\item{lambda}{ A roughness penalty.}	
\item{basis.x}{ Basis used for \code{fdata} or \code{fd} covariates.}
\item{basis.b}{ Basis used for beta parameter estimation.}
\item{beta.l}{ List of estimated beta parameter of functional covariates.}
\item{data}{ List that containing the variables in the model.}
\item{formula}{ formula used in ajusted model.}
\item{formula.ini}{ formula in call.}
%\item{mean}{list of mean of predictors }
%\item{JJ}{JJ}
\item{XX}{ desing matrix }
\item{W}{ inverse of covariance matrix}
%\item{fdataob}{ }
%\item{rn}{ rn}
%\item{vs.list}{ }
\item{correlation}{ See glsObject for the components of the fit. }
}
\references{
 Oviedo de la Fuente, M., Febrero-Bande, M., Pilar Munoz, and Dominguez, A. Predicting seasonal influenza transmission using Functional Regression Models with Temporal Dependence. arXiv:1610.08718. \url{https://arxiv.org/abs/1610.08718}}

\examples{

\dontrun{ 
data(tecator)
x=tecator$absorp.fdata
x.d2<-fdata.deriv(x,nderiv=)
tt<-x[["argvals"]]
dataf=as.data.frame(tecator$y)

# plot the response
plot(ts(tecator$y$Fat))

nbasis.x=11;nbasis.b=7
basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)
basis.x=list("x.d2"=basis1)
basis.b=list("x.d2"=basis2)
ldata=list("df"=dataf,"x.d2"=x.d2)
res.gls=fregre.gls(Fat~x.d2,data=ldata, correlation=corAR1(),
                   basis.x=basis.x,basis.b=basis.b)
summary(res.gls)                   
}

}
\keyword{regression,models}

