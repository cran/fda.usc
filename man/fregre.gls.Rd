% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fregre.gls.r
\name{fregre.gls}
\alias{fregre.gls}
\title{Fit Functional Linear Model Using Generalized Least Squares}
\usage{
fregre.gls(
  formula,
  data,
  correlation = NULL,
  basis.x = NULL,
  basis.b = NULL,
  rn,
  lambda,
  weights = NULL,
  subset,
  method = c("REML", "ML"),
  control = list(),
  verbose = FALSE,
  criteria = "GCCV1",
  ...
)
}
\arguments{
\item{formula}{a two-sided linear formula object describing the model, with
the response on the left of a \code{~} operator and the terms, separated by
\code{+} operators, on the right.}

\item{data}{an optional data frame containing the variables named in
\code{model}, \code{correlation}, \code{weights}, and \code{subset}. By
default the variables are taken from the environment from which \code{gls}
is called.}

\item{correlation}{an optional \code{\link{corStruct}} object describing the
within-group correlation structure. See the documentation of
\code{\link{corClasses}} for a description of the available \code{corStruct}
classes. If a grouping variable is to be used, it must be specified in the
\code{form} argument to the \code{corStruct} constructor. Defaults to
\code{NULL}, corresponding to uncorrelated errors.}

\item{basis.x}{List of basis for functional explanatory data estimation.}

\item{basis.b}{List of basis for \eqn{\beta(t)} parameter estimation.}

\item{rn}{List of Ridge parameter.}

\item{lambda}{List of Roughness penalty parameter.}

\item{weights}{an optional \code{\link{varFunc}} object or one-sided formula
describing the within-group heteroscedasticity structure. If given as a
formula, it is used as the argument to \code{\link{varFixed}}, corresponding
to fixed variance weights. See the documentation on \code{\link{varClasses}}
for a description of the available \code{\link{varFunc}} classes. Defaults
to \code{NULL}, corresponding to homoscedastic errors.}

\item{subset}{an optional expression indicating which subset of the rows of
\code{data} should be used in the fit. This can be a logical vector, or a
numeric vector indicating which observation numbers are to be included, or a
character vector of the row names to be included.  All observations are
included by default.}

\item{method}{a character string.  If \code{"REML"} the model is fit by
maximizing the restricted log-likelihood.  If \code{"ML"} the log-likelihood
is maximized.  Defaults to \code{"REML"}.}

\item{control}{a list of control values for the estimation algorithm to
replace the default values returned by the function
\code{\link{glsControl}}.  Defaults to an empty list.}

\item{verbose}{an optional logical value. If \code{TRUE} information on the
evolution of the iterative algorithm is printed. Default is \code{FALSE}.}

\item{criteria}{GCCV criteria, see \code{\link{GCCV.S}}.}

\item{\dots}{some methods for this generic require additional arguments.
None are used in this methodl.}
}
\value{
an object of class \code{"gls"} representing the functional linear
model fit. Generic functions such as \code{print}, \code{plot}, and
\code{summary} have methods to show the results of the fit.\cr 
See \code{\link{glsObject}} for the components of the fit. The functions
\code{\link{resid}}, \code{\link{coef}} and \code{\link{fitted}}, can be
used to extract some of its components.\cr 
Beside, the class(z) is "gls", "lm" and "fregre.lm" with the following
objects: 
\itemize{
\item \code{sr2:}{ Residual variance.} 
\item \code{Vp:}{ Estimated covariance matrix for the parameters.} 
\item \code{lambda:}{ A roughness penalty.}
\item \code{basis.x:}{ Basis used for \code{fdata} or \code{fd} covariates.}
\item \code{basis.b:}{ Basis used for beta parameter estimation.} 
\item \code{beta.l:}{ List of estimated beta parameter of functional covariates.} 
\item \code{data:}{ List that containing the variables in the model.} 
\item \code{formula:}{ formula used in ajusted model.} 
\item \code{formula.ini:}{ formula in call.}
\item \code{W:}{ inverse of covariance matrix} 
\item \code{correlation:}{ See glsObject for the components of the fit. }
}
}
\description{
This function fits a functional linear model using generalized least
squares. The errors are allowed to be correlated and/or have unequal
variances.
}
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
\references{
Oviedo de la Fuente, M., Febrero-Bande, M., Pilar Munoz, and
Dominguez, A.  (2018). Predicting seasonal influenza transmission using 
functional regression models with temporal dependence. PloS one, 13(4), e0194250.
\doi{10.1371/journal.pone.0194250}
}
\keyword{models}
\keyword{regression}
