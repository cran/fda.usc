% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.fregre.gkam.R, R/predict.fregre.glm.R,
%   R/predict.fregre.gsam.R, R/predict.fregre.lm.R, R/predict.fregre.plm.R
\name{predict.fregre.gkam}
\alias{predict.fregre.gkam}
\alias{predict.fregre.glm}
\alias{predict.fregre.gsam}
\alias{predict.fregre.lm}
\alias{predict.fregre.plm}
\title{Predict method for functional linear model}
\usage{
\method{predict}{fregre.gkam}(object, newx = NULL, type = "response", ...)

\method{predict}{fregre.glm}(object, newx = NULL, type = "response", ...)

\method{predict}{fregre.gsam}(object, newx = NULL, type = "response", ...)

\method{predict}{fregre.lm}(
  object,
  newx = NULL,
  type = "response",
  se.fit = FALSE,
  scale = NULL,
  df = df,
  interval = "none",
  level = 0.95,
  weights = 1,
  pred.var = res.var/weights,
  ...
)

\method{predict}{fregre.plm}(object, newx = NULL, ...)
}
\arguments{
\item{object}{\code{fregre.lm}, \code{fregre.plm}, \code{fregre.glm},
\code{fregre.gsam}\cr or \code{fregre.gkam} object.}

\item{newx}{An optional data list in which to look for variables with which
to predict. If omitted, the fitted values are used. List of new explanatory
data.}

\item{type}{a character vector, Type of prediction: (\code{response}, \code{terms} for model terms or \code{effects}  for model terms where
the partial effects are summarized for each functional variable.}

\item{\dots}{Further arguments passed to or from other methods.}

\item{se.fit}{=TRUE (not default) standard error estimates are returned for
each prediction.}

\item{scale}{Scale parameter for std.err. calculation.}

\item{df}{Degrees of freedom for scale.}

\item{interval}{Type of interval calculation.}

\item{level}{Tolerance/confidence level.}

\item{weights}{variance weights for prediction. This can be a numeric vector
or a one-sided model formula. In the latter case, it is interpreted as an
expression evaluated in newdata}

\item{pred.var}{the variance(s) for future observations to be assumed for
prediction intervals. See \code{link{predict.lm}} for more details.}
}
\value{
Return the predicted values and optionally:
\itemize{
\item {predict.lm,predict.glm,predict.gam}{ produces a vector of predictions
or a matrix of predictions and bounds with column names fit, lwr, and upr if
interval is set. If se.fit is TRUE, a list with the following components is
returned: fit vector or matrix as above.} 
\item {se.fit}{ standard error of predicted means.} 
\item {residual.scale}{ residual standard deviations.}
\item {df}{ degrees of freedom for residual.}
}
}
\description{
Computes predictions for regression between functional (and non functional)
explanatory variables and scalar response. 
\itemize{ 
\item \code{predict.fregre.lm}, Predict method for functional linear model of
\code{\link{fregre.lm}} fits object using basis or principal component
representation.
\item \code{predict.fregre.plm}, Predict method for
semi-functional linear regression model of \code{\link{fregre.plm}} fits
object using using asymmetric kernel estimation. 
\item \code{predict.fregre.glm}, Predict method for functional generalized linear
model of \code{\link{fregre.glm}} fits object using basis or principal
component representation. 
\item \code{predict.fregre.gsam}, Predict method for functional generalized 
spectral additive model of \code{\link{fregre.gsam}} fits object using basis 
or principal component representation.
\item \code{predict.fregre.gkam}, Predict method for functional generalized 
kernel additive model of \code{\link{fregre.gkam}} fits object using 
backfitting algorithm. 
}

These functions use the model fitting function \code{\link{lm}},
\code{\link{glm}} or \code{\link{gam}} properties.\cr If using functional
data derived, is recommended to use a number of bases to represent beta
lower than the number of bases used to represent the functional data. \cr
The first item in the \code{data} list of \code{newx} argument is called
\emph{"df"} and is a data frame with the response and non functional
explanatory variables, as \code{\link{lm}}, \code{\link{glm}} or
\code{\link{gam}}. Functional variables (\code{fdata} and \code{fd} class)
are introduced in the following items in the \code{data} list of \code{newx}
argument.
}
\examples{
\dontrun{
data(tecator)
ind <- 1:129
x <- tecator$absorp.fdata
x.d2 <- fdata.deriv(x,nderiv=2)
tt <- x[["argvals"]]
dataf <- as.data.frame(tecator$y)
ldat <- ldata("df" = dataf[ind,], "x.d2" = x.d2[ind])
basis.x <- list("x.d2" = create.pc.basis(ldat$x.d2))
res <- fregre.gsam(Fat ~  s(x.d2,k=3),
                   data=ldat, family = gaussian(),
                   basis.x = basis.x)
newldat <- ldata("df" = dataf[-ind,], "x.d2" = x.d2[-ind])
pred <- predict(res, newldat)
plot(pred,tecator$y$Fat[-ind])
res.glm <- fregre.glm(Fat  ~  x.d2, data = ldat,
                  family = gaussian(),basis.x = basis.x)
pred.glm <- predict(res.glm, newldat)
newy <- tecator$y$Fat[-ind]
points(pred.glm,tecator$y$Fat[-ind],col=2)

# Time-consuming 
res.gkam <- fregre.gkam(Fat ~ x.d2, data = ldat)
pred.gkam <- predict(res.gkam, newldat)
points(pred.gkam,tecator$y$Fat[-ind],col = 4)

((1/length(newy)) * sum((drop(newy)-pred)^2)) / var(newy)
((1/length(newy)) * sum((newy-pred.glm)^2)) / var(newy)    
((1/length(newy)) * sum((newy-pred.gkam)^2)) / var(newy)    
}                                                                                                              
}
\references{
Febrero-Bande, M., Oviedo de la Fuente, M. (2012).
\emph{Statistical Computing in Functional Data Analysis: The R Package
fda.usc.} Journal of Statistical Software, 51(4), 1-28.
\url{https://www.jstatsoft.org/v51/i04/}
}
\seealso{
See Also as: \code{\link{fregre.lm}}, \code{\link{fregre.plm}},
\code{\link{fregre.glm}}, \code{\link{fregre.gsam}} and
\code{\link{fregre.gkam}}.
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente
\email{manuel.oviedo@udc.es}
}
\keyword{regression}
