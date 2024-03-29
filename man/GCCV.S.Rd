% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GCCV.S.R
\name{GCCV.S}
\alias{GCCV.S}
\title{The generalized correlated cross-validation (GCCV) score.}
\usage{
GCCV.S(
  y,
  S,
  criteria = "GCCV1",
  W = NULL,
  trim = 0,
  draw = FALSE,
  metric = metric.lp,
  ...
)
}
\arguments{
\item{y}{Response vectorith length \code{n} or Matrix of set cases with
dimension (\code{n} x \code{m}), where \code{n} is the number of curves and
\code{m} are the points observed in each curve.}

\item{S}{Smoothing matrix, see \code{\link{S.NW}}, \code{\link{S.LLR}} or
\eqn{S.KNN}.}

\item{criteria}{The penalizing function. By default \emph{"Rice"} criteria.
"GCCV1","GCCV2","GCCV3","GCV") Possible values are \emph{"GCCV1"},
\emph{"GCCV2"}, \emph{"GCCV3"}, \emph{"GCV"}.}

\item{W}{Matrix of weights.}

\item{trim}{The alpha of the trimming.}

\item{draw}{=TRUE, draw the curves, the sample median and trimmed mean.}

\item{metric}{Metric function, by default \code{\link{metric.lp}}.}

\item{\dots}{Further arguments passed to or from other methods.}
}
\value{
Returns GCCV score calculated for input parameters.
}
\description{
The generalized correlated cross-validation (GCV) score.
}
\details{
\deqn{ }{\sum(y-y.fit)^2 / (1-tr(C)/n)^2}\deqn{GCCV=\frac{\sum_{i=1}^n
{y_{i}-\hat{y}_{i,b}}^2}{1-\frac{tr(C)}{n}^2} }{\sum(y-y.fit)^2 /
(1-tr(C)/n)^2} %where \eqn{C=2{S\Sigma(\theta)}-{S\Sigma(\theta)S'}} \cr
%and \eqn{\Sigma=\sigma C} is the n x n covariance matrix with
\eqn{cor(\epsilon_i,\epsilon_j ) =\sigma}\cr

where \eqn{S} is the smoothing matrix \eqn{S} and:\cr A.-If \eqn{C=2S\Sigma
- S\Sigma S} \cr B.-If \eqn{C=S\Sigma} \cr C.-If \eqn{C=S\Sigma S'} \cr with
\eqn{\Sigma} is the n x n covariance matrix with
\eqn{cor(\epsilon_i,\epsilon_j ) =\sigma}
}
\note{
Provided that \eqn{C = I} and the smoother matrix S is symmetric and
idempotent, as is the case for many linear fitting techniques, the trace
term reduces to \eqn{n - tr[S]}, which is proportional to the familiar
denominator in GCV.
}
\examples{
\dontrun{
data(tecator)
x=tecator$absorp.fdata
x.d2<-fdata.deriv(x,nderiv=)
tt<-x[["argvals"]]
dataf=as.data.frame(tecator$y)
y=tecator$y$Fat
# plot the response
plot(ts(tecator$y$Fat))

nbasis.x=11;nbasis.b=7
basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)
basis.x=list("x.d2"=basis1)
basis.b=list("x.d2"=basis2)
ldata=list("df"=dataf,"x.d2"=x.d2)
# No correlation
res.gls=fregre.gls(Fat~x.d2,data=ldata, 
                   basis.x=basis.x,basis.b=basis.b)
# AR1 correlation                   
res.gls=fregre.gls(Fat~x.d2,data=ldata, correlation=corAR1(),
                   basis.x=basis.x,basis.b=basis.b)
GCCV.S(y,res.gls$H,"GCCV1",W=res.gls$W)
res.gls$gcv
}
}
\references{
Carmack, P. S., Spence, J. S., and Schucany, W. R. (2012).
Generalised correlated cross-validation. Journal of Nonparametric
Statistics, 24(2):269--282.

Oviedo de la Fuente, M., Febrero-Bande, M., Pilar Munoz, and Dominguez, A.
Predicting seasonal influenza transmission using Functional Regression
Models with Temporal Dependence. arXiv:1610.08718.
\url{https://arxiv.org/abs/1610.08718}
}
\seealso{
See Also as \code{\link{optim.np}}. \cr Alternative method
(independent case): \code{\link{GCV.S}}
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente
\email{manuel.oviedo@udc.es}
}
\keyword{utilities}
