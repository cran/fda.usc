% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/semimetric.NPFDA.r
\name{semimetric.NPFDA}
\alias{semimetric.NPFDA}
\alias{semimetric.deriv}
\alias{semimetric.fourier}
\alias{semimetric.hshift}
\alias{semimetric.mplsr}
\alias{semimetric.pca}
\title{Proximities between functional data (semi-metrics)}
\source{
\url{https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/}
}
\usage{
semimetric.deriv(
  fdata1,
  fdata2 = fdata1,
  nderiv = 1,
  nknot = ifelse(floor(ncol(DATA1)/3) > floor((ncol(DATA1) - nderiv - 4)/2),
    floor((ncol(DATA1) - nderiv - 4)/2), floor(ncol(DATA1)/3)),
  ...
)

semimetric.fourier(
  fdata1,
  fdata2 = fdata1,
  nderiv = 0,
  nbasis = ifelse(floor(ncol(DATA1)/3) > floor((ncol(DATA1) - nderiv - 4)/2),
    floor((ncol(DATA1) - nderiv - 4)/2), floor(ncol(DATA1)/3)),
  period = NULL,
  ...
)

semimetric.hshift(fdata1, fdata2 = fdata1, t = 1:ncol(DATA1), ...)

semimetric.mplsr(fdata1, fdata2 = fdata1, q = 2, class1, ...)

semimetric.pca(fdata1, fdata2 = fdata1, q = 1, ...)
}
\arguments{
\item{fdata1}{Functional data 1 or curve 1. \code{DATA1} with dimension
(\code{n1} x \code{m}), where \code{n1} is the number of curves and \code{m}
are the points observed in each curve.}

\item{fdata2}{Functional data 2 or curve 2. \code{DATA1} with dimension
(\code{n2} x \code{m}), where \code{n2} is the number of curves and \code{m}
are the points observed in each curve.}

\item{nderiv}{Order of derivation, used in \code{semimetric.deriv} and \cr
\code{semimetric.fourier}}

\item{nknot}{semimetric.deriv argument: number of interior knots (needed for
defining the B-spline basis).}

\item{\dots}{Further arguments passed to or from other methods.}

\item{nbasis}{\code{semimetric.fourier}: size of the basis.}

\item{period}{\code{semimetric.fourier}:allows to select the period for the
fourier expansion.}

\item{t}{\code{semimetric.hshift}: vector which defines \code{t} (one can
choose \code{1,2,...,nbt} where \code{nbt} is the number of points of the
discretization)}

\item{q}{If \code{semimetric.pca}: the retained number of principal
components.\cr If \code{semimetric.mplsr}: the retained number of factors.}

\item{class1}{\code{semimetric.mplsr}: vector containing a categorical
response which corresponds to class number for units stored in \code{DATA1}.}
}
\value{
Returns a proximities matrix between two functional datasets.
}
\description{
Computes semi-metric distances of functional data based on Ferraty F and
Vieu, P. (2006).
}
\details{
\code{semimetric.deriv}: approximates \eqn{L_2} metric
between derivatives of the curves based on ther B-spline representation. The
derivatives set with the argument \code{nderiv}.\cr
\code{semimetric.fourier}: approximates \eqn{L_2} metric between the curves
based on ther B-spline representation. The derivatives set with the argument
\code{nderiv}.\cr \code{semimetric.hshift}: computes distance between curves
taking into account an horizontal shift effect.\cr \code{semimetric.mplsr}:
computes distance between curves based on the partial least squares
method.\cr \code{semimetric.pca}: computes distance between curves based on
the functional principal components analysis method.

In the next semi-metric functions the functional data \eqn{X} is
approximated by \eqn{k_n} elements of the Fourier, B--spline, PC or PLS basis
using, \eqn{\hat{X_i} =\sum_{k=1}^{k_n}\nu_{k,i}\xi_k}, where \eqn{\nu_k}
are the coefficient of the expansion on the basis function
\eqn{\left\{\xi_k\right\}_{k=1}^{\infty}}.\cr The distances between the q-order derivatives of two curves \eqn{X_{1}} and
\eqn{X_2} is,
\deqn{d_{2}^{(q)}\left(X_1,X_2\right)_{k_n}=\sqrt{\frac{1}{T}\int_{T}\left(X_{1}^{(q)}(t)-X_{2}^{(q)}(t)\right)^2
dt}} where \eqn{X_{i}^{(q)}\left(t\right)} denot the \eqn{q} derivative of
\eqn{X_i}.

\code{semimetric.deriv} and \code{semimetric.fourier} function use a
B-spline and Fourier approximation respectively for each curve and the
derivatives are directly computed by differentiating several times their
analytic form, by default \code{q=1} and \code{q=0} respectively.
\code{semimetric.pca} and \code{semimetric.mprls} function compute
proximities between curves based on the functional principal components
analysis (FPCA) and the functional partial least square analysis (FPLS),
respectively. The FPC and FPLS reduce the functional data in a reduced
dimensional space (q components). \code{semimetric.mprls} function requires
a scalar response.

\deqn{d_{2}^{(q)}\left(X_1,X_2\right)_{k_n}\approx\sqrt{\sum_{k=1}^{k_n}\left(\nu_{k,1}-\nu_{k,2}\right)^2\left\|\xi_k^{(q)}\right\|dt}}
\code{semimetric.hshift} computes proximities between curves taking into
account an horizontal shift effect.

\deqn{d_{hshift}\left(X_1,X_2\right)=\min_{h\in\left[-mh,mh\right]}d_2(X_1(t),X_2(t+h))}
where \eqn{mh} is the maximum horizontal shifted allowed.
}
\examples{
\dontrun{ 
#	INFERENCE PHONDAT
data(phoneme)
ind=1:100 # 2 groups
mlearn<-phoneme$learn[ind,]
mtest<-phoneme$test[ind,]
n=nrow(mlearn[["data"]])
np=ncol(mlearn[["data"]])
mdist1=semimetric.pca(mlearn,mtest)
mdist2=semimetric.pca(mlearn,mtest,q=2)
mdist3=semimetric.deriv(mlearn,mtest,nderiv=0)
mdist4=semimetric.fourier(mlearn,mtest,nderiv=2,nbasis=21)
#uses hshift function
#mdist5=semimetric.hshift(mlearn,mtest) #takes a lot
glearn<-phoneme$classlearn[ind]
#uses mplsr function
mdist6=semimetric.mplsr(mlearn,mtest,5,glearn)
mdist0=metric.lp(mlearn,mtest)
b=as.dist(mdist6)
c2=hclust(b)
plot(c2)
memb <- cutree(c2, k = 2)
table(memb,phoneme$classlearn[ind])
 } 
  
}
\references{
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional
data analysis.} Springer Series in Statistics, New York.

Ferraty, F. and Vieu, P. (2006). \emph{NPFDA in practice}.  Free access on
line at \url{https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/}
}
\seealso{
See also \code{\link{metric.lp}} and \code{\link{semimetric.basis}}
}
\keyword{cluster}
