\name{summary.classif.fd}
\Rdversion{1.1}
\alias{summary.classif.fd}
\alias{print.classif.fd}
\title{ Summarizes information from kernel classification methods. }
\description{
	Summary function for \code{\link{classif.knn.fd}}, \code{\link{classif.kernel.fd}} or \code{\link{classif.kernel.fb}}.
}
\usage{
%summary.classif.fd(object, ...)
%print.classif.fd(object,digits = max(3, getOption("digits") - 3), ...)
\method{summary}{classif.fd}(object,\dots)
\method{print}{classif.fd}(x,digits = max(3, getOption("digits") - 3),\dots)
}
\arguments{
  \item{object}{ Estimated by kernel classification.}
  \item{x}{ Estimated by kernel classification.}
  \item{digits}{a non-null value for digits specifies the minimum number of significant digits to be printed in values. The default, NULL, uses \code{\link{getOption}}(digits).}
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{
\code{object} from one of the following functions:\cr
\tabular{ll}{
\tab \code{\link{classif.knn.fd}}\cr
\tab \code{\link{classif.kernel.fd}}\cr
\tab \code{\link{classif.kernel.fb}}
}
}
\value{

Shows:\cr
\tabular{ll}{
\tab -Probability of correct classification by group \code{prob.classification}.\cr
\tab -Confusion matrix between the theoretical groups  and estimated groups.\cr
\tab -Highest probability of correct classification \code{max.prob}.\cr 
}
If the object is returned from the function: \code{\link{classif.knn.fd}}\cr
\tabular{ll}{
\tab -Vector of probability of correct classification by number of neighbors \code{knn}.\cr
\tab -Optimal number of neighbors: \code{knn.opt}.\cr
}
If the object is returned from the function: \code{\link{classif.kernel.fd}}\cr
\tabular{ll}{
\tab -Vector of probability of correct classification by banwidth \code{h}.\cr
\tab -Functional measure of closeness (optimal distance, \code{h.opt}).\cr
}
If the object is returned from the function: \code{\link{classif.kernel.fb}}\cr
\tabular{ll}{
\tab -Matrix of probability of correct classification by number of basis \code{num.basis} and bandwidth \code{h}.\cr
\tab -Functional measure of proximity(optimal distance, \code{h.opt}).\cr
\tab -Optimal number of basis: \code{basis.opt} and optimal bandwidth: \code{h.opt}.\cr
}
  \item{object}{ Estimated by kernel classification.}
}  
\author{Febrero-Bande, M. and Oviedo de la Fuente, M.}
\seealso{
See Also as:
\code{\link{classif.knn.fd}}, \code{\link{classif.kernel.fd}}, \cr
\code{\link{classif.kernel.fb}} and \code{\link{summary.classif.fd}}}
\examples{
data(phoneme)
mlearn<-phoneme[["learn"]]
glearn<-phoneme[["classlearn"]]
out=classif.knn.fd(mlearn,glearn,knn=c(3,5,7))
summary.classif.fd(out)
out2=classif.kernel.fd(mlearn,glearn,h=2^(0:5))
summary.classif.fd(out2)
out3=classif.kernel.fb(mlearn,glearn,h=2^(0:5),
par.basis = list(nbasis=c(7,17,27,37)))
summary.classif.fd(out3)
}
\keyword{print}