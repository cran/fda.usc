\name{summary.classif}
\Rdversion{1.1}
\alias{summary.classif}
\alias{print.classif}
\title{ Summarizes information from kernel classification methods. }
\description{
	Summary function for \code{\link{classif.knn}} or \code{\link{classif.kernel}}.
}
\usage{
%summary.classif.fd(object, ...)
%print.classif.fd(object,digits = max(3, getOption("digits") - 3), ...)
\method{summary}{classif}(object,\dots)
\method{print}{classif}(x,digits = max(3, getOption("digits") - 3),\dots)
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
\tab \code{\link{classif.knn}}\cr
\tab \code{\link{classif.kernel}}\cr
}
}
\value{

Shows:\cr
\tabular{ll}{
\tab -Probability of correct classification by group \code{prob.classification}.\cr
\tab -Confusion matrix between the theoretical groups  and estimated groups.\cr
\tab -Highest probability of correct classification \code{max.prob}.\cr 
}
If the object is returned from the function: \code{\link{classif.knn}}\cr
\tabular{ll}{
\tab -Vector of probability of correct classification by number of neighbors \code{knn}.\cr
\tab -Optimal number of neighbors: \code{knn.opt}.\cr
}
If the object is returned from the function: \code{\link{classif.kernel}}\cr
\tabular{ll}{
\tab -Vector of probability of correct classification by banwidth \code{h}.\cr
\tab -Functional measure of closeness (optimal distance, \code{h.opt}).\cr
}
  \item{object}{ Estimated by kernel classification.}
}  
\author{Febrero-Bande, M. and Oviedo de la Fuente, M.}
\seealso{
See Also as:
\code{\link{classif.knn}}, \code{\link{classif.kernel}} \cr
 and \code{\link{summary.classif}}}
\examples{
data(phoneme)
mlearn<-phoneme[["learn"]]
glearn<-phoneme[["classlearn"]]
out=classif.knn(glearn,mlearn,knn=c(3,5,7))
summary.classif(out)
out2=classif.kernel(glearn,mlearn,h=2^(0:5))
summary.classif(out2)
}
\keyword{print}