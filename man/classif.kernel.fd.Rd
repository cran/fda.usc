\name{classif.kernel.fd}
\Rdversion{1.1}
\alias{classif.kernel.fd}
\title{Kernel  Classifier from  Functional Data }
\description{
Fits Nonparametric Supervised Classification for Functional Data.  \cr}
\usage{
classif.kernel.fd(fdataobj, group, h = NULL, 
metric = metric.lp, ...)
}

\arguments{
 \item{fdataobj}{ \code{\link{fdata}} class object.} 
  \item{group}{ Factor of length \code{n}.}
	\item{h}{ Vector of smoothing parameter or bandwidth.}
  \item{metric}{ Metric function, by default \code{\link{metric.lp}}.}
  \item{\dots}{  Arguments to be passed  for \code{\link{metric.lp}} o other metric function and \code{\link{Kernel}} function.}
   }
  \details{
 Make the group classification of a training dataset using kernel estimation: \code{\link{Kernel}}.\cr  
  Different types of metric funtions can be used. }

\value{
 \item{fdataobj}{ \code{\link{fdata}} class object.}  
	\item{group}{  Factor of length \code{n}.}
  \item{group.est}{ Estimated vector groups}
  \item{prob.group}{ Matrix of predicted class probabilities. For each functional point shows the probability of each possible group membership.}
	\item{max.prob}{ Highest probability of correct classification.}
  \item{h.opt}{ Optimal smoothing parameter or bandwidht estimated.}   
	\item{D}{ Matrix of distances of the optimal quantile distance \code{hh.opt}.}	  
  \item{prob.classification}{ Probability of correct classification by group.}  	
    \item{misclassification}{ Vector of probability of misclassification by number of neighbors \code{knn}.}
  \item{h}{ Vector of smoothing parameter or bandwidht.}
   \item{C}{ A call of function \code{classif.kernel.fd}.}
}
\references{
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional data analysis.}
Springer Series in Statistics, New York.

Ferraty, F. and Vieu, P. (2006). \emph{NPFDA in practice}. Free access on line at
\url{http://www.lsp.ups-tlse.fr/staph/npfda/}
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}

\seealso{
See Also as \code{\link{predict.classif.fd}}}

\examples{

data(phoneme)
mlearn<-phoneme[["learn"]]
glearn<-phoneme[["classlearn"]]

h=seq(7,19,len=10)
out=classif.kernel.fd(mlearn,glearn,h=h)
summary.classif.fd(out)
#round(out$prob.group,4)

}

\keyword{classif}