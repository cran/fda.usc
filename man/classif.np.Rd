\name{classif.np}
\Rdversion{1.1}
\alias{classif.np}
\alias{classif.kernel}
\alias{classif.knn}
\title{Kernel  Classifier from  Functional Data }
\description{
Fits Nonparametric Supervised Classification for Functional Data.  }
\usage{
classif.np(group,fdataobj,h=NULL,Ker=AKer.norm,metric=metric.lp,
type.CV = GCV.S,type.S=S.NW,par.CV=list(trim=0),par.S=list(),...)
classif.knn(group,fdataobj,knn=NULL,metric=metric.lp,
type.CV = GCV.S,par.CV=list(trim=0),par.S=list(),...)
classif.kernel(group,fdataobj,h=NULL,Ker=AKer.norm,metric=metric.lp,
type.CV = GCV.S,par.CV=list(trim=0),par.S=list(),...)
}
\arguments{
 \item{group}{ Factor of length \emph{n}}
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
%   \item{w}{ Optional case weights.}
	\item{h}{ Vector of smoothing parameter or bandwidth.}
	\item{knn}{ Vector of number of nearest neighbors considered.}
	\item{Ker}{ Type of kernel used.}	
  \item{metric}{ Metric function, by default \code{\link{metric.lp}}.}
   \item{type.CV}{ Type of cross-validation. By default generalized cross-validation \code{\link{GCV.S}} method.} 
  \item{type.S}{ Type of smothing matrix \code{S}. By default \code{S} is calculated by Nadaraya-Watson kernel estimator (\code{S.NW}).}
   \item{par.CV}{ List of parameters for \code{type.CV}: \code{trim}, the alpha of the trimming and \code{draw=TRUE}.}	  
  \item{par.S}{ List of parameters for \code{type.S}: \code{w}, the weights.}	      
   \item{\dots}{  Arguments to be passed  for \code{\link{metric.lp}} o other metric function and \code{\link{Kernel}} function.}
   }
  \details{
 Make the group classification of a training dataset using kernel or KNN estimation: \code{\link{Kernel}}.\cr  
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
   \item{C}{ A call of function \code{classif.kernel}.}
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
See Also as \code{\link{predict.classif}}}

\examples{
data(phoneme)
mlearn<-phoneme[["learn"]]
glearn<-phoneme[["classlearn"]]

h=9:19
out=classif.np(glearn,mlearn,h=h)
summary.classif(out)
#round(out$prob.group,4)
}
\keyword{classif}