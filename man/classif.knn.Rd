\name{classif.knn}
\Rdversion{1.1}
\alias{classif.knn}
\title{ k-Nearest Neighbor Classifier from Functional Data}
\description{
Performs k-nearest neighbor classification of a test set using a training set.   
For each row of the test set, the k nearest training set vectors are found.\cr \cr
Make supervised classification for training functional data via k-nearest neighbors method. \cr \cr
The method classifies the functional data to the group with the highest number of nearest neighbors.  
In case of tie the data is classified in the group with a shorter distance.
}
\usage{
classif.knn(group, fdataobj, w = NULL, 
knn = seq(3, floor(min(table(group))/3),by = 2),metric=metric.lp,...) 
}
\arguments{
 \item{group}{ Factor of length \emph{n}}
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
  \item{w}{ Optional case weights.}
   \item{knn}{ Vector of number of nearest neighbors considered }
  \item{metric}{ Metric function, by default \code{\link{metric.lp}}}
  \item{\dots}{ Further arguments passed to or from other methods. Arguments to be passed  for  metric function.}
}
\details{
The nearest neighbor method \code{\link{classif.knn}} to estimate the number of k nearest neighbors to classify a sample of best practice. In the training sample is known for each data functional group (group). You can use different metric functions or types of distance metric by changing the parameters: \emph{p} and \emph{w} of \code{\link{metric.lp}} function.
}

\value{
 \item{fdataobj}{ \code{\link{fdata}} class object.}  
  \item{group}{ Factor of length \emph{n}}
  \item{group.est}{ Vector of classes of the k nearest neighbors.}
  \item{prob.group}{ Matrix of predicted class probabilities. For each functional point shows the probability of each possible group membership.}  
  \item{prob.classification}{ Probability of correct classification by group for the estimated number of k neighbors \code{knn.opt}.}
  \item{knn.opt}{ Number of neighbors with highest probability of correct classification.}
  \item{D}{ Matrix of distances of the optimal k nearest neighbors \code{knn.opt}.}
  \item{max.prob}{ Highest probability of correct classification.}
  \item{misclassification}{ Vector of probability of misclassification by number of neighbors \code{knn}.}
  \item{knn}{ Vector of number of k nearest neighbors considered.}
  \item{C}{ A call of function \code{classif.knn}.}
  }
\references{
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional data analysis.}
Springer Series in Statistics, New York.\cr

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

out=classif.knn(glearn,mlearn)
summary.classif(out)
}
\keyword{classif}