\name{classif.kernel.fb}
\Rdversion{1.1}
\alias{classif.kernel.fb}
\title{	Kernel classifier from Functional Data Training by basis representation}
\description{
 Fits Nonparametric Supervised Classification for Functional Data Training using a basis representation.  \cr
}
\usage{
classif.kernel.fb(fdataobj,group,h=NULL,metric=metric.lp,
type.basis="bspline",par.basis=list(rangeval=NULL,nbasis=NULL),...)
}
\arguments{
\item{fdataobj}{ \code{\link{fdata}} class object.}  
   \item{group}{ Vector group of length n. The vector should be a factor, if not a factor becomes a factor.}
  \item{h}{ Vector of smoothing parameter or bandwidht.}
  \item{metric}{ Metric function, by default \code{\link{metric.lp}}.}
  \item{type.basis}{ Character string which determines type of basis. By default "bspline".}
  \item{par.basis}{  List of parameters for \code{\link{create.basis}}}
  \item{\dots}{ Further arguments passed to or from other methods. \cr
   Arguments to be passed  for \code{\link{metric.lp}} o other metric function.}
  }
\details{
Classification is done using a kernel estimation after a representation in a basis using \link[fda]{Data2fd} function of \code{fda} package.\cr
Different types of metric funtions and type of basis (see \link[fda]{create.basis}) can be used.
}
\value{	
 \item{fdataobj}{ \code{\link{fdata}} class object.}  
	\item{group}{ Vector group of length \code{n} (as a factor).}
  \item{group.est}{ Estimated vector groups.}
  \item{max.prob}{ Highest probability of correct classification.}  
  \item{h.opt}{ Optimal smoothing parameter or bandwidht estimated.}  
  \item{prob.classification}{ Probability of correct classification by group.}
  \item{prob.group}{ Matrix of estimated class probabilities. For each functional point shows the probability of each possible group membership.}  
  \item{basis.opt}{ Optimal number of basis.}
  \item{basis.obj.opt}{ Optimal basis object.}
	\item{D}{ Matrix of distances between objects.}  
  \item{coefs}{ Coefficient array  of Functional Data  object, \link[fda]{Data2fd}.}
	\item{misclassification}{ Vector of probability of misclassification by number of neighbors \code{knn}.}
  \item{num.basis}{ Vector of number of basis considered.}
  \item{h}{ Vector of smoothing parameter or bandwidth.}
  \item{C}{ A call of function \code{classif.kernel.fb}.}
}
           

\references{
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional data analysis.}
Springer Series in Statistics, New York.

Ferraty, F. and Vieu, P. (2006). \emph{NPFDA in practice}. Free access on line at
\url{http://www.lsp.ups-tlse.fr/staph/npfda/}

Ramsay, James O., and Silverman, Bernard W. (2006), \emph{ Functional Data Analysis}, 2nd ed., Springer, New York. 

Ramsay, James O., and Silverman, Bernard W. (2002), \emph{ Applied Functional Data Analysis}, Springer, New York. 
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
vbasis<-seq(7,27,len=3)
h<-2^(seq(0,3))
out<-classif.kernel.fb(mlearn,glearn,h=h,par.basis=list(nbasis=vbasis))
out$prob.classification
1-out$misclassification
}

\keyword{classif}