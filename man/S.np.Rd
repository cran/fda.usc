\name{S.np}
\Rdversion{1.1}
\alias{S.np}
\alias{S.LLR}
\alias{S.NW}
\alias{S.KNN}
\title{ Smoothing matrix by nonparametric methods.}
\description{
 Provides the smoothing matrix S for the discretization points \code{tt} by: \cr
   Nadaraya-Watson kernel estimator (S.NW) with bandwidth parameter \code{h}. \cr
   Local Linear Smoothing (S.LLR) with bandwidth parameter \code{h} .\cr
   K nearest neighbors estimator (S.KNN) with parameter \code{knn} .\cr
}
\usage{
S.LLR(tt, h, Ker = Ker.norm,w=NULL,cv=FALSE)
S.NW(tt, h, Ker = Ker.norm,w=NULL,cv=FALSE)
S.KNN(tt,h=NULL,Ker=Ker.unif,w=NULL,cv=FALSE)

}
\arguments{
  \item{tt}{ Vector of discretization points or distance matrix \code{mdist}}
  \item{h}{ Smoothing parameter or bandwidth. In S.KNN, number of k-nearest neighbors.}
  \item{Ker}{ Type of kernel used, by default normal kernel. }
   \item{w}{ Optional case weights.}
  \item{cv}{ =TRUE cross-validation is done.}
  }
\value{
If \code{S.LLR} return the smoothing matrix by  Local Linear Smoothing.
If \code{S.NW} return the smoothing matrix by Nadaraya-Watson kernel estimator.
If \code{S.KNN} return the smoothing matrix by k nearest neighbors estimator.
}
\references{
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional data analysis.}
Springer Series in Statistics, New York. \cr

Wasserman, L. \emph{All of Nonparametric Statistics}. Springer Texts in Statistics, 2006.
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as  \code{\link{S.basis}}}

\examples{

tt=1:101
S=S.LLR(tt,h=5)
S2=S.LLR(tt,h=10,Ker=Ker.tri)
S3=S.NW(tt,h=10,Ker=Ker.tri)
S4=S.KNN(tt,h=5,Ker=Ker.tri)
par(mfrow=c(2,2))
image(S)
image(S2)
image(S3)
image(S4)

}

\keyword{smooth} 