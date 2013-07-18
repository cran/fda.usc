\name{S.basis}
\Rdversion{1.1}
\alias{S.basis}
\title{ Smoothing matrix with roughness penalties by basis representation.}
\description{
Provides the smoothing matrix \code{S}  with roughness penalties. 
}
\usage{
S.basis(tt,basis,lambda=0,Lfdobj=vec2Lfd(c(0,0)),w=NULL,...)
}
\arguments{
   \item{tt}{ Discretization points.}
   \item{basis}{ Basis to use. See \link[fda]{create.basis}.}
   \item{lambda}{ A roughness penalty. By default, no penalty \code{lambda}=0.}  
   \item{Lfdobj}{ See \link[fda]{eval.penalty}.}
   \item{w}{ Optional case weights.}
  \item{\dots}{ Further arguments passed to or from other methods. Arguments to be passed by default to \link[fda]{create.basis}}
}
\details{
Provides the smoothing matrix S for the discretization points \code{tt} and b\code{basis} with roughness penalties. If \code{lambda=0} is not used penalty, else a basis roughness penalty matrix is caluclated using \link[fda]{getbasispenalty}.

} 
\value{
   Return the smoothing matrix \code{S}. 
}
\references{

Ramsay, James O. and Silverman, Bernard W. (2006). \emph{ Functional Data Analysis}, 2nd ed., Springer, New York. 

Wasserman, L. \emph{All of Nonparametric Statistics}. Springer Texts in Statistics, 2006.
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{
See Also as \code{\link{S.np}}}

\examples{

np=101
tt=seq(0,1,len=np)

nbasis=11
base1 <- create.bspline.basis(c(0, np), nbasis)
base2 <- create.fourier.basis(c(0, np), nbasis)

S1<-S.basis(tt,basis=base1,lambda=3)
image(S1) 
S2<-S.basis(tt,basis=base2,lambda=3)
image(S2)


}

\keyword{smooth} 