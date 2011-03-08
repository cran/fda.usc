\name{Kernel.integrate}
\Rdversion{1.1}
\alias{Kernel.integrate}
\alias{IKer.norm}
\alias{IKer.cos}
\alias{IKer.epa}
\alias{IKer.tri}
\alias{IKer.tri}
\alias{IKer.quar}
\alias{IKer.unif}
\title{ Integrate Smoothing Kernels.}
\description{ Represent   integrate  kernels: normal,  cosine,  triweight, quartic and uniform.
}
\usage{
Kernel.integrate(u,Ker=Ker.norm,a=-1)
IKer.norm(u)
IKer.cos(u)
IKer.epa(u)
IKer.tri(u)
IKer.quar(u)
IKer.unif(u)
}

\arguments{
  \item{u}{ data}
  \item{Ker}{ Type of Kernel. By default normal kernel.}
  \item{a}{ Lower  limit of integration.}
}
\details{  
Type of integrate kernel:
\tabular{ll}{
\tab		 Integrate Normal Kernel: \code{IKer.norm}\cr
\tab		 Integrate Cosine Kernel: \code{IKer.cos}\cr
\tab     Integrate Epanechnikov Kernel: \code{IKer.epa}\cr
\tab     Integrate Triweight Kernel: \code{IKer.tri}\cr
\tab     Integrate Quartic Kernel: \code{IKer.quar}\cr
\tab     Integrate Uniform Kernel: \code{IKer.unif}\cr
 } 
}
\value{
\item{res}{ Returns integrate kernel.} 
}
\references{
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional data analysis.}
Springer Series in Statistics, New York. \cr

Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University Press, 1994.
}\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:  \code{\link{Kernel}} and  \link[stats]{integrate}. 
}
\examples{


y=qnorm(seq(.1,.9,len=100))
d=IKer.tri(y)
e=IKer.cos(y)
e2=Kernel.integrate(u=y,Ker=Ker.cos)
e-e2
f=IKer.epa(y)
f2=Kernel.integrate(u=y,Ker=Ker.epa)
f-f2


plot(d,type="l",ylab="Integrate Kernel")
lines(e,col=2,type="l")
lines(f,col=4,type="l")


}

\keyword{ kernel } 