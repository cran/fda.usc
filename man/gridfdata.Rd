\name{gridfdata, rcombfdata}
\Rdversion{1.1}
\alias{gridfdata}
\alias{rcombfdata}
\title{ Utils for generate functional data }
\description{

\code{gridfdata} generates \code{n} curves as lineal combination of the original curves  \code{fdataobj} plus a functional trend \code{mu}.

%Por ejemplo se pueden generar proyeccionnes a partir de las componentes principales funcionales.

\code{rcombfdata} generates \code{n} random linear combinations of the \code{fdataobj} curves plus a functional trend \code{mu}. The coefficients of the combinations follows a normal distribution with zero mean and standard deviation \code{sdarg}. 

}
\usage{
gridfdata(coef,fdataobj,mu)

rcombfdata(n = 10, fdataobj, mu, 
           sdarg = rep(1,nrow(fdataobj)), norm = 1)
}

\arguments{
\item{coef}{Coefficients of the combination. A matrix with number of columns equal to number of curves in \code{fdataobj}}

\item{fdataobj}{ \code{\link{fdata}} class object.}       

\item{mu}{ Functional trend, by default \code{mu}=\eqn{\mu(t)=0}. An object of class \code{\link{fdata}}.
 %If \code{mu} is a \code{fdata} class object, \code{t}\eqn{=}\code{argvals(mu)}.
 }  
 
\item{n}{ Number of  curves to be generated}       

\item{sdarg}{ Standard deviation of the coefficients.}       

\item{norm}{ Norm of the coefficients. The norm is adjusted before the transformation for \code{sdarg} is performed.} 

}

\value{
Return the functional trajectories as a \code{fdata} class object.
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{
See Also as  \code{\link{rproc2fdata}}
} 
\examples{
tt=seq(0,1,len=51)
fou3=create.fourier.basis(c(0,1),nbasis=3)
fdataobj=fdata(t(eval.basis(tt,fou3)),argvals=tt)

coef=expand.grid(0,seq(-1,1,len=11),seq(-1,1,len=11))
grid=gridfdata(coef,fdataobj)
plot(grid,lty=1)

rcomb=rcombfdata(n=51,fdataobj,mu=fdata(30*tt*(1-tt),tt))
plot(rcomb,lty=1)
}
\keyword{generation} 