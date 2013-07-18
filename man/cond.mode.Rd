\name{cond.mode}
\Rdversion{1.1}
\alias{cond.mode}
\title{ Conditional mode}
\description{
 Computes the mode for conditional distribution function.
}
\usage{
 cond.mode(Fc,method="monoH.FC")
}
\arguments{
  \item{Fc}{ Object estimated by  \code{cond.F} function.}
  \item{method}{ Specifies the type of spline to be used. Possible values are \emph{"diff"}, \emph{"fmm"}, \emph{"natural"}, \emph{"periodic"} and \emph{"monoH.FC"}.}
}
\details{ The conditional mode is calculated as the maximum argument of the derivative of the conditional distribution function.}
\value{
Return the mode for conditional distribution function.
 \item{mode.cond}{ Conditional mode.}
 \item{ind.mode.cond}{ Index of Conditional mode.} 
 }
\references{
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional data analysis.}
Springer Series in Statistics, New York.
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{
See Also as:  \code{\link{cond.F}},  \code{\link{cond.quantile}}  and  \link[stats]{splinefun} . 
}
\examples{

n= 500
t= seq(0,1,len=101)
beta = t*sin(2*pi*t)^2
x = matrix(NA, ncol=101, nrow=n)
y=numeric(n)
x0<-rproc2fdata(n,seq(0,1,len=101),sigma="wiener")
x1<-rproc2fdata(n,seq(0,1,len=101),sigma=0.1)
x<-x0*3+x1
fbeta = fdata(beta,t)
y<-inprod.fdata(x,fbeta)+rnorm(n,sd=0.1)

prx=x[1:100];pry=y[1:100]
ind=101;ind2=101:110
pr0=x[ind];pr10=x[ind2]
ndist=161
gridy=seq(-1.598069,1.598069, len=ndist)
## Conditional Function
I=5
## Time consuming
# res = cond.F(pr10[I], gridy, prx, pry, h=1)

# mcond=cond.mode(res)
# mcond2=cond.mode(res,method="diff")



}
\keyword{ distribution }