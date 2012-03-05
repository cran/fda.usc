\name{cond.quantile}
\Rdversion{1.1}
\alias{cond.quantile}
\title{ Conditional quantile}
\description{
 Computes the quantile for conditional distribution function.
}
\usage{
cond.quantile(qua=0.5,fdata0,fdataobj,y,fn,a=min(y),b=max(y),
tol=10^floor(log10(max(y)-min(y))-3),iter.max=100,...)
}
\arguments{
  \item{qua}{ Quantile value, by default the median (\code{qua}=0.5).}
  \item{fdata0}{ Conditional functional explanatory data of \code{\link{fdata}} class object.}  
  \item{fdataobj}{ Functional explanatory data of \code{\link{fdata}} class object.}  
   \item{y}{ Scalar Response.}	
  \item{fn}{ Conditional distribution function.}	
  \item{a}{ Lower limit.}
  \item{b}{ Upper limit.}
  \item{tol}{ Tolerance.}
  \item{iter.max}{ Maximum iterations allowed, by default \code{100}.}
  \item{\dots}{ Further arguments passed to or from other methods.}  
}
\value{
 Return the quantile for conditional distribution function.
 }
\references{
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional data analysis.}
Springer Series in Statistics, New York.
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:  \code{\link{cond.F}} and  \code{\link{cond.mode}}. 
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
ind4=5
y0 = gridy[ind4]

#Conditional median an CI 95 conditional
med=cond.quantile(qua=0.5,fdata0=pr0,fdataobj=prx,y=pry,fn=cond.F,h=1)
lo=cond.quantile(qua=0.025,fdata0=pr0,fdataobj=prx,y=pry,fn=cond.F,h=1)
up=cond.quantile(qua=0.975,fdata0=pr0,fdataobj=prx,y=pry,fn=cond.F,h=1)
print(c(lo,med,up))
}
\keyword{ distribution }