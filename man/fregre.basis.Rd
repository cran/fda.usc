\name{fregre.basis}
\Rdversion{1.1}
\alias{fregre.basis}
\title{ Functional Regression with scalar response using basis representation.}
\description{
	Computes functional regression between functional explanatory variable and scalar response using basis representation.
}
\usage{
fregre.basis(fdataobj,y,basis.x=NULL,basis.b=NULL,
lambda=0,Lfdobj=vec2Lfd(c(0,0),rtt),...)
}
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
	\item{y}{ Scalar response with length \code{n}.}
	\item{basis.x}{ Basis for functional explanatory data \code{fdataobj}.}
	\item{basis.b}{  Basis for functional beta parameter.}
	\item{lambda}{ A roughness penalty. By default, no penalty \code{lambda=0}.}	
  \item{Lfdobj}{ See \link[fda]{eval.penalty}.}   	
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{ If \code{basis=}\emph{NULL} creates the bspline basis. \cr
If the functional covariate \code{fdataobj} is in a format \emph{"raw data"}, such as matrix or data.frame, creates an object of class "fdata" with default attributes, see \code{\link{fdata}}.
}
\value{
Return:\cr
\item{call}{ The matched call.}
\item{beta.est}{ beta coefficient estimated of class \code{fd}}
\item{a.est}{ Intercept parameter estimated}
\item{fitted.values}{  Estimated scalar response.}
\item{H}{  Hat matrix.}
\item{residuals}{ \code{y} minus \code{fitted values}.}
\item{df}{ The residual degrees of freedom.}
\item{r2}{ Coefficient of determination.}
\item{sr2}{ Residual variance.}
\item{y}{ Response.}
\item{fdataobj}{ Functional explanatory data of class \code{fdata}.}
\item{basis.b}{  Basis used for beta parameter estimation.}
\item{lm}{  Return \code{lm} object }
}
\references{
Ramsay, James O., and Silverman, Bernard W. (2006), \emph{ Functional Data Analysis}, 2nd ed., Springer, New York. \cr
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:   \code{\link{fregre.basis.cv}}, \code{\link{summary.fregre.fd}} and \code{\link{predict.fregre.fd}}.\cr
Alternative method: \code{\link{fregre.pc}} and   \code{\link{fregre.np}}.
}
\examples{

# fregre.basis
data(tecator)
names(tecator)
absorp=tecator$absorp.fdata
ind=1:129
x=absorp[ind,]
y=tecator$y[ind,1]
tt=absorp[["argvals"]]


res1=fregre.basis(x,y)
summary(res1)
basis1=create.bspline.basis(rangeval=range(tt),nbasis=19)
basis2=create.bspline.basis(rangeval=range(tt),nbasis=9)
res5=fregre.basis(x,y,basis1,basis2)
summary(res5)
x.d2=fdata.deriv(x,nbasis=19,nderiv=1,method="bspline",class.out="fdata")
#error: no tiene argvals!
res7=fregre.basis(x.d2,y,basis1,basis2)
summary(res7)
}
\keyword{regression}