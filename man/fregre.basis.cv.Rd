\name{fregre.basis.cv}
\Rdversion{1.1}
\alias{fregre.basis.cv}
\title{ Cross-validation Functional Regression with scalar response using basis representation. }
\description{
	Computes functional regression between functional explanatory variables and scalar response using basis representation by cross-validation method.
}
\usage{
fregre.basis.cv(fdataobj,y,basis.x=NULL,basis.b=NULL,
type.basis=NULL,lambda=0,Lfdobj=vec2Lfd(c(0,0),rtt),
type.CV=GCV.S,par.CV=list(trim=0),...)
}
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.} 
	\item{y}{ Scalar response with length \code{n}.}
	\item{basis.x}{ Basis for functional explanatory data \code{fdataobj}.}
	\item{basis.b}{  Basis for functional beta parameter.}
	 \item{type.basis}{ A vector of character string which determines type of basis. By default \emph{"bspline"}. It is  only used when
		\code{basis.x} or \code{basis.b}  are a vector of number of basis considered.}	
  \item{lambda}{ A roughness penalty. By default, no penalty \code{lambda=0}.}	
   \item{Lfdobj}{ See \link[fda]{eval.penalty}.}     
  \item{type.CV}{ Type of cross-validation. By default generalized cross-validation \code{\link{GCV.S}} method.}  
  \item{par.CV}{ List of parameters for \code{type.CV}: \code{trim}, the alpha of the trimming and \code{draw}.}	  
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{ If \code{basis = NULL} creates  bspline basis.\cr 
If the functional covariate \code{fdataobj} is in a format \code{raw data}, such as matrix or data.frame, creates an object of class \code{fdata} with default attributes, see \code{\link{fdata}}.
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
\item{y}{ Scalar response.}
\item{fdataobj}{ Functional explanatory data of class \code{fdata}.}
\item{lambda.opt}{  \code{lambda} value that minimizes CV or GCV method.}
\item{gcv.opt}{ Minimum value of CV or GCV method. If \code{lambda=0}, \code{gcv.opt=NA}.}
\item{basis.x.opt}{  Basis used for  functional explanatory data estimation  \code{fdata}.}
\item{basis.b}{  Basis used for for functional \code{beta} parameter estimation.}
\item{lm}{  Return \code{lm} object.}
}
\references{
Ramsay, James O. and Silverman, Bernard W. (2006), \emph{Functional Data Analysis}, 2nd ed., Springer, New York. \cr
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:   \code{\link{fregre.basis}}, \code{\link{summary.fregre.fd}} and \code{\link{predict.fregre.fd}} .\cr
Alternative method: \code{\link{fregre.pc.cv}} and   \code{\link{fregre.np.cv}}.
}
\examples{

data(tecator)
x<-tecator$absorp.fdata
y=tecator$y[,1]
tt=x[["argvals"]]
a=seq(21,77,by=8)
b=seq(5,17,by=4)
l=2^(-4:15)
res1=fregre.basis.cv(x,y)
summary(res1)
res2=fregre.basis.cv(x,y,a,b,lambda=l,type.CV=GCV.S,par.CV=list(trim=0.05))
ind=c(7,8,9,13,14)
rbind(c(res1[ind],res1$basis.x.opt$type,res1$basis.x.opt$nbasis,
res1$basis.x.opt$type,res1$basis.b.opt$nbasis),
c(res2[ind],res2$basis.x.opt$type,res2$basis.x.opt$nbasis,res2$basis.x.opt$type
,res2$basis.b.opt$nbasis))

}
\keyword{regression}