\name{fregre.lm}
\Rdversion{1.1}
\alias{fregre.lm}
\title{ Functional Regression with scalar response using basis representation.}
\description{
	Computes functional regression between functional (and non functional) explanatory variables and scalar response using basis representation.
}
\usage{
fregre.lm(formula,data,basis.x=NULL,basis.b=NULL,w,...)}
\arguments{
 \item{formula}{ an object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under \code{Details}.}
  \item{data}{ List that containing the variables in the model. }
	\item{basis.x}{ An optional list of basis for functional explanatory data estimation.}
	\item{basis.b}{ List of basis for functional beta parameter estimation.}
  \item{w}{ An optional vector of weights to be used in the fitting process.}
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{ 
The first item in the \code{data} list  is called \emph{"df"} and  is a data frame  with the response and non functional explanatory variables, as \code{\link{lm}}.\cr

Functional variables o class  \code{fdata} or \code{fd} are introduced in the following items in the \code{data} list.\cr

\code{basis.x} is a list of basis  for convert \code{fdata} class data to \code{fd} class data. Only used for \code{fdata}  class data.\cr

\code{basis.b} is a list of basis  for beta parameter. 

When using functional data derived recommend using a number of basis to represent beta lower than the number of basis used to represent the functional data. 
}
\value{
Return \code{lm} object plus:\cr
\item{basis.x}{ Basis used for \code{fdata} variables.}
\item{basis.b}{ Basis used for beta parameter estimation.}
\item{beta.l}{ List of estimated beta parameter of functional covariates.}
\item{data}{ List that containing the variables in the model.}
\item{formula}{ formula.}
}
\references{
Ramsay, James O., and Silverman, Bernard W. (2006), \emph{ Functional Data Analysis}, 2nd ed., Springer, New York. \cr
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:  \code{\link{predict.fregre.lm}} and \code{\link{summary.lm}}.\cr
Alternative method: \code{\link{fregre.pc}}, \code{\link{fregre.basis}} and   \code{\link{fregre.np}}.
}
\examples{

data(tecator)
x=tecator$absorp.fdata
y=tecator$y[,1]
tt=x[["argvals"]]
dataf=as.data.frame(tecator$y)

nbasis.x=11
nbasis.b=7
basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)
 
f=Fat~Protein+x
basis.x=list("x"=basis1)
basis.b=list("x"=basis2)
ldata=list("df"=dataf,"x"=x)
res=fregre.lm(f,ldata,basis.x=basis.x,basis.b=basis.b)
summary(res)

f2=Fat~Protein+xd
xd=fdata.deriv(x,nderiv=2,class.out='fdata',nbasis=nbasis.x)
ldata2=list("df"=dataf,"xd"=xd)
basis.x2=list("xd"=basis1)
basis.b2=list("xd"=basis2)
res2=fregre.lm(f2,ldata2,basis.x=basis.x2,basis.b=basis.b2)
summary(res2)

par(mfrow=c(2,1))
plot(res$beta.l$x,main="functional beta estimation")
plot(res2$beta.l$xd,col=2)

}
\keyword{regression}