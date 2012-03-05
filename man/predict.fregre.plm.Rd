\name{predict.fregre.plm}
\Rdversion{1.1}
\alias{predict.fregre.plm}
\title{ Predict method for semi-functional  linear regression model.}
\description{
 	 Computes predictions for  functional (and non functional) explanatory variables and scalar response using asymmetric kernel estimation.
}
\usage{
%predict.fregre.plm(object,newx=NULL,...)
\method{predict}{fregre.plm}(object,newx=NULL,\dots)
}
\arguments{
   \item{object}{ \code{fregre.plm} object.}
  \item{newx}{ An optional data list in which to look for variables with which to predict. If omitted, the fitted values are used. List of new explanatory data.}
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{ 
Predicts from a fitted \code{fregre.plm} object.\cr

The first item in the \code{data} list of \code{newx} argument is called \emph{"df"} and  is a data frame  with the response and non functional explanatory variables, as \code{\link{lm}}. Functional variable (\code{fdata} and \code{fd} class) is introduced in the second item in the \code{data} list of \code{newx} argument.\cr
}
\value{
Return:\cr
\item{fitted.values}{  Predited scalar response.}
}
\references{
Aneiros-Perez G. and Vieu P. (2005). \emph{Semi-functional partial linear regression}.  Statistics {&} Probability Letters, 76:1102-1110.\cr
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional data analysis.} Springer Series in Statistics, New York. \cr
Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University Press, 1994.
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:  \code{\link{fregre.plm}}. \cr
}
\examples{

data(tecator)
ind=1:129
x=tecator$absorp.fdata[ind,]
newx=tecator$absorp.fdata[-ind,]
y=tecator$y$Fat[ind]
newy=tecator$y$Fat[-ind]
dataf=as.data.frame(tecator$y[ind,])

f=Fat~Water+x
ldata=list("df"=dataf,"x"=x)
res.plm=fregre.plm(f,ldata)
ndataf=as.data.frame(tecator$y[-ind,])
nldata=list("df"=ndataf,"x"=newx)
pred.plm=predict.fregre.plm(res.plm,nldata)

# with 2nd derivative of functional data
x.fd=fdata.deriv(x,nderiv=2)
f2=Fat~Protein+x.fd
ldata2=list("df"=dataf,"x.fd"=x.fd)
res.plm2=fregre.plm(f2,ldata2)
res.plm.trim2=fregre.plm(f2,ldata2,par.CV=list(trim=0.01))
newx.fd=fdata.deriv(newx,nderiv=2)
nldata2=list("df"=ndataf,"x.fd"=newx.fd)
pred.plm2=predict.fregre.plm(res.plm2,nldata2)
pred.plm.trim2=predict.fregre.plm(res.plm.trim2,nldata2)

MEP.plm=((1/length(newy))*sum((newy-pred.plm)^2))/var(newy)
MEP.plm2=((1/length(newy))*sum((newy-pred.plm2)^2))/var(newy)
MEP.plm.trim2=((1/length(newy))*sum((newy-pred.plm.trim2)^2))/var(newy)
par(mfrow=c(1,3))
boxplot(newy-pred.plm,ylim=c(-5,5))
boxplot(newy-pred.plm2,ylim=c(-5,5))
boxplot(newy-pred.plm.trim2,ylim=c(-5,5))

}
\keyword{regression}