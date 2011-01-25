\name{predict.fregre.lm}
\Rdversion{1.1}
\alias{predict.fregre.lm}
\title{ Predict method for functional linear model of fregre.lm fits object}
\description{
 	 Computes predictions for regression between functional (and no functional) explanatory variables and scalar response using: basis representation.
}
\usage{
predict.fregre.lm(object,newx=NULL,...)
%\method{predict}{fregre.lm}(object,newx=NULL,\dots)
}
\arguments{
   \item{object}{ \code{fregre.lm} object.}
  \item{newx}{ An optional data list in which to look for variables with which to predict. If omitted, the fitted values are used. List of new explanatory data.}
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{ 
Predicts from a fitted \code{fregre.lm} object.\cr
This function uses the model fitting function \code{\link{lm}} properties.\cr
If using functional data derived, is recommended to use a number of bases to represent beta lower than the number of bases used to represent the functional data. 

The first item in the \code{data} list of \code{newx} argument is called \emph{"df"} and  is a data frame  with the response and non functional explanatory variables, as \code{\link{lm}}. Functional variables (\code{fdata} and \code{fd} class) are introduced in the following items in the \code{data} list of \code{newx} argument.\cr


}
\value{
Return:\cr
\item{predict.lm}{ produces a vector of predictions or a matrix of predictions and bounds with column names fit, lwr, and upr if interval is set. If se.fit is TRUE, a list with the following components is returned: fit vector or matrix as above.} 
\item{se.fit}{ standard error of predicted means.} 
\item{residual.scale}{ residual standard deviations.}
\item{df}{ degrees of freedom for residual.} 
}
\references{
Ramsay, James O., and Silverman, Bernard W. (2006), \emph{ Functional Data Analysis}, 2nd ed., Springer, New York. \cr
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:  \code{\link{fregre.lm}}. \cr
}
\examples{

data(tecator)
dataf=as.data.frame(tecator$y)
x=tecator$absorp.fdata

ind=1:165
y=tecator$y[ind,1]
newy=tecator$y[-ind,1]
dataf=as.data.frame(tecator$y[ind,])
newdataf=as.data.frame(tecator$y[-ind,])

tt=x[["argvals"]]
x.d2=fdata.deriv(x[ind,],nderiv=2)
newx=x[-ind,]
newx.d2=fdata.deriv(newx,nderiv=2)
ldata=list("df"=dataf,"x"=x[ind,],"x.d2"=x.d2)
newldata=list("df"=newdataf,"x"=newx,"x.d2"=newx.d2)

basis1=create.bspline.basis(rangeval=range(tt),nbasis=25)
basis2=create.bspline.basis(rangeval=range(tt),nbasis=13)
basis.x=list("x"=basis1)
basis.b=list("x"=basis2)

#ex 1.
f=Fat~Protein+x
res.lm=fregre.lm(f,ldata,basis.x=basis.x,basis.b=basis.b)
summary(res.lm)
pred.lm=predict.fregre.lm(res.lm,newldata)
MEP.lm=((1/length(newy))*sum((newy-pred.lm)^2))/var(newy)
MEP.lm

# with 2nd derivative of functional data
f2=Fat~Protein+x.d2
basis.x2=list("x.d2"=basis1)
basis.b2=list("x.d2"=basis2)
res.lm2=fregre.lm(f2,ldata,basis.x=basis.x2,basis.b=basis.b2)
summary(res.lm2)
pred.lm2=predict.fregre.lm(res.lm2,newldata)
MEP.lm2=((1/length(newy))*sum((newy-pred.lm2)^2))/var(newy)
MEP.lm2

plot(newy,pred.lm)
points(newy,pred.lm2,col=2)

}
\keyword{regression}