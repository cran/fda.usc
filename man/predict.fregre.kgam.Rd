\name{predict.fregre.kgam}
\Rdversion{1.1}
\alias{predict.fregre.kgam}
\title{ Predict method for functional generalized additive model of fregre.kgam fits object}
\description{
 	 Computes predictions for regression between functional  explanatory variables and scalar response using: basis representation.
}
\usage{
%predict.fregre.kgam(object, newx = NULL, type = "response", ...) 
\method{predict}{fregre.kgam}(object, newx = NULL, type = "response",\dots)
}
\arguments{
   \item{object}{ \code{fregre.kgam} object.}
  \item{newx}{ An optional data list in which to look for variables with which to predict. If omitted, the fitted values are used. List of new explanatory data.}
  \item{type}{ Type of prediction (response or model term).}
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{ 
Predicts from a fitted \code{fregre.kgam} object.\cr
This function uses the model fitting function \code{\link{gam}} properties.\cr
If using functional data derived, is recommended to use a number of bases to represent beta lower than the number of bases used to represent the functional data. 

The first item in the \code{data} list of \code{newx} argument is called \emph{"df"} and  is a data frame  with the response and non functional explanatory variables, as \code{\link{lm}}. Functional variables (\code{fdata} and \code{fd} class) are introduced in the following items in the \code{data} list of \code{newx} argument.\cr


}
\value{
Return:\cr
\item{predict.kgam}{ produces a vector of predictions or a matrix of predictions and bounds with column names fit, lwr, and upr if interval is set. If se.fit is TRUE, a list with the following components is returned: fit vector or matrix as above.} 
\item{se.fit}{ standard error of predicted means.} 
\item{residual.scale}{ residual standard deviations.}
\item{df}{ degrees of freedom for residual.} 
}
\references{
Febrero-Bande M. and Gonzalez-Maneiga W. (2011). \emph{Generalized Additive Models for Functional Data}. Recent Advances in functional Data Analyisis and Related Topics. Physica-Velag. \cr
McCullagh and Nelder (1989), \emph{Generalized Linear Models} 2nd ed. Chapman and Hall. \cr
Opsomer J.D. and Ruppert D.(1997). \emph{Fitting a bivariate additive model by local polynomial regression}.Annals of Statistics, \code{25}, 186-211.
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:  \code{\link{fregre.gsam}}. \cr
}
\examples{
data(tecator)
ind<-1:129
data(tecator)
x=tecator$absorp.fdata
x.d1<-fdata.deriv(x)
tt=x[["argvals"]]
dataf=as.data.frame(tecator$y)
f=Fat~x
ldata=list("df"=dataf[ind,],"x"=x[ind],"x.d1"=x.d1[ind])
res=fregre.kgam(Fat~x,data=ldata,family=gaussian())
p<-predict.fregre.kgam(res,ldata)
#res2=fregre.kgam(Fat~x+x.d1,data=ldata,family=gaussian())
#p2<-predict.fregre.kgam(res2,ldata)
# plot(p-ldata$df$Fat)
# points(p2-ldata$df$Fat,col=2)

}
\keyword{regression}