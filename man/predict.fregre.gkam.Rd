\name{predict.fregre.gkam}
\Rdversion{1.1}
\alias{predict.fregre.gkam}
\title{ Predict method for functional generalized kernel additive model of fregre.gkam fits object}
\description{
 	 Computes predictions for regression between functional explanatory variables and scalar response using backfitting algorithm.
}
\usage{
%predict.fregre.kgam(object, newx = NULL, type = "response", ...) 
\method{predict}{fregre.gkam}(object, newx = NULL, type = "response",\dots)
}
\arguments{
   \item{object}{ \code{fregre.gkam} object.}
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
Return: \item{predict.kgam}{ produces a vector of predictions} 
}
\references{
Febrero-Bande M. and Gonz\'{a}lez-Manteiga W. (2011). \emph{Generalized Additive Models for Functional Data}. Recent Advances in functional Data Analysis and Related Topics. Physica-Velag. \cr
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
 ab=tecator$absorp.fdata
 ab2=fdata.deriv(ab,2)

## Example: Binomial family 
 ycat=ifelse(tecator$y[,"Fat"]<15,0,1)
 yfat.cat<-ycat[ind]
 f<-yfat.cat~ab+ab2
 xlist=list("df"=data.frame(yfat.cat),"ab"=ab[ind],"ab2"=ab2[ind])

## Time-consuming 
# res=fregre.gkam(f,family=binomial(),data=xlist,control=list(maxit=2))
# newxlist=list("ab"=ab[-ind],"ab2"=ab2[-ind])

# pred=predict(res,newxlist)
# table(ycat[-ind],ifelse(pred<.5,0,1))


}
\keyword{regression}