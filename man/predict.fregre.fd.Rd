\name{predict.fregre.fd}
\Rdversion{1.1}
\alias{predict.fregre.fd}
\title{ Predict method for functional linear model (fregre.fd class)}
\description{
 	 Computes predictions for  regression between functional explanatory variables and scalar response using: basis representation, Principal Components Analysis, 
 	 Partial least squares or nonparametric kernel estimation.
}
\usage{
predict.fregre.fd(object,new.fdataobj=NULL,...)
%\method{predict}{fregre.fd}(object,new.fdataobj=NULL,\dots)
}
\arguments{
   \item{object}{ \code{fregre.fd} object.}
  \item{new.fdataobj}{ New functional explanatory data of \code{fdata} class.}
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{ 
Predicts from a fitted \code{fregre.basis} object,see  \code{\link{fregre.basis}} or \code{\link{fregre.basis.cv}}\cr
Predicts from a fitted \code{fregre.pc} object,see  \code{\link{fregre.pc}} or \code{\link{fregre.pc.cv}}\cr
Predicts from a fitted \code{fregre.pls} object,see  \code{\link{fregre.pls}} or \code{\link{fregre.pls.cv}}\cr
Predicts from a fitted \code{fregre.np} object, see  \code{\link{fregre.np}} or \code{\link{fregre.np.cv}}.
}
\value{
Return:\cr
\item{fitted.values}{  Predited scalar response.}
}
\references{
Cai TT, Hall P. 2006. \emph{Prediction in functional linear regression}. Annals of Statistics 34: 2159{-}2179.\cr 

Cardot H, Ferraty F, Sarda P. 1999. \emph{Functional linear model}. Statistics and Probability Letters 45: 11{-}22.\cr 

Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional data analysis.}
Springer Series in Statistics, New York. \cr

Hall P, Hosseini{-}Nasab M. 2006. \emph{On properties of functional principal components analysis}. Journal of the Royal Statistical
Society B 68: 109{-}126.\cr

Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University Press, 1994.
Ramsay, James O., and Silverman, Bernard W. (2006), \emph{ Functional Data Analysis}, 2nd ed., Springer, New York. \cr

}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:   \code{\link{fregre.basis}}, \code{\link{fregre.basis.cv}}, \code{\link{fregre.np}}, \code{\link{fregre.np.cv}}, \cr 
\code{\link{fregre.pc}}, \code{\link{fregre.pc.cv}}, \code{\link{fregre.pls}}, \code{\link{fregre.pls.cv}} \cr
 and \code{\link{summary.fregre.fd}}.\cr
}
\examples{

data(tecator)
absorp=tecator$absorp.fdata
ind=1:129
x=absorp[ind,]
y=tecator$y$Fat[ind]

newx=absorp[-ind,]
newy=matrix(tecator$y$Fat[-ind],ncol=1)

# Functional PC regression
res.pc=fregre.pc(x,y,1:6)
summary(res.pc)
pred.pc=predict.fregre.fd(res.pc,newx)

# Functional PLS regression
res.pls=fregre.pls(x,y,1:6)
summary(res.pls)
pred.pls=predict.fregre.fd(res.pls,newx)

# Functional nonparametric regression
# Functional nonparametric regression with other semimetric.
res.np=fregre.np(x,y,Ker=AKer.tri,metric=semimetric.deriv)
summary(res.np)
pred.np=predict.fregre.fd(res.np,newx)


# Functional regression with basis representation

res.basis=fregre.basis.cv(x,y)
summary(res.basis)
pred.basis=predict.fregre.fd(res.basis,newx)

#x.d=fdata.deriv(x,nbasis=19,nderiv=1)
#res.basis2=fregre.basis.cv(x.d,y)
#summary(res.basis2)
#newx.d=fdata.deriv(newx,nbasis=19,nderiv=1)
#pred.basis2=predict.fregre.fd(res.basis2,newx.d)

dev.new()
plot(pred.pc-newy)
points(pred.pls-newy,col=2,pch=2)
points(pred.np-newy,col=3,pch=3)
points(pred.basis-newy,col=4,pch=4)
sum((pred.pc-newy)^2,na.rm=TRUE)/sum((newy-mean(newy))^2,na.rm=TRUE)
sum((pred.pls-newy)^2,na.rm=TRUE)/sum((newy-mean(newy))^2,na.rm=TRUE)
sum((pred.np-newy)^2,na.rm=TRUE)/sum((newy-mean(newy))^2,na.rm=TRUE)
sum((pred.basis-newy)^2,na.rm=TRUE)/sum((newy-mean(newy))^2,na.rm=TRUE)
}
\keyword{regression}