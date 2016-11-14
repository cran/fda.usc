\name{predict.fregre.gls}
\title{Predictions from a functional gls object}
\usage{
\method{predict}{fregre.gls}(object,  newx = NULL, type = "response", 
          se.fit = FALSE, scale = NULL, df = df, interval = "none", 
          level = 0.95, weights = 1, pred.var = res.var/weights, 
          n.ahead = 1,  \dots)
\method{predict}{fregre.igls}(object,new.fdataobj=NULL,data,se.fit=FALSE,
scale = NULL,df=df,interval = "none", level = 0.95,weights = 1, 
pred.var = res.var/weights,\ldots)
}
\alias{predict.fregre.gls}
\alias{predict.fregre.igls}
\arguments{
  \item{object}{ \code{fregre.gls} object.}
  \item{newx,new.fdataobj}{ An optional data list in which to look for variables with which to predict. If omitted, the fitted values are used. List of new explanatory data.}
  \item{type}{ Type of prediction (response or model term).}
  \item{se.fit}{=TRUE (not default) standard error estimates are returned for each prediction.}
 \item{scale}{	Scale parameter for std.err. calculation.}
  \item{df}{Degrees of freedom for scale.}
  \item{interval}{Type of interval calculation.}
  \item{level}{Tolerance/confidence level.}
  \item{weights}{variance weights for prediction. This can be a numeric vector or a one-sided model formula. In the latter case, it is interpreted as an expression evaluated in newdata}  
  \item{pred.var}{the variance(s) for future observations to be assumed for prediction intervals. See \code{link{predict.lm}} for more details.}
    \item{data}{Data frame with the time or spatinal index}
  \item{n.ahead}{number of steps ahead at which to predict.}
  \item{\dots}{ Further arguments passed to or from other methods.}
}
% \item{object}{an object inheriting from class \code{"\link{gls}"}, representing  a functional generalized least squares fitted linear model.}
\description{
  The predictions for the functional generalized least squares fitted linear model represented by \code{object} are   obtained at the covariate values defined in \code{newx}.
}
\value{
  a vector with the predicted values.
}
\references{
 Oviedo de la Fuente, M., Febrero-Bande, M., Pilar Munoz, and Dominguez, A. Predicting seasonal influenza transmission using Functional Regression Models with Temporal Dependence. arXiv:1610.08718. \url{https://arxiv.org/abs/1610.08718}
}

\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{\code{\link{fregre.gls}}}

\examples{
\dontrun{
data(tecator)
ind<-1:129
x <-fdata.deriv(tecator$absorp.fdata,nderiv=1)
dataf=as.data.frame(tecator$y)
dataf$itime <- 1:nrow(x)
ldata=list("df"=dataf[ind,],"x"=x[ind])
newldata=list("df"=dataf[-ind,],"x"=x[-ind])
newy <- tecator$y$Fat[-ind]
ff <- Fat~x
res.gls=fregre.gls(ff,data=ldata, correlation=corAR1())
pred.gls <- predict.fregre.gls(res.gls,newldata)
par.cor <- list("cor.ARMA"=list("index"=c("itime"),"p"=1))
res.igls <- fregre.igls(ff,data=ldata,correlation=par.cor) 
pred.igls<-predict.fregre.igls(res.igls,newldata,data=newldata$df["index",drop=F])
}
}
\keyword{models}
