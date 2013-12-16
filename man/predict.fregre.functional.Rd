\name{predict.functional.response}
\Rdversion{1.1}
\alias{predict.fregre.fr}
%\alias{predict.fregre.lm.func}
%\alias{predict.fregre.glm.func}
%\alias{predict.fregre.gsam.func}
%\alias{predict.fregre.gkam.func}
\title{ Predict method for functional response model}
\description{
 Computes predictions for regression between functional explanatory variables and functional response.
%\itemize{      
%\item \code{predict.fregre.lm.func},	Predict method for functional linear model of \code{\link{fregre.lm}} fits object using basis or principal component representation.
%\item \code{predict.fregre.glm.func},	Predict method for functional generalized  linear model of \code{\link{fregre.glm}} fits object using basis or principal component representation.

%\item \code{predict.fregre.gsam.func},	Predict method for functional generalized spectral additive model of \code{\link{fregre.gsam}} fits object using basis or principal component representation.
%\item \code{predict.fregre.gkam.func}, Predict method for functional generalized kernel additive model of \code{\link{fregre.gkam}} fits object using backfitting algorithm.
%}
}
\usage{
\method{predict}{fregre.fr}(object,new.fdataobj=NULL,\dots)
}
%\method{predict}{fregre.lm.func}(object,newx=NULL,type="response",\dots)
%\method{predict}{fregre.glm.func}(object,newx=NULL,type="response",\dots)
%\method{predict}{fregre.gsam.func}(object, newx = NULL, type = "response",\dots)
%\method{predict}{fregre.gkam.func}(object, newx = NULL, type = "response",\dots)
\arguments{
  %\item{object}{ \code{fregre.lm.func},  \code{fregre.glm.func},  \code{fregre.gsam.func} or \code{fregre.gkam.func} object.}
    \item{object}{ \code{fregre.fr} object.}
  \item{new.fdataobj}{ New functional explanatory data of \code{fdata} class.}
  %\item{newx}{ An optional data list in which to look for variables with which to predict. If omitted, the fitted values are used. List of new explanatory data.}
  %\item{type}{ Type of prediction (response or model term).}
  \item{\dots}{ Further arguments passed to or from other methods.}
}
%\details{ 
%These functions use the model fitting function \code{\link{lm}}, \code{\link{glm}} or \code{\link{gam}} properties.\cr
%If using functional data derived, is recommended to use a number of bases to represent beta lower than the number of bases used to represent the functional data. \cr                                                                                                                                                                
%The   item called \emph{"df"} in the \code{data} list of \code{newx} argument is a data frame  with non functional explanatory variables, as \code{\link{lm}}, \code{\link{glm}} or \code{\link{gam}}. Functional variables and functional response are introduced in the following items of \code{newx} argument.
%}
\value{Return the predicted functional data.}
%\references{
%Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
%Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
%}
\author{ Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}}
\seealso{
See Also as:  \code{\link{fregre.basis.fr}}%, \code{\link{fregre.lm.func}},  \code{\link{fregre.glm.func}},  \code{\link{fregre.gsam.func}} and \code{\link{fregre.gkam.func}}. \cr
}
\examples{
\dontrun{

# CV prediction for CandianWeather data
rtt<-c(0, 365)
basiss  <- create.bspline.basis(rtt,7)
basist  <- create.bspline.basis(rtt,9)
nam<-dimnames(CanadianWeather$dailyAv)[[2]]

# fdata class (raw data)
tt<-1:365
tempfdata<-fdata(t(CanadianWeather$dailyAv[,,1]),tt,rtt)
log10precfdata<-fdata(t(CanadianWeather$dailyAv[,,3]),tt,rtt)
rng<-range(log10precfdata) 
for (ind in 1:35){
 res1<-  fregre.basis.fr(tempfdata[-ind], log10precfdata[-ind],
 basis.s=basiss,basis.t=basist)
 pred1<-predict.fregre.fr(res1,tempfdata[ind])
 plot( log10precfdata[ind],col=1,ylim=rng,main=nam[ind])
 lines(pred1,lty=2,col=2)
 Sys.sleep(1)
}

# fd class  (smooth data)
basis.alpha  <- create.constant.basis(rtt)
basisx  <- create.bspline.basis(rtt,65)

dayfd<-Data2fd(day.5,CanadianWeather$dailyAv,basisx)
tempfd<-dayfd[,1]
log10precfd<-dayfd[,3]
for (ind in 1:35){
 res2 <-  fregre.basis.fr(tempfd[-ind], log10precfd[-ind],
 basis.s=basiss,basis.t=basist)
 pred2<-predict.fregre.fr(res2,tempfd[ind])
 plot( log10precfd[ind],col=1,ylim=range(log10precfd$coef),main=nam[ind]) 
 lines(pred2,lty=2,col=2)
 Sys.sleep(1)
}
}
}
\keyword{regression}