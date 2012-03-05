\name{predict.fregre.glm}
\Rdversion{1.1}
\alias{predict.fregre.glm}
\title{ Predict method for functional linear model of fregre.glm fits object}
\description{
 	 Computes predictions for regression between functional (and no functional) explanatory variables and scalar response using: basis representation.
}
\usage{
%predict.fregre.glm(object,newx=NULL,type="response",...)
\method{predict}{fregre.glm}(object,newx=NULL,type="response",\dots)
}
\arguments{
   \item{object}{ \code{fregre.glm} object.}
  \item{newx}{ An optional data list in which to look for variables with which to predict. If omitted, the fitted values are used. List of new explanatory data.}
  \item{type}{ Type of prediction (response or model term).}
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{ 
Predicts from a fitted \code{fregre.glm} object.\cr
This function uses the model fitting function \code{\link{glm}} properties.\cr
If using functional data derived, is recommended to use a number of bases to represent beta lower than the number of bases used to represent the functional data. 

The first item in the \code{data} list of \code{newx} argument is called \emph{"df"} and  is a data frame  with the response and non functional explanatory variables, as \code{\link{lm}}. Functional variables (\code{fdata} and \code{fd} class) are introduced in the following items in the \code{data} list of \code{newx} argument.\cr


}
\value{
Return:\cr
\item{predict.glm}{ produces a vector of predictions or a matrix of predictions and bounds with column names fit, lwr, and upr if interval is set. If se.fit is TRUE, a list with the following components is returned: fit vector or matrix as above.} 
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
ind<-1:129                                                                    
tecator$y$Fat<-ifelse(tecator$y$Fat<20,0,1)                                   
x<-tecator$absorp.fdata
tt=x[["argvals"]]                                                             
dataf=as.data.frame(tecator$y[ind,])                                                                               
x.d2=fdata.deriv(x,nderiv=2,method="fmm")                                           
ldata=list("df"=dataf,"x.d2"=x.d2[ind])  
basis.pc2=create.pc.basis(x.d2[ind],c(1:2))     
basis.x=list("x.d2"=basis.pc2)                 
                             
res.pc=fregre.glm(Fat~x.d2,ldata,family=binomial,basis.x=basis.x)     

basis1=create.fdata.basis(x.d2[ind],type.basis="fourier",l=1:11)              
basis2=create.fdata.basis(x.d2[ind],type.basis="fourier",l=1:5)               
basis.x=list("x.d2"=basis1)                                                   
basis.b=list("x.d2"=basis2)                                      
res.basis=fregre.glm(Fat~x.d2,ldata,family=binomial,basis.x=basis.x,basis.b=basis.b)     

newy<-floor(tecator$y$Fat[-ind])                                              
newx<-tecator$absorp.fdata[-ind,]                                             
newx.d2<-fdata.deriv(newx,nderiv=2,method="fmm")                                    
newdataf=as.data.frame(tecator$y[-ind,])                                      
newldata=list("df"=newdataf,"x.d2"=newx.d2)                                                                                    

pred.pc<-predict.fregre.glm(res.pc,newldata)                                      
y.pred1<-ifelse (pred.pc>0.5,1,0)                                               
table(y.pred1,newy)                                                           
                                                                              
pred.basis<-predict.fregre.glm(res.basis,newldata)                                      
y.pred2<-ifelse (pred.basis>0.5,1,0)                                               
table(y.pred2,newy)                                                                                                                       

}
\keyword{regression}