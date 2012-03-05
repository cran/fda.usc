\name{fregre.plm}
\Rdversion{1.1}
\alias{fregre.plm}
\title{ Semi-functional linear regression with scalar response.}
\description{
	Computes functional regression between functional (and non functional) explanatory variables and scalar response using asymmetric kernel estimation.
}
\usage{
fregre.plm(formula,data,h=NULL,Ker=AKer.norm,metric=metric.lp,
type.CV = GCV.S,type.S=S.NW,par.CV=list(trim=0,draw=FALSE),
par.S=list(w=1),...)}
\arguments{
 \item{formula}{ an object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under \code{Details}.}
  \item{data}{ List that containing the variables in the model. }
	\item{Ker}{ Type of asymmetric kernel used, by default asymmetric normal kernel.}	
	\item{h}{ Smoothing parameter or bandwidth.}
  \item{metric}{ Metric function, by default \code{\link{metric.lp}}.}
  \item{type.CV}{ Type of cross-validation. By default generalized cross-validation \code{\link{GCV.S}} method.} 
  \item{type.S}{ Type of smothing matrix \code{S}. By default \code{S} is calculated by Nadaraya-Watson kernel estimator (\code{S.NW}).}   
  \item{par.CV}{ List of parameters for \code{type.CV}: \code{trim}, the alpha of the trimming and \code{draw=TRUE}.}	  
   \item{par.S}{ List of parameters for \code{type.S}: \code{w}, the weights.}	      
   \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{ 
The first item in the \code{data} list  is called \emph{"df"} and  is a data frame  with the response and non functional explanatory variables, as \code{link{lm}} . If non  functional data into the formula then \code{\link{lm}} regression is performed.\cr
Functional variable (\code{fdata} or \code{fd} class) is introduced in the second item in the \code{data} list.  If only  functional variable into the formula then \code{\link{fregre.np.cv}} is performed.\cr
 
The function estimates the value of smoothing parameter or the bandwidth \code{h} through  Generalized Cross-validation \code{GCV} criteria.
It computes the distance between curves using the \code{\link{metric.lp}}, although you can also use other metric funtions. \cr
Different asymmetric kernels can be used, see   \code{\link{Kernel.asymmetric}}.\cr
}
\value{
\item{call}{ The matched call.}
\item{fitted.values}{  Estimated scalar response.}
\item{residuals}{ \code{y} minus \code{fitted values}.}
\item{df}{ The residual degrees of freedom.}
\item{H}{  Hat matrix.}
\item{r2}{ Coefficient of determination.}
\item{sr2}{ Residual variance.}
\item{y}{ Scalar response.}
\item{fdataobj}{ Functional explanatory data.}
\item{XX}{ Non functional explanatory data.}
\item{mdist}{ Distance matrix between curves.}
\item{betah}{ beta coefficient estimated}
\item{data}{  List that containing the variables in the model.}
\item{Ker}{ Asymmetric kernel used.}
\item{h.opt}{ Value that minimizes CV or GCV method.}
\item{h}{ Smoothing parameter or bandwidth.}\item{data}{ List that containing the variables in the model.}
\item{gcv}{  GCV values.}
\item{formula}{ formula.} 
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
See Also as:  \code{\link{predict.fregre.plm}} and  \code{\link{summary.fregre.fd}}\cr
Alternative methods: \code{\link{fregre.lm}}, \code{\link{fregre.np}} and   \code{\link{fregre.np.cv}}
}
\examples{
data(tecator)
x=tecator$absorp.fdata
newx=tecator$absorp.fdata
y=tecator$y$Fat
dataf=as.data.frame(tecator$y)

f=Fat~Water+x
ldata=list("df"=dataf,"x"=x)
res.plm=fregre.plm(f,ldata)
summary(res.plm)

# with 2nd derivative of functional data
x.fd=fdata.deriv(x,nderiv=2)
f2=Fat~Water+x.fd
ldata2=list("df"=dataf,"x.fd"=x.fd)
res.plm2=fregre.plm(f2,ldata2)#,Ker=AKer.tri)
summary(res.plm)

par(mfrow=c(2,1))
plot(y,res.plm$fitted.values)
points(y,res.plm2$fitted.values,col=2)
plot(res.plm$residuals)
points(res.plm2$residuals,col=2)

# fregre.plm as fregre.np.cv
#f=Fat~+x
#res.np=fregre.plm(f,ldata)
#res.np2=fregre.np.cv(x,y)
#res.np
#res.np2

# fregre.plm as lm
f=Fat~+Water
res.lm=fregre.plm(f,ldata)
res.lm2=lm(f,dataf)
res.lm
res.lm2
}
\keyword{regression}