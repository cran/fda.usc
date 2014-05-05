\name{influnce.fdata}
\Rdversion{1.1}
\alias{influence.fdata}
\title{ Functional influence measures}
\description{
	  Once estimated the functional regression model with scalar response, influence.fdata function is used to obtain the functional influence measures.
	  }
\usage{
%influence.fdata(model,...)
\method{influence}{fdata}(model,...)
}
\arguments{
  \item{model}{ \code{fregre.pc}, \code{fregre.basis} or \code{fregre.basis.cv} object.}	
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{
Identify influential observations in the functional linear model in which the predictor is functional and the response is scalar.\cr
Three statistics are introduced for measuring the influence:   Distance Cook Prediction \code{DCP}, Distance Cook Estimation \code{DCE} and
 Distance \eqn{\mbox{pe}\tilde{\mbox{n}}\mbox{a} }{} \code{DP} respectively.
}
\value{
Return:\cr 
\item{DCP}{ Cook's Distance for Prediction.}
\item{DCE}{ Cook's Distance for Estimation.}
\item{DP}{  \eqn{\mbox{Pe}\tilde{\mbox{n}}\mbox{a's} }{} Distance.}
%\deqn{\mbox{pe}\tilde{\mbox{n}}\mbox{a} }{}
}

\references{
Febrero-Bande, M., Galeano, P. and Gonzalez-Manteiga, W. (2010). \emph{Measures of influence for the functional linear model with scalar response}. Journal of Multivariate Analysis 101, 327-339.

Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{
See Also as:  \code{\link{fregre.pc}}, \code{\link{fregre.basis}},  \code{\link{influence.quan}}
}
\examples{
\dontrun{
data(tecator)
x=tecator$absorp.fdata[1:129]
y=tecator$y$Fat[1:129]

res1=fregre.pc(x,y,1:5)  
# time consuming
res.infl1=influence.fdata(res1)  
res2=fregre.basis(x,y)  
res.infl2=influence.fdata(res2)  

res<-res1
res.infl<-res.infl1
mat=cbind(y,res$fitted.values,res.infl$DCP,res.infl$DCE,res.infl$DP)
colnames(mat)=c("Resp.","Pred.","DCP","DCE","DP")
pairs(mat)
}

}
\keyword{ outliers }