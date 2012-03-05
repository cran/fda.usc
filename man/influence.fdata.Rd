\name{influnce.fdata}
\Rdversion{1.1}
\alias{influence.fdata}
\title{ Functional influence measures}
\description{
	  Once estimated the functional regression model with scalar response, influence.fd function is used to obtain the functional influence measures.
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
\item{DCE}{ Cook's Distance for  Estimation.}
\item{DP}{  \eqn{\mbox{Pe}\tilde{\mbox{n}}\mbox{a's} }{} Distance.}
%\deqn{\mbox{pe}\tilde{\mbox{n}}\mbox{a} }{}
}

\references{
Febrero-Bande, M., Galeano, P. and Gonzalez-Manteiga, W. (2010). \emph{Measures of influence for the functional linear model with scalar response}. Journal of Multivariate Analysis 101, 327-339.
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:  \code{\link{fregre.pc}}, \code{\link{fregre.basis}},  \code{\link{influence.quan}}
}
\examples{

data(tecator)
x=tecator$absorp.fdata
y=tecator$y$Fat
#NO RUN
#res1=fregre.pc(x,y,1:3)  
#res.infl1=influence.fdata(res1)  
#res2=fregre.basis(x,y)  
#res.infl2=influence.fdata(res2)  

#res<-res1
#res.infl<-res.infl1
#mat=cbind(y,res$fitted.values,res.infl$DCP,res.infl$DCE,res.infl$DP)
#colnames(mat)=c("Resp.","Pred.","DCP","DCE","DP")
#pairs(mat)


}
\keyword{ outliers }