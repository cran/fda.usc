\name{influence.quan}
\Rdversion{1.1}
\alias{influence.quan}
\title{  Quantile for  influence measures }
\description{Estimate the quantile of measures of influence for each observation.
}
\usage{
\method{influence}{quan}(model,out.influ,mue.boot=500,
smo=0.1,smoX=0.05,alpha=0.95,kmax.fix=FALSE,\dots)
%influence.quan(model,out.influ,mue.boot=500,
                 %smo=0.1,smoX=0.05,alpha=0.95,kmax.fix=FALSE,...)
                 }
\arguments{
  \item{model}{ \code{fregre.pc}, \code{fregre.basis} or \code{fregre.basis.cv} object.}	
 \item{out.influ}{ \code{inflluence.fd} bject}	
 \item{mue.boot}{ Number of bootstrap samples }	
 \item{smo}{ Smoothing parameter as a proportion of response variance.}	
 \item{smoX}{  Smoothing parameter for \code{fdata} object as a proportion of variance-covariance matrix of the explanatory functional variable.}	
 \item{alpha}{ Significance level. }	
 \item{kmax.fix}{ The maximum number of principal comoponents or number of basis is fixed by \code{model} object.}	 
 \item{\dots}{ Further arguments passed to or from other methods.}	    
}
\details{
Compute the quantile of measures of influence estimated in \code{\link{influence.fdata}} for functional\cr
regression using principal components representation (\code{\link{fregre.pc}}) or basis representation\cr
(\code{\link{fregre.basis}} or \code{\link{fregre.basis.cv}}).\cr
A smoothed bootstrap method is used to estimate the quantiles of the influence measures, which allows  to point out which observations have the larger influence on estimation and prediction. 
}
\value{
Return:\cr
\item{quan.cook.for}{ Distance Cook Prediction Quantile.}
\item{quan.cook.est}{ Distance Cook Estimation Quantile.}
\item{quan.cook.\eqn{\mbox{pe}\tilde{\mbox{n}}\mbox{a} }{}}{ \eqn{\mbox{Pe}\tilde{\mbox{n}}\mbox{a} }{} Distance Quantile.}
\item{mues.est}{  Sample Cook generated.}
\item{mues.pena}{ Sample \eqn{\mbox{Pe}\tilde{\mbox{n}}\mbox{a} }{} generated.}
\item{beta.boot}{ Functional beta estimated by bootstrap method.}
}
\references{
Febrero-Bande, M., Galeano, P. and Gonzalez-Manteiga, W. (2010). \emph{Measures of influence for the functional linear model with scalar response}. Journal of Multivariate Analysis 101, 327-339.
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{
See Also as:  \code{\link{influence.fdata}}, \code{\link{fregre.basis}}, \code{\link{fregre.pc}}. 
}
\examples{
\dontrun{
data(tecator)
x=tecator$absorp.fdata
y=tecator$y$Fat
res=fregre.pc(x,y,1:6)

#time consuming
res.infl=influence.fdata(res)
resquan=influence.quan(res,res.infl,4,0.01,0.05,0.95)
plot(res.infl$betas,type="l",col=2)
lines(res$beta.est,type="l",col=3)
lines(resquan$betas.boot,type="l",col="gray")

res=fregre.basis(x,y)
res.infl=influence.fdata(res)
resquan=influence.quan(res,res.infl,mue.boot=4,kmax.fix=T)
plot(resquan$betas.boot,type="l",col=4)
lines(res.infl$betas,type="l",col=2)
lines(resquan$betas.boot,type="l",col="gray")
}
}
\keyword{ outliers }