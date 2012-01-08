\name{fregre.bootstrap}
\Rdversion{1.2}
\alias{fregre.bootstrap}
\title{  Bootstrap regression}
\description{Estimate the beta parameter by wild or smoothed bootstrap procedure }
\usage{
fregre.bootstrap(model, nb = 500, wild = TRUE,type.wild="golden",newX = NULL,  smo = 0.1, 
    smoX = 0.05, alpha = 0.95, kmax.fix = FALSE, draw = TRUE,...)                 
}
\arguments{
  \item{model}{ \code{fregre.pc}, \code{fregre.pls} or  \code{fregre.basis} object.}	
 \item{nb}{ Number of bootstrap samples }	
 \item{wild}{=TRUE, wild bootstrap. =FALSE, smoothed bootstrap.}
  \item{type.wild}{Type of distribution of V in wild bootstrap procedure, see \code{\link{rwild}}.}
 \item{smo}{ (only for smoothed bootstrap)  Smoothing parameter as a proportion of response variance.}	
 \item{smoX}{ (only for smoothed bootstrap)  Smoothing parameter for \code{fdata} object as a proportion of variance-covariance matrix of the explanatory functional variable.}	
 \item{newX}{ (optional, only for smoothed bootstrap), new functional explanatory data of \code{fdata} class .} 
 \item{alpha}{ Significance level.}	
% \item{kmax.fix}{ The maximum number of principal comoponents,  partial least squares components  or number of basis is fixed by \code{model} object.}	 
 \item{kmax.fix}{=TRUE, the bootstrap procedure considers the same number of components used in the previous fitted model.
   =FALSE, the bootstrap procedure estimates the best components in each iteration.}
   \item{draw}{=TRUE, plot the bootstrap estimated beta, and (optional) the CI for the predicted response values.}
 \item{\dots}{ Further arguments passed to or from other methods.}	    
}
\details{
Compute the beta estimated in \code{\link{influence.fdata}} for functional regression using principal components representation \code{\link{fregre.pc}}, Partial least squares components (PLS) representation  \code{\link{fregre.pls}} or basis representation \code{\link{fregre.basis}}.\cr
If a new curves are  in newX argument the bootstrap method estimates the response using the bootstrap resamples.   

If the model exhibits heteroskedasticity, the use of wild bootstrap procedure is recommended (by default).

}
\value{
Return:\cr
\item{beta.boot}{ Functional beta estimated by bootstrap method.}
\item{norm.boot}{  norm of diferences beetween the nboot betas estimated by bootstrap and beta estimated by regression model.}
\item{coefs.boot}{ nboot X nbasis of bootstrap estimated coefficients of basis}
\item{knn.fix}{nb basis, CP or PLS factors selected in each regression.}
\item{y.pred}{ nboot X n matrix of predicted response values.}
\item{y.boot}{ nboot X n matrix of bootstrap predicted response values.}
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
See Also as: \code{\link{fregre.pc}}, \code{\link{fregre.pls}}, \code{\link{fregre.basis}}, . 
}
\examples{

data(tecator)
x=tecator$absorp.fdata[1:129]
y=tecator$y$Fat[1:129]

res.basis=fregre.basis(x,y)
newx=tecator$absorp.fdata[-c(1:129)]
newy=tecator$y$Fat[-c(1:129)]
nb<-10 # Time-consuming
res.boot1=fregre.bootstrap(res.basis,nb=nb,newX=newx,kmax.fix=TRUE)
predy=predict(res.basis,newx)
plot(predy-newy,col=2)

## predicted responses by bootstrap
# res.boot1$ypred

## With FPC and FPLS basis
# res.pc=fregre.pc(x,y,1:6)
# res.pls=fregre.pls(x,y)
# res.boot2=fregre.bootstrap(res.pls,nb=5,newX=newx,kmax.fix=TRUE)
# res.boot3=fregre.bootstrap(res.pc,nb=5,newX=newx,kmax.fix=TRUE)

}
\keyword{ regression }