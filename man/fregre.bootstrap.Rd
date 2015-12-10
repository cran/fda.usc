\name{fregre.bootstrap}
\Rdversion{1.2}
\alias{fregre.bootstrap}
\title{  Bootstrap regression}
\description{Estimate the beta parameter by wild or smoothed bootstrap procedure }
\usage{
fregre.bootstrap(model, nb = 500, wild = TRUE,type.wild="golden",
                 newX = NULL, smo = 0.1, smoX = 0.05, alpha = 0.95, 
                 kmax.fix = FALSE, draw = TRUE,...)                 
}
\arguments{
  \item{model}{ \code{fregre.pc}, \code{fregre.pls} or  \code{fregre.basis} object.}	
 \item{nb}{ Number of bootstrap samples.}	
 \item{wild}{Naive or smoothed bootstrap depending of the \code{smo} and \code{smoX} parameters.}
 %=TRUE, wild bootstrap. =FALSE, smoothed bootstrap.}
  \item{type.wild}{Type of distribution of V in wild bootstrap procedure, see \code{\link{rwild}}.}
 \item{smo}{ If \eqn{>0}, smoothed bootstrap on the residuals (proportion of response variance).}	
 \item{smoX}{ If \eqn{>0}, smoothed bootstrap on the explanatory functional variable \code{fdata} (proportion of variance-covariance matrix of \code{fdata} object.}	
 \item{newX}{ A \code{fdata} class containing the values of the model covariates at which predictions are required (only for smoothed bootstrap).}
  \item{kmax.fix}{ The number of maximum components to consider in each bootstrap iteration.  =TRUE, the bootstrap procedure considers the same number of components used in the previous fitted model.
   =FALSE, the bootstrap procedure estimates the best components in each iteration.}
 \item{alpha}{ Significance level used for graphical option, \code{draw=TRUE}.}  
 \item{draw}{=TRUE, plot the bootstrap estimated beta, and (optional) the CI for the predicted response values.}
 \item{\dots}{ Further arguments passed to or from other methods.}	    
}
\details{
Estimate the  beta parameter by wild or smoothed bootstrap procedure  using principal components representation \code{\link{fregre.pc}}, Partial least squares components (PLS) representation  \code{\link{fregre.pls}} or basis representation \code{\link{fregre.basis}}.\cr
If a new curves are in \code{newX} argument the bootstrap method estimates the response using the bootstrap resamples.   

If the model exhibits heteroskedasticity, the use of wild bootstrap procedure is recommended (by default).

}
\value{
Return:\cr
 \item{model}{ \code{fregre.pc}, \code{fregre.pls} or  \code{fregre.basis} object.}  
\item{beta.boot}{functional beta estimated by the \code{nb} bootstrap regressions.}
\item{norm.boot}{norm of diferences beetween the nboot betas estimated by bootstrap and beta estimated by regression model.}
\item{coefs.boot}{matrix with the bootstrap estimated  basis coefficients.}
%\item{kn.boot}{If \code{kmax.fix=TRUE} is the index of the factors used in the \code{nb} models. Otherwise, list of length \code{nb} with index of the basis, PC or PLS factors selected in each bootstrap regression.}
\item{kn.boot}{vector or list of length \code{nb} with index of the basis, PC or PLS factors selected in each bootstrap regression.}
\item{y.pred}{predicted response values using \code{newX} covariates.}
\item{y.boot}{matrix of bootstrap predicted response values using \code{newX} covariates.}
\item{newX}{ a \code{fdata} class containing the values of the model covariates at which predictions are required (only for smoothed bootstrap).}
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

\dontrun{ 
data(tecator)
iest<-1:129
x=tecator$absorp.fdata[iest]
y=tecator$y$Fat[iest]
nb<-5
## Time-consuming
res.pc=fregre.pc(x,y,1:6)
# Fix the compontents used in the each regression
res.boot1=fregre.bootstrap(res.pc,nb=nb,wild=FALSE,kmax.fix=TRUE)
# Select the "best" compontents used in the each regression
res.boot2=fregre.bootstrap(res.pc,nb=nb,wild=FALSE,kmax.fix=FALSE) 
res.boot3=fregre.bootstrap(res.pc,nb=nb,wild=FALSE,kmax.fix=10) 
## predicted responses and bootstrap confidence interval
newx=tecator$absorp.fdata[-iest]
res.boot4=fregre.bootstrap(res.pc,nb=nb,wild=FALSE,newX=newx,draw=TRUE)
}
}
\keyword{ regression }