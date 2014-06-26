\name{dcor.xy}
\alias{dcor.test}
\alias{dcor.dist}       
\alias{bcdcor.dist}
\alias{dcor.xy}
\title{ Distance Correlation Statistic and t-Test}
\description{
 Distance correlation t-test of multivariate and functional independence (wrapper functions of energy package).}
\usage{
dcor.xy(x,y,test=TRUE,metric.x,metric.y,par.metric.x,par.metric.y,n)
dcor.test(D1,D2,n)
bcdcor.dist(D1,D2,n)
dcor.dist(D1, D2)
}
\arguments{
  \item{x}{ data (fdata, matrix or data.frame class) of first sample.} 
  \item{y}{ data (fdata, matrix or data.frame class) of second sample.} 
 \item{test}{ if TRUE, compute bias corrected distance correlation statistic and the corresponding t-test, else compute distance correlation statistic.}
  \item{metric.x,metric.y}{Name of metric or semi-metric function used for compute the 
  distances of \code{x} and \code{y} object respectively. By default, \code{\link{metric.lp}} for functional data and \code{\link{metric.dist}} for multivariate data.}
  \item{par.metric.x,par.metric.y}{ List of parameters for the corresponding metric function.}
  \item{n}{ The sample size used in bias corrected version of distance correlation, by default is the number of rows of \code{x}.} 
  \item{D1}{ Distances of first sample.}
  \item{D2}{ Distances of second sample.}  
  }
\details{
These wrapper functions extend the functions of the \code{energy} package for multivariate data to functional data. 
Distance correlation is a measure of dependence between random vectors introduced by Szekely, Rizzo, and Bakirov (2007).

 \code{dcor.xy}  performs a nonparametric t-test of 
 multivariate or functional  independence in high dimension. The distribution of
 the test statistic is approximately Student t with \eqn{n(n-3)/2-1}
 degrees of freedom and for \eqn{n \geq 10} the statistic is approximately   
 distributed as standard normal.  Wrapper function of \code{energy:::dcor.ttest}.   The t statistic is a transformation of a bias corrected version of distance correlation (see SR 2013 for details). Large values (upper tail) of the t statistic are significant.
 
  \code{dcor.test} similar to  \code{dcor.xy}  but only for distance matrix.
    
 \code{dcor.dist}  compute distance correlation statistic.  Wrapper function of \code{energy:::dcor}  but only for distance matrix.
 
 \code{bcdcor.dist} compute bias corrected distance correlation statistic.  Wrapper function of \code{energy:::bcdcor}  but only for distance matrix.
}
\value{

\code{dcor.test} returns a list with class \code{htest} containing
   \item{     method}{ description of test}
   \item{  statistic}{ observed value of the test statistic}
   \item{  parameter}{ degrees of freedom}
   \item{   estimate}{ bias corrected distance correlation \code{bcdcor(x,y)}}
   \item{    p.value}{ p-value of the t-test}
   \item{  data.name}{ description of data}

\code{dcor.xy} returns the previoius list with class \code{htest} and
   \item{  D1}{ the distance matrix of \code{x}}   
   \item{  D2}{ the distance matrix of \code{y}}   
   
 \code{dcor.dist} returns the distance correlation statistic.
 
 \code{bcdcor.dist} returns the bias corrected distance correlation statistic.
   
}
\seealso{
 \code{\link{metric.lp}} amd \code{\link{metric.dist}}.
}

\references{
 Szekely, G.J. and Rizzo, M.L. (2013). The distance correlation t-test of  independence in high dimension. \emph{Journal of Multivariate Analysis},  Volume 117, pp. 193-213. \cr
 \url{http://dx.doi.org/10.1016/j.jmva.2013.02.012}

Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007), 
 Measuring and Testing Dependence by Correlation of Distances, 
 \emph{Annals of Statistics}, Vol. 35 No. 6, pp. 2769-2794.
 \cr \url{http://dx.doi.org/10.1214/009053607000000505}
} 
\author{ 
Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es} and
Manuel Febrero Bande
}
\examples{
x<-rproc2fdata(100,1:50)
y<-rproc2fdata(100,1:50)
dcor.xy(x, y,test=TRUE)
dx <- metric.lp(x)
dy <- metric.lp(y)
dcor.test(dx, dy)
bcdcor.dist(dx, dy)
dcor.xy(x, y,test=FALSE)
dcor.dist(dx, dy)
}   
\keyword{ htest }
\keyword{ multivariate }
\keyword{ nonparametric }
\concept{ independence }
\concept{ multivariate }
\concept{ functional data }
\concept{ distance correlation }

