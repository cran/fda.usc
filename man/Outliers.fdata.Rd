\name{Outliers.fdata}
\Rdversion{1.1}
\alias{Outliers.fdata}
\alias{outliers.depth.pond}
\alias{outliers.depth.trim}
\alias{outliers.thres.lrt}
\alias{outliers.lrt}
\alias{quantile.outliers.trim}
\alias{quantile.outliers.pond}
\title{ Detecting outliers for functional dataset
}
\description{
Procedure for detecting funcitonal outliers.
}
\usage{
outliers.thres.lrt(fdataobj,nb=200,smo=0.05,trim=0.10,...)
outliers.lrt(fdataobj,nb=200,smo=0.05,trim=0.10,...)
outliers.depth.trim(fdataobj,nb=200,smo=0.05,trim=0.01,
dfunc=depth.mode,...)
outliers.depth.pond(fdataobj,nb=200,smo=0.05,dfunc=depth.mode,...)
%quantile.outliers.pond(fdataobj,dfunc=depth.mode,
%nb=200,smo=0.05,ns=0.01,...)
%quantile.outliers.trim(fdataobj,dfunc=depth.mode,trim=0.1,
%nb=200,smo=0.05,ns=0.01,...)
\method{quantile}{outliers.pond}(x,dfunc=depth.mode,
nb=200,smo=0.05,ns=0.01,\dots)
\method{quantile}{outliers.trim}(x,dfunc=depth.mode,trim=0.1,
nb=200,smo=0.05,ns=0.01,\dots)
}
\arguments{
  \item{fdataobj,x}{ \code{\link{fdata}} class object.}  
  \item{nb}{ The number of bootstrap samples.}
  \item{smo}{ The smoothing parameter for the bootstrap samples.}
	\item{trim}{ The alpha of the trimming.}  
	\item{dfunc}{ Type of depth measure, by default \code{depth.mode}.}  
	\item{ns}{ Significance level, by defaul 1\%. }
  \item{\dots}{ Further arguments passed to or from other methods.}
  
    
  }
  \details{ 
Outlier detection in functional data by likelihood ratio test (\code{outliers.depth.lrt}). The threshold for outlier detection is given by the 
 \code{outliers.lrt}.\cr

Outlier detection in functional data by depth measures: \cr
i.-\code{outliers.depth.pond} function weights the data according to depth.\cr
ii.-\code{outliers.depth.trim} function uses trimmed data.\cr

quantile.outliers.pond and quantile.outliers.trim functions provides the quantiles of the bootstrap samples for functional outlier detection by data ponderation amd trimmed data respectively. Bootstrap smoothing function (\code{\link{fdata.bootstrap}} with \code{nb} resamples) is applied to these weighted or trimmed data. If \code{smo=0} smoothed bootstrap is not performed.  The function returns a vector of size \code{1}x\code{nb} with quantiles of the bootstrap samples.
  }
\value{
  \item{outliers}{ Indexes of functional outlier.}
  \item{dep.out}{ Depth value of functional outlier.}      
  \item{dep.out}{ Iteration in which the  functional outlier is detected.}        
  \item{quantile}{ Threshold for outlier detection.}  
  \item{dep}{ Depth value of functional data.}
}
\references{
Cuevas A, Febrero M, Fraiman R. 2006.  \emph{On the use of bootstrap for estimating functions with functional data.} Computational Statistics and Data Analysis 51: 1063{-}1074.\cr

Febrero-Bande, M., Galeano, P., and Gonzalez-Manteiga, W. (2008).  \emph{Outlier detection in functional data by depth measures with application to identify abnormal NOx levels}. Environmetrics 19, 4, 331{-}345. 

Febrero-Bande, M., Galeano, P. and Gonzalez-Manteiga, W.	 (2007). \emph{A functional analysis of NOx levels: location and scale estimation and outlier detection}. Computational Statistics 22, 3, 411{-}427.

Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}

}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{
See Also: \code{\link{fdata.bootstrap}}, \code{\link{Depth}}.
}

\examples{
data(aemet)
nb=200
# NOT RUN
#out.trim<-outliers.depth.trim(aemet$temp,dfunc=depth.FM,nb=nb)
#plot(aemet$temp,col=1,lty=1)
#lines(aemet$temp[out.trim[[1]]],col=2)
}
\keyword{ outliers }