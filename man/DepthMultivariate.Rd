\name{Depth for a multivariate dataset}
\Rdversion{1.1}
\alias{Depth.Multivariate}
\alias{mdepth.HS}
\alias{mdepth.MhD}
\alias{mdepth.SD}
\alias{mdepth.LD}
\alias{mdepth.TD}
\alias{mdepth.RP}
%\alias{mdepth.MB} 
\title{ Provides the depth measure for multivariate data}
\description{
 Compute measure of centrality of the multivariate data. Type of depth function:  simplicial depth (SD), Mahalanobis depth (MhD), Random Half--Space depth (HS), random projection depth (RP) and Likelihood Depth (LD).
 \itemize{
   \item  The \code{\link{mdepth.SD}} function provides the simplicial depth measure for bivariate data.  
   \item The \code{\link{mdepth.MhD}}function implements a Mahalanobis depth measure.
   \item The \code{\link{mdepth.RP}} function provides the depth measure using random projections for multivariate data.
   \item  The \code{\link{mdepth.LD}} function provides the Likelihood depth measure for multivariate data.  
   \item  The \code{\link{mdepth.TD}} function provides the Tukey depth measure for multivariate data.     
%   \item  The \code{\link{mdepth.MB}} function provides the Modified Band depth.  
      }
}

\usage{
mdepth.SD(x,  xx = NULL, scale=FALSE) 
mdepth.HS(x, xx=x,proj=50,scale=FALSE,xeps=1e-15,random=FALSE)
mdepth.MhD(x,xx=x,scale=FALSE)
mdepth.RP(x,  xx = x, proj = 50, scale=FALSE) 
mdepth.TD(x,xx=x,xeps=1e-15,scale=FALSE)
mdepth.LD(x,xx=x,metric=metric.dist,h=NULL,scale=FALSE,...)  
}
\arguments{
 \item{x}{is a set of points, a d-column matrix.}    
 \item{xx}{  is a d-dimension multivariate sample, a d-column matrix.}  
 \item{proj}{ are the directions for random projections, by default 500 random projections generated from a scaled \code{runif(500,-1,1)}.}  
 \item{scale}{ =TRUE, scale the depth, see \link[base]{scale}.}   
\item{metric}{ Metric function, by default \code{\link{metric.dist}}. Distance matrix between \code{x} and \code{xx} is computed.} 
 \item{xeps}{ Accuracy. The left limit  of the empirical distribution function.}
  \item{random}{=TRUE for random projections. =FALSE for deterministic projections.}
  \item{h}{ Bandwidth, \code{h>0}. Default argument values are provided as the 15\%--quantile of the distance between \code{x} and \code{xx}.}
   \item{\dots}{ Further arguments passed to or from other methods.}  
} 

\details{
Type of depth measures,
 \itemize{
\item The \code{\link{mdepth.SD}} calculates the simplicial depth (HD) of the points in \code{x} w.r.t. \code{xx}. 
\item The \code{\link{mdepth.HS}} function calculates the random half--space depth (HS) of the points in \code{x} w.r.t. \code{xx} based on random projections \code{proj}.
 \item The \code{\link{mdepth.MhD}} function  calculates the Mahalanobis depth (MhD) of the points in \code{x} w.r.t. \code{xx}.
 \item The  \code{\link{mdepth.RP}} calculates the random projection depth (RP) of the points in \code{x} w.r.t. \code{xx} based on random projections \code{proj}.
 \item The \code{\link{mdepth.LD}} calculates the Likelihood depth (LD) of the points in \code{x} w.r.t. \code{xx}.   
% \item The \code{\link{mdepth.MB}} calculates the Modified Band depth (MB) of the points in \code{x} w.r.t. \code{xx}.   
 }
}

\value{
     \item{lmed}{ Index of deepest element \code{median} of \code{xx}.}
    \item{ltrim}{ Index of set of points  \code{x} with trimmed mean \code{mtrim}. }      
    \item{dep}{ Depth of each point \code{x} w.r.t. \code{xx}.}
    \item{proj}{ The projection value of each point on set of points. }
} 

\references{
Liu, R. Y., Parelius, J. M., and Singh, K. (1999). Multivariate analysis by data depth: descriptive statistics, graphics and inference,(with discussion and a rejoinder by Liu and Singh). \emph{The Annals of Statistics}, 27(3), 783-858.
%Li, J., P.C., Cuesta--Albertos, J.A. and Liu, R. (2012). DD--Classifier: Nonparametric Classification Procedure Based on DD--plot.  \emph{Journal of the American Statistical Association} , 107, 737--753. 
%Mahalanobis, P. (1936). \emph{On the generalized distance in statistics}. Proceedings of the National Academy India, 12, 49-55.
%Liu, R.Y. (1990). {On a notion of data depth based on random simplices}. \emph{Annals of Statistics}, 18, 405-414.
%Zuo, Y.J. (2003). Projection based depth functions and associated medians. Annals of Statistics, 31, 1460_1490.
%Lopez--Pintado, S. and Romo, J. (2009). On the concept of depth for functional data.\emph{Journal of the American Statistical Association}, 104(486):718--734.  
%Tukey, J. (1975). \emph{Mathematics and picturing data.} Proceedings of the 1975 International
%Congress of Mathematics, 2, 523--531.  
}
\author{
 \code{\link{mdepth.RP}}, \code{\link{mdepth.MhD}} and \code{\link{mdepth.HS}} are  versions created by Manuel Febrero Bande and Manuel Oviedo de la Fuente of the original version created by Jun Li, Juan A. Cuesta Albertos and Regina Y. Liu for polynomial classifier.
}

\seealso{
Functional depth functions: \code{\link{depth.FM}}, \code{\link{depth.mode}},  \code{\link{depth.RP}},  \code{\link{depth.RPD}} and \code{\link{depth.RT}}.
}
\examples{
data(iris)
group<-iris[,5]
x<-iris[,1:2]
                                  
MhD<-mdepth.MhD(x)
PD<-mdepth.RP(x)
HD<-mdepth.HS(x)
SD<-mdepth.SD(x)

x.setosa<-x[group=="setosa",]
x.versicolor<-x[group=="versicolor",] 
x.virginica<-x[group=="virginica",]
d1<-mdepth.SD(x,x.setosa)$dep
d2<-mdepth.SD(x,x.versicolor)$dep
d3<-mdepth.SD(x,x.virginica)$dep


}
\keyword{descriptive} 