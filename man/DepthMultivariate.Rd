\name{Depth Multivariate}
\Rdversion{1.1}
\alias{Depth.Multivariate}
\alias{depth.PD}
\alias{depth.HD}
\alias{depth.MhD}
\alias{depth.SD}
\title{ Provides the depth measure for multivariate data}
\description{
 Compute measure of centrality of the multivariate data. Type of depth function:  Mahalanobis depth (MhD), Half-space depth (HD), projection depth (PD) and simplicial depth (SD).
 \itemize{
   \item The \code{\link{depth.HD}} function implements a half-space depth measure based on random projections.
   \item The \code{\link{depth.MhD}}function implements a Mahalanobis depth measure.
   \item The \code{\link{depth.PD}} function provides the depth measure using random projections for multivariate data.
   \item  The \code{\link{depth.SD}} function provides the simplicial depth measure for bivariate data.  
   }
}

\usage{
depth.HD(x, xx = x, trim = 0.25,  proj = NULL, draw = FALSE) 
depth.MhD(x, xx = x, trim = 0.25, draw = FALSE) 
depth.PD(x,  xx = x, trim = 0.25, proj = NULL, draw = FALSE) 
depth.SD(x,  xx = NULL, trim = 0.25, draw = FALSE) 
}
\arguments{
 \item{x}{is a set of points, a d-column matrix.}    
 \item{xx}{  is a d-dimension multivariate sample, a d-column matrix.}  
 \item{trim}{ The alpha of the trimming.}
 \item{proj}{ are the directions for random projections, by default 500 random projections generated from a scaled \code{runif(500,-1,1)}.}  
 \item{draw}{ =TRUE, draw the curves, the sample median and trimmed mean.}

} 

\details{
Type of depth measures,
 \itemize{
  \item The \code{\link{depth.HD}} function calculates the half-space depth (HD) of the points in \code{x} w.r.t. \code{xx} based on random projections \code{proj}.
  \item The \code{\link{depth.MhD}} function  calculates the Mahalanobis depth (MhD) of the points in \code{x} w.r.t. \code{xx}.
  \item The  \code{\link{depth.PD}} calculates the projection depth (PD) of the points in \code{x} w.r.t. \code{xx} based on random projections \code{proj}.
  \item The \code{\link{depth.SD}} calculates the simplicial depth (HD) of the points in \code{x} w.r.t. \code{xx}. 
 }
 
Graphic options,
\itemize{
\item If \code{draw=TRUE} and 2-column matrix: draw scatterplot  of a multivariate data set specifying a color for one of points in a gray scale of it detph value, the median in red and the trimed mean in yellow. 
\item If \code{draw=TRUE} and 2<d<6  d-column matrix: draw scatterplot matrices of a multivariate data set specifying a color for one of points in a gray scale of it detph value.
\item If \code{draw=TRUE} and d>=6  d-column matrix: star plots of a multivariate data set specifying a color for one of the stars in a gray scale of it detph value.
\item If \code{draw=``pairs''}  draw pairs plot of a multivariate data set  based on its depth (in a gray scale).
\item If \code{draw=``stars''} draw star plots of a multivariate data set.  specifying a color for one of the stars in a gray scale of it detph value.
}
}

\value{
    \item{median}{ Deepest point of \code{xx}.}
    \item{lmed}{ Index of deepest element \code{median} of \code{xx}.}
    \item{mtrim}{ The average from the \code{(1-trim)\%} deepest set of points \code{xx}. }
    \item{ltrim}{ Index of set of points  \code{x} with trimmed mean \code{mtrim}. }      
    \item{dep}{ Depth of each point \code{x} w.r.t. \code{xx}.}
    \item{dep.ori}{ Depth of each point \code{xx} w.r.t. \code{xx}.}        
    \item{proj}{ The projection value of each point on set of points. }
} 

\references{
Li, J., P.C., Cuesta--Albertos, J.A. and Liu, R. DD--Classifier: Nonparametric
 Classification Procedure Based on DD--plot. 
 Journal of the American Statistical Association (2012), 107, 737--753. 

Mahalanobis, P. (1936). \emph{On the generalized distance in statistics}. Proceedings of the National
Academy India, 12, 49-55.\cr

Liu, R.Y. (1990). {On a notion of data depth based on random simplices}. Annals of Statistics,
 18, 405-414.
 
Zuo, Y.J. (2003). Projection based depth functions and associated medians. 
Annals of Statistics, 31, 1460_1490.
  
%Tukey, J. (1975). \emph{Mathematics and picturing data.} Proceedings of the 1975 International
%Congress of Mathematics, 2, 523--531.  
}
\author{
 \code{\link{depth.PD}}, \code{\link{depth.MhD}} and \code{\link{depth.HD}} are  versions created by Manuel Febrero Bande and Manuel Oviedo de la Fuente of the original version created by Jun Li, Juan A. Cuesta Albertos and Regina Y. Liu for polynomial classifier.
}

\seealso{
See Also as \code{\link{depth.FM}}, \code{\link{depth.mode}},  \code{\link{depth.RP}},  \code{\link{depth.RPD}} and \code{\link{depth.RT}}.
}
\examples{
# data(iris)
# group<-iris[,5]
# x<-iris[,1:2]
                                  
# MhD<-depth.MhD(x,trim=0.1,draw=TRUE)
# PD<-depth.PD(x,trim=0.1,draw=TRUE)
# HD<-depth.HD(x,trim=0.1,draw=TRUE)
# SD<-depth.SD(x,trim=0.1,draw=TRUE)

# x.setosa<-x[group=="setosa",]
# x.versicolor<-x[group=="versicolor",] 
# x.virginica<-x[group=="virginica",]
# d1<-depth.SD(x,x.setosa)$dep
# d2<-depth.SD(x,x.versicolor)$dep
# d3<-depth.SD(x,x.virginica)$dep
# plot(d1,d2,col=as.numeric(group),pch=as.numeric(group))
# plot(d1,d3,col=as.numeric(group),pch=as.numeric(group))
# plot(d2,d3,col=as.numeric(group),pch=as.numeric(group))

}
\keyword{descriptive} 