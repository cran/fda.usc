\name{depth.mode}
\Rdversion{1.1}
\alias{depth.mode}
\title{ Provides the depth measure (mode) for functional data}
\description{
The depth.mode function implements the modal depth that select the curve most densely surrounded by others as the deepest one. By default, the distance is calculated  using \code{\link{metric.lp}} function. 
}
\usage{
depth.mode(fdataobj,fdataori=fdataobj,trim=0.25,
metric=metric.lp,h=NULL,scale=FALSE,draw=FALSE,...)
}

\arguments{
 \item{fdataobj}{  A set of new curves to evaluate the depth. \code{\link{fdata}} class object. }    
 \item{fdataori}{  A set of original curves where the depth is computed.  \code{\link{fdata}} class object.}  
 \item{trim}{ The alpha of the trimming. }
  \item{h}{ Bandwidth, \code{h>0}.}
  \item{metric}{ Metric function, by default \code{\link{metric.lp}}.}
  \item{scale}{ =TRUE, scale the matrix, see \link[base]{scale}.}
  \item{draw}{ =TRUE, draw the curves, the sample median and trimmed mean.}
  \item{\dots}{ Further arguments passed to or from other methods.}
}

\details{ The mode of a functional data will be the deepest curve.
}
\value{
    \item{median}{ Deepest curve.}
    \item{lmed}{ Index deepest element \code{median}.}
		\item{mtrim}{ code{fdata} class object with the average from the \code{(1-trim)\%} deepest curves. }
    \item{ltrim}{ Index of curves with trimmed mean \code{mtrim}. }    \item{dep}{ Depth of each curve. }
        \item{dep}{ Depth of each curve.}
    \item{dist}{ Distance matrix between curves or functional data.}
}
\references{
Cuevas, A.; Febrero-Bande, M. and Fraiman, R. (2007). \emph{Robust estimation and classification for functional data via projection-based depth notions.} Computational Statistics 22, 3, 481{-}496.  \

Febrero, M., Galeano, P., and Gonzalez-Manteiga, W. (2008).  \emph{Outlier detection in functional data by depth measures with application to identify abnormal NOx levels}. Environmetrics 19, 4, 331{-}345.  

Cuevas A, Febrero M, Fraiman R. 2006.  \emph{On the use of bootstrap for estimating functions with functional data.} Computational Statistics and Data Analysis 51: 1063{-}1074.

Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}

}
\author{
Febrero-Bande, M., and Oviedo de la Fuente, M.
}
\seealso{
See Also as \code{\link{depth.FM}}, \code{\link{depth.RP}}, and \code{\link{depth.RPD}}.
}

\examples{


#Ex: CanadianWeather data
t=1:365
fdataobj<-fdata(t(CanadianWeather$dailyAv[,,1]),t)
#Modal Depth
out.mode=depth.mode(fdataobj)
out.mode2=depth.mode(fdataobj,trim=0.10,draw=TRUE)
out.mode3=depth.mode(fdataobj,metric=semimetric.basis, draw=TRUE)
plot(out.mode$mtrim)
lines(out.mode2$mtrim,col=2)
lines(out.mode3$mtrim,col=3)

}

\keyword{descriptive} 