\name{depth.RP}
\Rdversion{1.1}
\alias{depth.RP}
\title{ Provides the depth measure using random projections for functional data}
\description{
 The depth.RP function implements a depth measure based on random projections.
}
\usage{
depth.RP(fdataobj,trim=0.25,nproj=50,xeps=0.0000001,draw=FALSE,...)
}

\arguments{
 \item{fdataobj}{ \code{\link{fdata}} class object.}  
   \item{trim}{ The alpha of the trimming.}
  \item{nproj}{ The number projection.}
  \item{xeps}{ Accuracy. The left limit  of the empirical distribution function.}
  \item{draw}{ =TRUE, draw the curves, the sample median and trimmed mean.}
  \item{\dots}{ Further arguments passed to or from other methods.}
}

\details{ It builds random projections  and calculates the functional depth by the  Tukey method combining the information of all projections..
}
\value{
    \item{median}{ Deepest curve.}
    \item{lmed}{ Index deepest element \code{median}.}
    \item{mtrim}{ code{fdata} class object with the average from the \code{(1-trim)\%} deepest curves. }
    \item{ltrim}{ Index of curves with trimmed mean \code{mtrim}. }    \item{dep}{ Depth of each curve. }
    \item{proj}{ The projection value of each point on the curves. }
}
\references{
Cuevas, A., Febrero-Bande, M. and Fraiman, R. (2007). \emph{Robust estimation and classification for functional data via projection-based depth notions.} Computational Statistics 22, 3, 481{-}496. \cr 

Febrero-Bande, M., Galeano, P., and Gonzalez-Manteiga, W. (2008).  \emph{Outlier detection in functional data by depth measures with application to identify abnormal NOx levels}. Environmetrics 19, 4, 331{-}345. \cr 

J. A. Cuesta-Albertos, R. Fraiman and T. Ransford. \emph{Random projections and goodness-of-fit tests in infinite-dimensional spaces.} Boletim da Sociedade Brasileira de Matematica (2006), 37(4), 1{-}25. 
}
\author{
Febrero-Bande, M., and Oviedo de la Fuente, M.
}
\seealso{
See Also as \code{\link{depth.RPD}}, \code{\link{depth.FM}} or \code{\link{depth.mode}}.
}
\examples{
#Ex: CanadianWeather data
fdataobj<-fdata(t(CanadianWeather$dailyAv[,,1]))

# Random Projections
t=1:365
out.RP=depth.RP(fdataobj,draw=TRUE)
out.RP2=depth.RP(fdataobj,trim=0.1,draw=TRUE)
out.RP3=depth.RP(fdataobj,nproj=5,draw=TRUE)
out.RP4=depth.RP(fdataobj,nproj=30,draw=TRUE)
plot(out.RP$mtrim,type="l",lwd=2)
lines(out.RP2$mtrim,col=2,lwd=2,lty=2)
lines(out.RP3$mtrim,col=3,lwd=2,lty=3)
lines(out.RP4$mtrim,col=4,lwd=2,lty=4)
}
\keyword{descriptive} 