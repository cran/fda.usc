\name{depth.FM}
\Rdversion{1.1}
\alias{depth.FM}
\title{ Fraiman-Muniz depth measure}
\description{
  The Fraiman and Muniz (FM) depth. 
}
\usage{
depth.FM(fdataobj,fdataori=fdataobj,trim=0.25,xeps=0.00000001,draw=FALSE,...)
}

\arguments{
 \item{fdataobj}{  A set of new curves to evaluate the depth. \code{\link{fdata}} class object. }    
 \item{fdataori}{  A set of original curves where the depth is computed.  \code{\link{fdata}} class object.}     
   \item{trim}{ The alpha of the trimming.}
  \item{xeps}{ Accuracy. The left limit  of the empirical distribution function.}
  \item{draw}{ =TRUE, draw the curves, the sample median and trimmed mean.}
  \item{\dots}{ Further arguments passed to or from other methods.}
}

\details{ Depth.FM measures the conditional quantile  on all points and then integrated. The depth.FM method sorts the depths along the axis x, ie, performs a ranking of depths.
}
\value{
    \item{median}{ Deepest curve.}
    \item{lmed}{ Index deepest element \code{median}.}
		\item{mtrim}{ \code{fdata} class object with the average of the \code{(1-trim)\%} deepest curves. }
    \item{ltrim}{ Index of curves with trimmed mean \code{mtrim}. } 
     \item{dep}{ Depth of each curve. }
   }
\references{
Cuevas, A.; Febrero-Bande, M. and Fraiman, R. (2007). \emph{Robust estimation and classification for functional data via projection-based depth notions.} Computational Statistics 22, 3, 481{-}496.

Febrero, M., Galeano, P., and Gonzalez-Manteiga, W. (2008).  \emph{Outlier detection in functional data by depth measures with application to identify abnormal NOx levels}. Environmetrics 19, 4, 331{-}345.

Fraiman R, Muniz G. (2001).  \emph{Trimmed means for functional data.} Test 10: 419{-}440.

Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{
See Also as \code{\link{depth.RP}}, \code{\link{depth.RPD}} and \code{\link{depth.mode}}.
}


\examples{

#Ex: CanadianWeather data
t=1:365
fdataobj<-fdata(t(CanadianWeather$dailyAv[,,1]),t)

# Fraiman-Muniz Depth
out.FM=depth.FM(fdataobj,draw=TRUE)
out.FM2=depth.FM(fdataobj,trim=0.1,draw=TRUE)
plot(out.FM$mtrim,type="l",lwd=2)
lines(out.FM2$mtrim,col=2,lwd=2,lty=2)
}

\keyword{descriptive} 