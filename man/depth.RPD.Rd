\name{depth.RPD}
\Rdversion{1.1}
\alias{depth.RPD}
\title{ Provides the depth measure by random projections using derivatives}
\description{
The depth.RPD function implements a depth measure based on random projections and using several derivatives. 
}
\usage{
depth.RPD(fdataobj,nproj=50,proj=1,deriv=c(0,1),trim=0.25,
dfunc2=depth.mode,method="fmm",draw=FALSE,...)
}

\arguments{
   \item{fdataobj}{ \code{\link{fdata}} class object.}    
   \item{nproj}{ Number of random projection.}
  \item{deriv}{  Number of derivatives described in integer vector \code{deriv}. \code{=0} means no derivative.}
  \item{proj}{  if is a character: create the random projection using a covariance matrix by process indicated in the argument (by default, proj=1, sigma=diag(ncol(fdataobj))), else if is a  matrix of random projection provided by the user.} 
   \item{trim}{ The alpha of the trimming.}
  \item{dfunc2}{by default,  \code{\link{depth.mode}}.}
  \item{method}{  Type of derivative method. See   \code{\link{fdata.deriv}} for more details.}
  \item{draw}{ =TRUE, draw the curves, the sample median and the trimmed mean.}
  \item{\dots}{ Further arguments passed to or from other methods.}
}

\details{ Calculate the depth in two steps. It builds random projections for the funcion and their derivatives (indicated in the parameter \code{deriv}) without looking at the functional data . Then it applyes a depth function  (by default \code{\link{depth.mode}}) to the set of random projections  created for calculate the functional depth by the Tukey method.
}
\value{
    \item{median}{ Deepest curve.}
    \item{lmed}{ Index deepest element \code{median}.}
		\item{mtrim}{ \code{fdata} class object with the average from the \code{(1-trim)\%} deepest curves. }
    \item{ltrim}{ Index of curves with trimmed mean \code{mtrim}. }    \item{dep}{ Depth of each curve. }
    \item{dep}{ Depth of each curve.}
    \item{proj}{ The projection value of each point on the curves.}
 }
\references{
Cuevas, A., Febrero-Bande, M. and Fraiman, R. (2007). \emph{Robust estimation and classification for functional data via projection-based depth notions.} Computational Statistics 22, 3, 481{-}496.\cr 

Febrero-Bande, M., Galeano, P., and Gonzalez-Manteiga, W. (2008).  \emph{Outlier detection in functional data by depth measures with application to identify abnormal NOx levels}. Environmetrics 19, 4, 331{-}345. \cr 

J. A. Cuesta-Albertos, R. Fraiman and T. Ransford. \emph{Random projections and goodness-of-fit tests in infinite-dimensional spaces.} Boletim da Sociedade Brasileira de Matematica (2006), 37(4), 1{-}25. \cr
}
\author{
Febrero-Bande, M., and Oviedo de la Fuente, M.
}
\seealso{
See Also as \code{\link{depth.RP}}, \code{\link{depth.FM}} and \code{\link{depth.mode}}.
}
\examples{

#Ex: CanadianWeather data
fdataobj<-fdata(t(CanadianWeather$dailyAv[,,1]))
t=1:365
# Double Random Projections
out.RPD=depth.RPD(fdataobj,deriv=c(0,1),dfunc2=depth.RP,draw=TRUE)
out.RPD2=depth.RPD(fdataobj,deriv=c(0,1),dfunc2=depth.RP,trim=0.1)
out.RPD3=depth.RPD(fdataobj,deriv=c(0,1),dfunc2=depth.mode,draw=TRUE,
method="fmm")
out.RPD4=depth.RPD(fdataobj,deriv=c(0,1),dfunc2=depth.FM,draw=TRUE,
method="fmm")

plot(out.RPD$mtrim,type="l",lwd=2)
lines(out.RPD2$mtrim,col=2,lwd=2,lty=2)
lines(out.RPD3$mtrim,col=3,lwd=2,lty=3)
lines(out.RPD4$mtrim,col=4,lwd=2,lty=4)
}
\keyword{descriptive} 