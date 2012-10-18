\name{depth.RT}
\Rdversion{1.1}
\alias{depth.RT}
\title{ The random Tukey depth for functional data}
\description{
 The depth.RT function implements a depth measure based on random projections using a half-space Tukey method.
}
\usage{
depth.RT(fdataobj, fdataori = fdataobj, trim = 0.25, nproj = 10, 
    proj = 1, xeps = 1e-07, draw = FALSE, ...)
}

\arguments{
 \item{fdataobj}{  A set of new curves to evaluate the depth. \code{\link{fdata}} class object. }    
 \item{fdataori}{  A set of original curves where the depth is computed.  \code{\link{fdata}} class object.}  
  \item{trim}{ The alpha of the trimming.}
  \item{nproj}{ The number projection.}
    \item{proj}{  if is a character: create the random projection using a covariance matrix by process indicated in the argument (by default, proj=1, sigma=diag(ncol(fdataobj))), else if is a  matrix of random projection provided by the user.} 
  \item{xeps}{ Accuracy. The left limit  of the empirical distribution function.}
  \item{draw}{ =TRUE, draw the curves, the sample median and trimmed mean.}
  \item{\dots}{ Further arguments passed to or from other methods.}
}

\details{ It builds random projections  and calculates the functional depth by the Tukey method combining the information of all projections.
}
\value{
    \item{median}{ Deepest curve.}
    \item{lmed}{ Index deepest element \code{median}.}
    \item{mtrim}{ code{fdata} class object with the average from the \code{(1-trim)\%} deepest curves. }
    \item{ltrim}{ Index of curves with trimmed mean \code{mtrim}. }    \item{dep}{ Depth of each curve. }
    \item{proj}{ The projection value of each point on the curves.}
}
\references{
Cuesta Albertos, J. A. and Nieto Reyes, A.  \emph{The Random Tukey Depth.}
Computational Statistics and Data Analysis (2008), Vol. 52, Issue 11, 4979{-}4988.
}
\author{
This version is created by Febrero-Bande, M., and Oviedo de la Fuente, M. modified the original version created by
Cuesta Albertos, J. A. and Nieto Reyes, A.
}
\seealso{
See Also as \code{\link{depth.RP}}, \code{\link{depth.RPD}}, \code{\link{depth.FM}} or \code{\link{depth.mode}}.
}
\examples{

#Ex: CanadianWeather data
fdataobj<-fdata(t(CanadianWeather$dailyAv[,,1]))

# Random Projections
t=1:365
out.RT=depth.RT(fdataobj,draw=TRUE)
out.RT2=depth.RT(fdataobj,trim=0.1,draw=TRUE)
out.RT3=depth.RT(fdataobj,nproj=5,draw=TRUE)
out.RT4=depth.RT(fdataobj,nproj=30,draw=TRUE)
plot(out.RT$mtrim,type="l",lwd=2)
lines(out.RT2$mtrim,col=2,lwd=2,lty=2)
lines(out.RT3$mtrim,col=3,lwd=2,lty=3)
lines(out.RT4$mtrim,col=4,lwd=2,lty=4)


}



\keyword{descriptive} 