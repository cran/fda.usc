\name{Depth}
\Rdversion{1.1}
\alias{Depth}
\title{ Provides the depth measure for functional data}
\description{
 Compute measure of centrality of the functional data. Type of depth function:  Fraiman and Muniz (FM) depth, modal depth, random project depth and double random project depth.
}
\arguments{
 \item{fdataobj}{ \code{\link{fdata}} class object.}   
}

\details{
The Fraiman and Muniz (FM)   depth \code{\link{depth.FM}} function is based on integrated an univariate depth measure along the axis x.   \cr
\cr
The modal depth \code{\link{depth.mode}} function calculates the depth of a datum accounting the number of curves in the neighbourhood. The \code{\link{depth.RP}} function calculates depth using univariate depth tools through random projection method (RP). \cr
The  \code{\link{depth.RPD}} function also involves the derivatives of each curve. 
}

\references{
Cuevas, A., Febrero-Bande, M. and Fraiman, R. (2007). \emph{Robust estimation and classification for functional data via projection-based depth notions.} Computational Statistics 22, 3, 481{-}496 \cr

Febrero-Bande, M., Galeano, P., and Gonzalez-Manteiga, W. (2008).  \emph{Outlier detection in functional data by depth measures with application to identify abnormal NOx levels}. Environmetrics 19, 4, 331{-}345 \cr

Fraiman R, Muniz G. 2001.  \emph{Trimmed means for functional data.} Test 10: 419{-}440.\cr
Cuesta Albertos, J. A. and Nieto Reyes, A.  \emph{The Random Tukey Depth.}
Computational Statistics and Data Analysis (2008), Vol. 52, Issue 11, 4979{-}4988.\cr
}
\author{
Febrero-Bande, M., and Oviedo de la Fuente, M.
}
\seealso{
See Also as 
\code{\link{depth.FM}}, \code{\link{depth.RP}}, \code{\link{depth.RPD}} and \code{\link{depth.mode}}.
}
\examples{
#Ex: CanadianWeather data
tt=1:365
fdataobj<-fdata(t(CanadianWeather$dailyAv[,,1]),tt)
# Fraiman-Muniz Depth
out.FM=depth.FM(fdataobj,trim=0.1,draw=TRUE)
#Modal Depth
out.mode=depth.mode(fdataobj,trim=0.1,draw=TRUE)
out.RP=depth.RP(fdataobj,trim=0.1,draw=TRUE)
out.RT=depth.RP(fdataobj,trim=0.1,draw=TRUE)
# Double Random Projections
out.RPD=depth.RPD(fdataobj,deriv=c(0,1),dfunc2=depth.FM,trim=0.1,draw=TRUE)
out<-c(out.FM$mtrim,out.mode$mtrim,out.RP$mtrim,out.RPD$mtrim)
plot(out)
}
\keyword{descriptive} 