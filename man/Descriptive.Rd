\name{Descriptive}
\Rdversion{1.1}
\alias{Descriptive}
\alias{func.mean}
\alias{func.trim.FM}
\alias{func.trim.mode}
\alias{func.trim.RP}
\alias{func.trim.RPD}
\alias{func.med.FM}
\alias{func.med.mode}
\alias{func.med.RP}
\alias{func.med.RPD}
\alias{func.var}
\alias{func.trimvar.FM}
\alias{func.trimvar.mode}
\alias{func.trimvar.RP}
\alias{func.trimvar.RPD}


\title{ Descriptive measures for functional data.}
\description{ Central and dispersion measures for functional data.}
\usage{
func.mean(fdataobj)
func.var(fdataobj)
func.trim.FM(fdataobj,...)
func.trim.mode(fdataobj,...)
func.trim.RP(fdataobj,...)
func.trim.RPD(fdataobj,...)
func.med.FM(fdataobj,...)
func.med.mode(fdataobj,...)
func.med.RP(fdataobj,...)
func.med.RPD(fdataobj,...)
func.trimvar.FM(fdataobj,...)
func.trimvar.mode(fdataobj,...)
func.trimvar.RP(fdataobj,...)
func.trimvar.RPD(fdataobj,...)

}
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}   
  \item{\dots}{ Further arguments passed to or from other methods.\cr
  			If the argument \code{p} is passed, it used \code{\link{metric.lp}} function, by default \code{p=2}.\cr
  			If the argument \code{trim} (alpha of the trimming) is passed, it used \code{\link{metric.lp}} function.\cr
  			If the argument \code{deriv} (number of derivatives to use) is passed. This parameter is used  in \code{\link{depth.RPD}} function, by default it uses \code{deriv =(0,1)}.}
 }   
\value{
\tabular{ll}{
\tab \code{\link{func.mean}} gives mean curve. \cr
\tab \code{\link{func.var}} gives variance curve. \cr
\tab \code{\link{func.trim.FM}} Returns the average from the \code{(1-trim)}\% deepest curves following FM criteria. \cr
\tab \code{\link{func.trim.mode}} Returns the average from the \code{(1-trim)}\% deepest curves following mode criteria. \cr
\tab \code{\link{func.trim.RP}} Returns the average from the \code{(1-trim)}\% deepest curves following RP criteria. \cr
\tab \code{\link{func.trim.RPD}} Returns the average from the \code{(1-trim)}\% deepest curves following RPD criteria. \cr
\tab \code{\link{func.med.FM}} Returns the deepest curve following FM criteria. \cr
\tab \code{\link{func.med.mode}} Returns the deepest curve following mode criteria. \cr
\tab \code{\link{func.med.RP}} Returns the deepest curve following RP criteria. \cr
\tab \code{\link{func.med.RPD}} Returns the deepest curve following RPD criteria. \cr
\tab \code{\link{func.trimvar.FM}}  Returns the marginal variance from the deepest curves followinng FM criteria. \cr
\tab \code{\link{func.trimvar.mode}}  Returns the marginal variance from the deepest curves followinng mode criteria. \cr
\tab \code{\link{func.trimvar.RP}}  Returns the marginal variance from the deepest curves followinng RP criteria. \cr
\tab \code{\link{func.trimvar.RPD}}  Returns the marginal variance from the deepest curves followinng RPD criteria. \cr

}
}
\references{
Cuevas, A., Febrero-Bande, M. and Fraiman, R. (2007). \emph{Robust estimation and classification for functional data via projection-based depth notions.} Computational Statistics 22, 3, 481{-}496. \cr

Cuevas, A., Febrero-Bande, M. and Fraiman, R.	 (2006). \emph{On the use of the bootstrap for estimating functions with functional data.} Computational Statistics and Data Analysis 51, n 2, 1063{-}107. \cr

Fraiman R, Muniz G. 2001.  \emph{Trimmed means for functional data.} Test 10: 419{-}440.\cr

Cuesta Albertos, J. A. and Nieto Reyes, A.  \emph{The Random Tukey Depth.}
Computational Statistics and Data Analysis (2008), Vol. 52, Issue 11, 4979{-}4988.\cr

Cuesta Albertos, J. A. and Nieto Reyes, A.  \emph{A Random Functional Depth.} Functional and Operational Statistics. 
Physica Verlag, (2008), 121-126. 
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\examples{

# Example with Montreal Daily Temperature (fda-package)
fdataobj<-fdata(MontrealTemp)
# Measures of central tendency
a1<-func.mean(fdataobj)
a2<-func.trim.FM(fdataobj)
a3<-func.trim.mode(fdataobj)
a4<-func.trim.RP(fdataobj)
# a5<-func.trim.RPD(fdataobj,deriv=c(0,1)) # Time-consuming
a6<-func.med.FM(fdataobj)
a7<-func.med.mode(fdataobj)
a8<-func.med.RP(fdataobj)
# a9<-func.med.RPD(fdataobj,deriv=c(0,1)) # Time-consuming

dev.new()
par(mfrow=c(1,2))
plot(c(a1,a2,a3,a4),ylim=c(-26,29),main="Central tendency: trimmed mean")
plot(c(a1,a6,a7,a8),ylim=c(-26,29),main="Central tendency: median")

## Not run
## Measures of dispersion
# b1<-func.var(fdataobj)
# b2<-func.trimvar.FM(fdataobj)
# b3<-func.trimvar.FM(fdataobj,trim=0.1)
# b4<-func.trimvar.mode(fdataobj)
# b5<-func.trimvar.mode(fdataobj,p=1)
# b6<-func.trimvar.RP(fdataobj)
# b7<-func.trimvar.RPD(fdataobj)
# b8<-func.trimvar.RPD(fdataobj)
# b9<-func.trimvar.RPD(fdataobj,deriv=c(0,1))
# dev.new()
# par(mfrow=c(1,2))
# plot(c(b1,b2,b3,b4,b5),ylim=c(0,79),main="Measures of dispersion I")
# plot(c(b1,b6,b7,b8,b9),ylim=c(0,79),main="Measures of dispersion II")

}
\keyword{descriptive} 