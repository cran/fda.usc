\name{Descriptive}
\Rdversion{1.1}
\alias{Descriptive}
\alias{func.mean}
\alias{func.trim.FM}
\alias{func.trim.mode}
\alias{func.trim.RP}
\alias{func.trim.RT}
\alias{func.trim.RPD}   
\alias{func.trim.HD}    
\alias{func.trim.PD}    
\alias{func.trim.SD}    
\alias{func.trim.MhD}   
\alias{func.med.FM}
\alias{func.med.mode}
\alias{func.med.RP}
\alias{func.med.RT}
\alias{func.med.RPD}
\alias{func.med.HD}
\alias{func.med.PD}
\alias{func.med.SD}
\alias{func.med.MhD}
\alias{func.var}
\alias{func.trimvar.FM}
\alias{func.trimvar.mode}
\alias{func.trimvar.RP}
\alias{func.trimvar.RT}
\alias{func.trimvar.RPD}


\title{ Descriptive measures for functional and multivariate data.}
\description{ Central and dispersion measures for functional and multivariate data.}
\usage{
func.mean(fdataobj)
func.var(fdataobj)
func.trim.FM(fdataobj,...)
func.trim.mode(fdataobj,...)
func.trim.RP(fdataobj,...)
func.trim.RT(fdataobj,...)
func.trim.RPD(fdataobj,...)
func.trim.HD(fdataobj,...) 
func.trim.PD(fdataobj,...) 
func.trim.SD(fdataobj,...) 
func.trim.MhD(fdataobj,...) 
func.med.FM(fdataobj,...)
func.med.mode(fdataobj,...)
func.med.RP(fdataobj,...)
func.med.RT(fdataobj,...)
func.med.RPD(fdataobj,...)
func.med.HD(fdataobj,...)
func.med.PD(fdataobj,...)
func.med.SD(fdataobj,...)
func.med.MhD(fdataobj,...)
func.trimvar.FM(fdataobj,...)
func.trimvar.mode(fdataobj,...)
func.trimvar.RP(fdataobj,...)
func.trimvar.RT(fdataobj,...)
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
\tab \code{\link{func.trim.RT}} Returns the average from the \code{(1-trim)}\% deepest curves following RT criteria. \cr
\tab \code{\link{func.trim.RPD}} Returns the average from the \code{(1-trim)}\% deepest curves following RPD criteria. \cr
\tab \code{\link{func.trim.HD}} Returns the average from the \code{(1-trim)}\% deepest curves following HD criteria. \cr         
\tab \code{\link{func.trim.PD}} Returns the average from the \code{(1-trim)}\% deepest curves following PD criteria. \cr         
\tab \code{\link{func.trim.SD}} Returns the average from the \code{(1-trim)}\% deepest curves following SD criteria. \cr        
\tab \code{\link{func.trim.MhD}} Returns the average from the \code{(1-trim)}\% deepest curves following MhD criteria. \cr       
\tab \code{\link{func.med.FM}} Returns the deepest curve following FM criteria. \cr
\tab \code{\link{func.med.mode}} Returns the deepest curve following mode criteria. \cr
\tab \code{\link{func.med.mode}} Returns the deepest curve following RP criteria. \cr
\tab \code{\link{func.med.RP}} Returns the deepest curve following RP criteria. \cr
\tab \code{\link{func.med.RPD}} Returns the deepest curve following RPD criteria. \cr
\tab \code{\link{func.med.HD}} Returns the deepest curve following HD criteria. \cr        
\tab \code{\link{func.med.PD}} Returns the deepest curve following PD criteria. \cr   
\tab \code{\link{func.med.SD}} Returns the deepest curve following SD criteria. \cr      
\tab \code{\link{func.med.MhD}} Returns the deepest curve following MhD criteria. \cr        
\tab \code{\link{func.trimvar.FM}}  Returns the marginal variance from the deepest curves followinng FM criteria. \cr
\tab \code{\link{func.trimvar.mode}}  Returns the marginal variance from the deepest curves followinng mode criteria. \cr
\tab \code{\link{func.trimvar.RP}}  Returns the marginal variance from the deepest curves followinng RP criteria. \cr
\tab \code{\link{func.trimvar.RP}}  Returns the marginal variance from the deepest curves followinng RT criteria. \cr
\tab \code{\link{func.trimvar.RPD}}  Returns the marginal variance from the deepest curves followinng RPD criteria. \cr
}
}
\references{
Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
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
# a10<-func.med.RT(fdataobj)

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
# b10<-func.trimvar.RT(fdataobj)
# dev.new()
# par(mfrow=c(1,2))
# plot(c(b1,b2,b3,b4,b5),ylim=c(0,79),main="Measures of dispersion I")
# plot(c(b1,b6,b7,b8,b9,10),ylim=c(0,79),main="Measures of dispersion II")

}
\keyword{descriptive} 