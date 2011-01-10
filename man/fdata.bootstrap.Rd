\name{fdata.bootstrap}
\Rdversion{1.1}
\alias{fdata.bootstrap}
\title{ Bootstrap samples of a functional statistic}
\description{
 fdata.bootstrap  provides  bootstrap samples for  functional data.
}
\usage{
fdata.bootstrap(fdataobj,statistic=func.mean,alpha=0.05,nb=200,
smo=0.0,draw=FALSE,draw.control=NULL,...)
}
\arguments{
	\item{fdataobj}{ \code{fdata} class object.}       
  \item{statistic}{ Sample statistic. It must be a function that returns an object of class \code{fdata}. By default, it uses sample mean \code{\link{func.mean}}. 
  See \code{\link{Descriptive}} for other statistics.}
  \item{alpha}{ Significance value.}
  \item{nb}{ Number of bootstrap resamples.}
  \item{smo}{ The smoothing parameter for the bootstrap samples as a proportion of the sample variance matrix.}
  \item{draw}{=TRUE, plot the bootstrap samples and the statistic.}
  \item{draw.control}{list that it specifies the \code{col}, \code{lty} and \code{lwd} for objects: \code{fdataobj}, \code{statistic}, \code{IN} and \code{OUT}.}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{
  \item{statistic}{  \code{fdata} class object with the statistic estimate from \code{nb} bootstrap samples.}
  \item{dband}{ Bootstrap estimate of \code{(1-alpha)\%} distance.}
  \item{rep.dist}{ Distance from every replicate.}
  \item{resamples}{  \code{fdata} class object with the bootstrap resamples.}
%  \item{center}{  \code{fdata} class object with center.}
  \item{fdataobj}{  \code{fdata} class object.}  
}
\references{
Cuevas A., Febrero-Bande, M. and Fraiman, R. (2007). \emph{Robust estimation and classification for functional data via projection-based depth notions.} Computational Statistics 22, 3: 481{-}496. \cr

Cuevas A., Febrero-Bande, M., Fraiman R. 2006.  \emph{On the use of bootstrap for estimating functions with functional data.} Computational Statistics and Data Analysis 51: 1063{-}1074.
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as \code{\link{Descriptive}}
}

\examples{

data(tecator)
absorp<-tecator$absorp.fdata
#Bootstrap for Trimmed Mean with depth mode
out.boot=fdata.bootstrap(absorp,statistic=func.trim.FM,nb=20,draw=TRUE)
names(out.boot)
#No run
#Bootstrap for Median with with depth mode
#control=list("col"=c("grey","blue","cyan"),"lty"=c(2,1,1),"lwd"=c(1,3,1))
#out.boot=fdata.bootstrap(absorp,statistic=func.med.mode,
#draw=TRUE,draw.control=control)
}
\keyword{datagen}