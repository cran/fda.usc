\name{fdata}
\Rdversion{1.1}
\alias{fdata}
\title{ Converts raw data or other functional data classes into fdata class.
}
\description{

Create a functional data object of class \code{fdata} from (\code{matrix, data.frame, numeric, integer, fd, fds, fts} or \code{sfts}) class data.
}
\usage{
fdata(mdata,argvals=NULL,rangeval=NULL,names=NULL)
}

\arguments{
  \item{mdata}{ Matrix  of  set cases with dimension (\code{n} x \code{m}), where \code{n} is the number of curves and \code{m} are the points observed in each curve.}
  \item{argvals}{ Argvals, by default:  \code{1:m}.}
  \item{rangeval}{ Range of discretization points, by default:  range(\code{argvals}).}
  \item{names}{ (optional) list with tree components:  \code{main} an overall title,	\code{xlab} title for \code{x} axis and \code{ylab} title for \code{y} axis.}
}

\value{
Return \code{fdata} class object with: i.- \code{"data"}: matrix  of  set cases with dimension (\code{n} x \code{m}), where \code{n} is the number of curves and \code{m} are the points observed in each curve, ii.- \code{"rangeval"}: the discretizations points values, if not provided:  \code{1:m},  iii.- \code{"rangeval"}: range of the discretizations points values, by default:  range(\code{argvals}) , iv.- \code{"names"} (optional) list with  \code{main} an overall title,	\code{xlab} title for \code{x} axis and \code{ylab} title for \code{y} axis.
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as \code{\link{plot.fdata}}}

\examples{



data(phoneme)
mlearn<-phoneme$learn[1:4,1:150]
# Center curves
fdata.c=fdata.cen(mlearn)$Xcen
par(mfrow=c(2,1))
plot.fdata(mlearn,type="l")
plot.fdata(fdata.c,type="l")

# Convert  from class fda to fdata
bsp1 <- create.bspline.basis(c(1,150),21)
fd1 <- Data2fd(1:150,y=t(mlearn$data),basisobj=bsp1)
fdataobj=fdata(fd1)

# Convert  from class fds, fts or sfts to fdata
#require(fds)
#a=fds(x = 1:20, y = Simulationdata$y, xname = "x", yname = "Simulated value")
#b=fts(x = 15:49, y = Australiasmoothfertility$y, xname = "Age",
#    yname = "Fertility rate")
#c=sfts(ts(as.numeric(ElNino$y), frequency = 12), xname = "Month",
#yname = "Sea surface temperature")
#class(a);class(b);class(c)
#fdataobj=fdata(b)

}
\keyword{manip}  