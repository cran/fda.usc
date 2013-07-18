\name{plot.fdata}
\Rdversion{1.1}
\alias{plot.fdata}
\alias{lines.fdata}
\alias{title.fdata}
\title{ Plot functional data: fdata.
}
\description{
Plot object of class \code{fdata}.
}
\usage{
\method{plot}{fdata}(x,type,main,xlab,ylab,\dots)
\method{lines}{fdata}(x,\dots)
\method{title}{fdata}(x,main=NULL,xlab=NULL,ylab=NULL,rownames=NULL) 
}
\arguments{
  \item{x}{ \code{fdata} class object with: 
\itemize{
		\item \code{"data"}: matrix  of  set cases with dimension (\code{n} x \code{m}), where \code{n} is the number of curves and \code{m} are the points observed in each curve
    \item \code{"rangeval"}: the discretizations points values, if not provided:  \code{1:m}
    \item \code{"rangeval"}: range of the discretizations points values, range(\code{argvals}) by default. 
    \item \code{"names"}: (optional) list with  \code{main} an overall title,	\code{xlab} title for \code{x} axis and \code{ylab} title for \code{y} axis.
    }  }
%   i.- \code{"data"}: matrix  of  set cases with dimension (\code{n} x \code{m}), where \code{n} is the number of curves and \code{m} are the discretization points observed in each curve, ii.- \code{"argvals"}: the discretizations points values, iii.- \code{"rangeval"}: range of the discretizations points values, iv.- \code{"names"} (optional) list with  \code{main} an overall title, \code{xlab} title for \code{x} axis and \code{ylab} title for \code{y} axis.}   
  \item{type}{1-character string giving the type of plot desired. The following values are possible, for details, see plot: "l" for lines (by default),"p" for points, , "o" for overplotted points and lines, "b", "c" for (empty if "c") points joined by lines, "s" and "S" for stair steps and "h" for histogram-like vertical lines. Finally, "n" does not produce any points or lines.}
  \item{main}{an overall title for the plot: see \code{\link{title}}.}
  \item{xlab}{xlab title for x axis, as in plot.}  
  \item{ylab}{ylab title for y axis, as in plot.}	
  \item{rownames}{row names.}		
  \item{\dots}{ Further arguments passed to \link[graphics]{matplot} function.}
  }
\details{  
Plot \code{fdata} object.
}
\author{
Manuel Febrero Bande and Manuel Oviedo de la Fuente <manuel.oviedo@usc.es>
}
\seealso{
See Also as \code{\link{fdata}} 
}

\examples{
a1<-seq(0,1,by=.01)
a2=rnorm(length(a1),sd=0.2)
f1<-(sin(2*pi*a1))+rnorm(length(a1),sd=0.2)
nc<-10
np<-length(f1)
tt=seq(0,1,len=101)
mdata<-matrix(NA,ncol=np,nrow=nc)
for (i in 1:nc) mdata[i,]<- (sin(2*pi*a1))+rnorm(length(a1),sd=0.2)
fdataobj<-fdata(mdata,tt)
res=plot.fdata(fdataobj,type="l",col=gray(1:nrow(mdata)/nrow(mdata)))
lines(func.mean(fdataobj),col=3,lwd=2) #original curve
}

\keyword{hplot} 