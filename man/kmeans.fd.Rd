\name{kmeans.fd}
\Rdversion{1.1}
\alias{kmeans.fd}
\alias{kmeans.center.ini}
\alias{kmeans.assig.groups}
\alias{kmeans.centers.update}
\title{ K-Means Clustering for functional data}
\description{
	 Perform k-means clustering on functional data. 
}
\usage{
kmeans.fd(fdataobj,ncl=2,metric=metric.lp,dfunc=func.trim.FM,
max.iter=100,par.metric=NULL,par.dfunc=list(trim=0.05),
par.ini=list(method="sample"),draw=TRUE,...)
kmeans.center.ini(fdataobj,ncl=2,metric=metric.lp,
draw=TRUE,method="sample",iter=100,par.metric=NULL,...)
 }
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}     
 	\item{ncl}{ See details section.}
  \item{metric}{ Metric function, by default \code{\link{metric.lp}}.}
	\item{dfunc}{ Type of depth measure, by default \code{depth.mode}}  
   \item{max.iter}{ Maximum number of iterations for the detection of centers. }  	
  \item{draw}{ =TRUE, draw the curves in the color of the centers.}		
  \item{par.dfunc}{ List of arguments to pass to the \code{dfunc} function .}			
  \item{par.ini}{ List of arguments to pass to the \code{kmeans.center.ini} function .}
	\item{method}{ Method for selecting initial centers. If \code{method}=\emph{"Sample"} (by default) takes \code{n} times a random selection by the \code{ncl} centers. The \code{ncl} curves  with  greater distance are the initial centers. If \code{method}=\emph{"Exact"} calculated all combinations of \code{ncl} centers. The \code{ncl} curves  with  greater distance are the initial centers (this method may be too slow).}  	
	\item{iter}{ Maximum number of random samples for the initial selection of centers.}	
  \item{par.metric}{ List of arguments to pass to the \code{metric} function.}
	\item{\dots}{ Further arguments passed to or from other methods.}
}
\details{ 
The method searches the locations around which are grouped data (for a predetermined number of groups).\cr

If \code{ncl=NULL}, randomizes the initial centers, \code{ncl=2} using \code{kmeans.center.ini} function.\cr
If \code{ncl} is an integer,  indicating the number of groups to classify,\cr
are selected \code{ncl} initial centers using \code{kmeans.center.ini} function.\cr
If \code{ncl} is a vector of integers, indicating the position of the initial centers with \code{length(ncl)} equal to number of groups.\cr
If \code{ncl} is a \code{fdata} class objecct, \code{ncl} are the initial centers curves with \code{nrow(ncl)} number of groups.\cr
}
\value{
Return:\cr
\item{cluster}{  Indexes of groups assigned.}
\item{centers}{ Curves centers.}
\item{lcenters}{  Indexes of initial curves centers.}
}
\references{
Hartigan, J. A. and Wong, M. A. (1979). \emph{A K{-}means clustering algorithm}. Applied Statistics 28, 100 \-108.}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{
See Also generic \link[stats]{kmeans} function.
}
\examples{

data(phoneme)
mlearn<-phoneme$learn[c(1:50,101:150,201:250),]

#Unsupervised classification
out.fd1=kmeans.fd(mlearn,ncl=3,draw=TRUE)
# NO RUN
#out.fd2=kmeans.fd(mlearn,ncl=3,draw=TRUE,par.ini=list(method="exact"))
# Different Depth function
#ind=c(17,77,126)
#out.fd3=kmeans.fd(mlearn,ncl=mlearn[ind,],draw=FALSE,
#dfunc=func.trim.FM,par.dfunc=list(trim=0.1))
#out.fd4=kmeans.fd(mlearn,ncl=mlearn[ind,],draw=FALSE,
#dfunc=func.med.FM)
#out.fd5=kmeans.fd(mlearn,ncl=3,dfunc=func.trim.RPD,
#max.iter=10,par.dfunc=list(dfunc="depth.FM",deriv=c(0,1,1)))
#group=c(rep(1,50),rep(2,50),rep(3,50))
#table(out.fd5$cluster,group)
}
\keyword{ cluster }