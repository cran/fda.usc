\name{fdata.methods}
\Rdversion{1.1}
\alias{Math.fdata}
\alias{Ops.fdata}
\alias{Summary.fdata}

\title{ fdata  S3 Group Generic Functions}
\description{fdata Group generic methods defined for four specified groups of functions, Math, Ops, Summary and Complex.}
\usage{
\method{Math}{fdata}(x, \dots)
\method{Ops}{fdata}(e1, e2)
\method{Summary}{fdata}(\dots, na.rm = FALSE)
}
\arguments{
 \item{x,e1,e2}{ \code{fdata} class object}. 
 \item{\dots}{ Turther arguments passed to methods. }  
 \item{na.rm}{ logical: should missing values be removed? }  
}
\seealso{
See  \link[base]{Ops}, \link[base]{Math}, \link[base]{Summary} and \link[base]{Complex}.
}
\examples{
#######################
#######################
data(tecator)
absor<-tecator$absorp.fdata
absor2<-fdata.deriv(absor,1)
absor<-absor2[1:5,1:4]
absor2<-absor2[1:5,1:4]
sum(absor)
round(absor,4)
log1<-log(absor)
}
\keyword{math}  