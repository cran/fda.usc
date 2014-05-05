\name{na.omit.fdata}
\Rdversion{1.1}
\alias{na.omit.fdata}
\alias{na.fail.fdata}
\title{ A wrapper for the na.omit and na.fail function for fdata object} 
\description{
\code{na.fail} returns the object if it does not contain any missing values, and signals an error otherwise. \code{na.omit} returns the object with incomplete cases removed.\cr
If \code{na.omit.fdata} removes cases, the row numbers of the cases form the \code{"na.action"} attribute of the result, of class \code{"omit"}, see generic function \code{\link{na.omit}}.
}
\usage{
\method{na.omit}{fdata}(object,\dots)
\method{na.fail}{fdata}(object,\dots)
}
\arguments{
\item{object}{ an \code{fdata} object.}
\item{\dots}{further potential arguments passed to methods.}
}   
\value{
The value returned from \code{omit} is a \code{fdata} object with incomplete cases removed.
}
\author{
Manuel Febrero Bande
} 
\examples{
fdataobj<-fdata(MontrealTemp)
fdataobj$data[3,3]<-NA
fdataobj$data[10,]<-NA
fdastaobj2<-na.omit.fdata(fdataobj)
} 
\keyword{descriptive} 