\name{Utilities}
\Rdversion{1.1}
\alias{split.fdata}
\alias{unlist.fdata}
\alias{split.fdata}
\title{ A wrapper for the split and unlist function for fdata object} 
%Divide into Groups the functional data and Flatten Lists
\description{\code{split.fdata} divides the data in the fdata object \code{x} into the groups defined by \code{f}. \cr
Given a list structure x, \code{unlist} simplifies it to produce a \code{fdata} class object which contains all the atomic components which occur in \code{x}.\cr
}
\usage{
\method{split}{fdata}(x,f,drop=FALSE,\dots)
\method{unlist}{fdata}(x, recursive = TRUE, use.names = TRUE)
}
\arguments{
\item{x}{ an \code{fdata} object containing values to be divided into groups or an \code{list} of \code{fdata} objects  containing values to be combine by rows in a to be flatten one \code{fdata} object}
\item{f}{	a factor in the sense that as.factor(f) defines the grouping, or a list of such factors in which case their interaction is used for the grouping.}
\item{drop}{logical indicating if levels that do not occur should be dropped (if f is a factor or a list).}
\item{\dots}{further potential arguments passed to methods.}
\item{recursive}{logical. Should unlisting be applied to list components of x?}
\item{use.names}{ logical. Should names be preserved?}
 }   
\value{
The value returned from \code{split} is a list of fdata objects containing the values for the groups. The components of the list are named by the levels of f (after converting to a factor, or if already a factor and drop = TRUE, dropping unused levels).\cr
The value returned from \code{unlist} is a \code{fdata} object containing the \code{fdata} components of the list.
}
\author{
Manuel Febrero Bande and Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
} 
\examples{
fdataobj<-fdata(MontrealTemp)
fac<-factor(c(rep(1,len=17),rep(2,len=17)))
a1<-split.fdata(fdataobj,fac)
dim(a1[[1]]);dim(a1[[2]])
a2<-unlist.fdata(a1)
a2==fdataobj
} 
\keyword{descriptive} 