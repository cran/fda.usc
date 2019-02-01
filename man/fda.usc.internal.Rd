\name{fda.usc.internal}
\Rdversion{1.1}
\alias{[.fdata}
\alias{[.fdist}
\alias{+.fdata}
\alias{-.fdata}
\alias{*.fdata}
\alias{/.fdata}
\alias{c.fdata}
\alias{==.fdata}
\alias{!=.fdata}
\alias{dim.fdata}
\alias{ncol.fdata}
\alias{nrow.fdata}
\alias{NROW.fdata}
\alias{NCOL.fdata}
\alias{length.fdata}
\alias{is.fdata}
\alias{*.fdata}
\alias{^.fdata}
\alias{argvals}
\alias{rangeval}
\alias{missing.fdata}
\alias{omit.fdata}
\alias{omit2.fdata}
\alias{is.na.fdata}
\alias{anyNA.fdata}
\alias{count.na.fdata}
\alias{colnames.fdata}
\alias{rownames.fdata}
\alias{argvals}
\alias{rangeval}

\title{ fda.usc internal functions}
\description{Internal undocumentation functions for fda.usc package.\cr
}
\usage{
\method{[}{fdata}(fdataobj, i = TRUE, j = TRUE,drop=FALSE)
\method{[}{fdist}(fdataobj, i = TRUE, j = TRUE,drop=FALSE)
\method{==}{fdata}(fdata1,fdata2)
\method{!=}{fdata}(fdata1,fdata2)
\method{+}{fdata}(fdata1,fdata2)
\method{-}{fdata}(fdata1,fdata2)
\method{*}{fdata}(fdata1,fdata2)
\method{/}{fdata}(fdata1,fdata2)
\method{^}{fdata}(fdataobj,pot)
\method{c}{fdata}(\dots) 
\method{dim}{fdata}(x)
\method{nrow}{fdata}(x)
\method{ncol}{fdata}(x)
\method{NROW}{fdata}(x)
\method{NCOL}{fdata}(x)
\method{colnames}{fdata}(x)
\method{rownames}{fdata}(x)
\method{length}{fdata}(x)
\method{is}{fdata}(fdataobj)
\method{omit}{fdata}(fdataobj,y=NULL)
\method{omit2}{fdata}(fdataobj,index.na=FALSE)
\method{missing}{fdata}(fdataobj,basis=NULL)
\method{is.na}{fdata}(x)
\method{anyNA}{fdata}(x, recursive = FALSE)
count.na.fdata( x )
argvals(fdataobj)
rangeval(fdataobj)
%"[.fdist"(fdataobj, i = TRUE, j = TRUE,drop=FALSE)
}
\arguments{
 \item{x,fdataobj,fdata1,fdata2}{ \code{fdata} class object}. %In \code{"[.fdata"} \code{fdataobj} is a object  from which to extract element(s) or in which to replace element(s).}  
 \item{i,j}{ Indices specifying elements to extract, replace. Indices are numeric or character vectors or empty}  
  \item{pot}{ Numeric value for exponentiation.}  
 \item{drop}{For \code{fdata} class object. If TRUE the result is coerced to the lowest possible dimension of  element \code{data}. This only works for extracting elements, not for the replacement.}
  \item{\dots}{ \code{fdata} objects to be concatenated. }
  \item{y}{ Vector}
  \item{basis}{ fd basis}
\item{index.na}{ Return the index of NA elements}  
\item{recursive}{ should \code{anyNA} be applied recursively to lists and pairlists?}

}
\note{ In "Ops" functions \code{"+.fdata"}, \code{"-.fdata"}, \code{"*.fdata"} and \code{"/.fdata"}: The lengths of the objects \code{fdata1} and \code{fdata2} may be different because operates recycled into minimum size as necessary.
}
\references{
Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
}
\keyword{math}  