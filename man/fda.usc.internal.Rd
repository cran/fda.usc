\name{fda.usc.internal}
\Rdversion{1.1}
\alias{[.fdata}
\alias{+.fdata}
\alias{-.fdata}
\alias{*.fdata}
\alias{/.fdata}
\alias{c.fdata}
\alias{==.fdata}
\alias{dim.fdata}
\alias{ncol.fdata}
\alias{nrow.fdata}
\alias{traza}
\alias{is.fdata}
\alias{count.na}

\title{ fda.usc internal functions}
\description{Internal undocumentation functions for fda.usc package}
\usage{
\method{[}{fdata}(fdataobj, i = TRUE, j = TRUE,drop=FALSE)
\method{==}{fdata}(fdata1,fdata2)
\method{+}{fdata}(fdata1,fdata2)
\method{-}{fdata}(fdata1,fdata2)
\method{*}{fdata}(fdata1,fdata2)
\method{/}{fdata}(fdata1,fdata2)
\method{c}{fdata}(\dots) 
\method{dim}{fdata}(x)
\method{nrow}{fdata}(fdataobj)
\method{ncol}{fdata}(fdataobj)
\method{is}{fdata}(fdataobj)
count.na(A)
traza(A) 
}
\arguments{
 \item{x,fdataobj,fdata1,fdata2}{ \code{fdata} class object}. %In \code{"[.fdata"} \code{fdataobj} is a object  from which to extract element(s) or in which to replace element(s).}  
  \item{i,j}{ Indices specifying elements to extract or replace. Indices are numeric or character vectors or empty}  
 \item{drop}{For \code{fdata} class object. If TRUE the result is coerced to the lowest possible dimension of  element \code{data}. This only works for extracting elements, not for the replacement.}
   \item{A}{  Matrix  of  set cases with dimension (\code{n} x \code{m}), where \code{n} is the number of curves and \code{m} are the points observed in each curve. \code{traza}, return the trace of a functional data matrix.}
  \item{\dots}{ \code{fdata} objects to be concatenated. }
}
\keyword{math}  