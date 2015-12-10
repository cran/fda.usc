\name{subset.fdata}
\Rdversion{1.1}
\alias{subset.fdata}
\title{ Subsetting   }
\description{
Return subsets of \code{fdata} which meet conditions.
}
\usage{
\method{subset}{fdata}(x, subset, select, drop = TRUE,\ldots)
%\method{subset}{lfdata}(x, subset, select, drop = TRUE,\ldots)
}
\arguments{  
 \item{x}{ object to be subsetted (\code{fdata} class).}  
 \item{subset}{  logical expression indicating elements or rows to keep.}  
 \item{select}{  logical expression indicating points or columns to keep.}  
 \item{drop}{ passed on to \code{[} indexing operator.}   
 \item{\ldots}{further arguments to be passed to or from other methods.}
}
\value{
 An object similar to \code{x} contain just the selected elements.
}
\seealso{
See \code{\link{subset}} and \code{\link{fdata}}.
}
%\keyword{multivariate} 
