% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ops.fda.usc.r
\name{ops.fda.usc}
\alias{ops.fda.usc}
\title{ops.fda.usc  Options Settings}
\usage{
ops.fda.usc(
  verbose = FALSE,
  trace = FALSE,
  warning = FALSE,
  ncores = NULL,
  int.method = "TRAPZ",
  reset = FALSE,
  eps = as.double(.Machine[[1]] * 10)
)
}
\arguments{
\item{verbose}{\code{logical}. Should R report extra information on progress? Set to \code{TRUE} by the command-line option --verbose.}

\item{trace}{\code{logical}. Show internal information of procedure.}

\item{warning}{\code{logical}: If true, warnings are shown.}

\item{ncores}{integer. Number of CPU cores on the current host.}

\item{int.method}{see \code{method} argument in \code{\link{int.simpson}} function.}

\item{reset}{\code{logical}.  If \code{TRUE} creates a new Parallel Socket Cluster (ncores>1) or a sequential parallel backend (ncores=1). It is useful when worker initialization failed or after a crush.}

\item{eps}{epsilon parameter.}
}
\description{
Set or query graphical and prompt output parameters. Allow the user to set and examine a variety of global or local options which affect the way in which fda.usc functions computes and displays its results.
}
\examples{
\dontrun{
# If worker initialization failed, please execute this code
 ncores <- max(parallel::detectCores() -1,1)
 if (ncores==1) {
     foreach::registerDoSEQ()
 }  else{
 cl <-  suppressWarnings(parallel::makePSOCKcluster(ncores ))
 doParallel::registerDoParallel(cl)
 }
 ops.fda.usc()
}

}
\author{
Manuel Oviedo de la Fuente (\email{manuel.oviedo@udc.es}).
}
