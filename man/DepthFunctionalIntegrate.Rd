\name{Depth for multivariate fdata}
\Rdversion{1.1}
\alias{Depth.pfdata}
\alias{depth.FMp}
\alias{depth.modep}
\alias{depth.RPp}
\title{ Provides the depth measure for a list of p--functional data objects}
\description{
 This function computes the depth measure for a list of p--functional data objects.
 The procedure extends the  Fraiman and Muniz (FM), modal, and random project depth functions from 1 functional dataset to p functional datasets.\cr
}
\usage{
depth.FMp(lfdata, lfdataref = lfdata, trim = 0.25, dfunc = "mdepth.MhD",
par.dfunc = list(scale = FALSE), draw = FALSE, ask = FALSE, ...) 

depth.RPp(lfdata, lfdataref = lfdata, nproj = 50, proj = 1, trim = 0.25, 
dfunc = "mdepth.MhD", par.dfunc = list(scale = TRUE), draw = FALSE,ask = FALSE)  
   
depth.modep(lfdata, lfdataref = lfdata, h = NULL, metric, par.metric = list(), 
method = "euclidean", scale = FALSE, trim = 0.25, draw = FALSE,ask = FALSE)
}

\arguments{
 \item{lfdata}{  A list of new curves (list of fdata ojects) to evaluate the depth.  }    
 \item{lfdataref}{  A set of reference curves (list of fdata ojects)  w.r.t. the depth of {lfdata} is computed.}  
  \item{trim}{ The alpha of the trimming.}
  \item{dfunc}{Type of multivariate depth (of order p) function used in Framiman and Muniz depth, \code{depth.FMp} or in Random Projection depth,\code{depth.FMp}: : 
  \itemize{
 \item  The \code{\link{mdepth.SD}} function provides the simplicial depth measure for bivariate data.  
  \item  The \code{\link{mdepth.LD}} function provides the Likelihood depth measure based on Nadaraya--Watson estimator of empirical density function.  
   \item The \code{\link{mdepth.HS}} function implements a half-space depth measure based on random projections.
   \item The \code{\link{mdepth.MhD}}function implements a Mahalanobis depth measure.
   \item The \code{\link{mdepth.RP}} function provides the depth measure using random projections for multivariate data.
      }
  }
  \item{par.dfunc}{list of parameters for the \code{dfunc} depth function, see \code{\link{Depth.Multivariate}}.}  
  \item{nproj}{ The number projection.}
  \item{proj}{  if is a character: create the random projection using a covariance matrix by process indicated in the argument (by default, proj=1, sigma=diag(ncol(fdataobj))), else if is a matrix of random projection provided by the user.} 
  \item{h}{ Bandwidth, \code{h>0}. Default argument values are provided as the 15\%--quantile of the distance between \code{fdataobj} and \code{fdataori}.}
  \item{metric}{ Metric or semi--metric function used for compute the distance between each element in \code{ldata} w.r.t. \code{ldataref}, by default \code{\link{metric.lp}}.}
  \item{par.metric}{list of parameters for the metric function.}  
  \item{method}{ Type of the distance measure (by default \code{euclidean}) to compute the metric between the p-distance matrix computed from the p functional data elements.}  
  \item{scale}{ =TRUE, scale the depth.}     
  \item{draw}{ =TRUE, draw the curves, the sample median and trimmed mean.}
  \item{ask}{logical. If \code{TRUE} (and the R session is interactive) the user is asked for input, before a new figure is drawn.} 
  \item{\dots}{ Further arguments passed to or from other methods.}
}

\details{

 \itemize{
  \item \code{\link{depth.FMp}}, this procedure suposes that each curve of the lfdataobj have the same support [0,T] (same argvals and rangeval).  The FMp depth is defined as:
\eqn{FM_i^p =\int_{0}^{T}Z_i^p(t)dt}{FM_i^p =\int Z_i^p(t)dt} where \eqn{Z_i^p(t)}{Z_i^p(t)} is a \eqn{p}--variate depth of the vector
 \eqn{(x_i^1(t),\ldots,x_i^p(t))}{(x_i^1(t),...,x_i^p(t))} w.r.t. the sample at \eqn{t}.
% This happens, for instance, when using the curves and its derivatives. In this case,note solo un dato funcional se reduce depth.FM=depth.FM1

 
 \item The  \code{\link{depth.RPp}} function calculates the depth in two steps. It builds random projections for the  each curve of the \code{lfdata} w.r.t. each curve of the \code{lfdataref} object. Then it applyes a multivariate depth function specified in \code{dfunc} argument to the set of random projections.
 This procedure is a generalization of Random Projection with derivatives (RPD)
 implemented in \code{\link{depth.RPD}} function.  Now, the procedure computes a p-variate depth with the projections using the \eqn{p} functional dataset.  
  
  \item The modal depth \code{\link{depth.modep}} function calculates the depth in three steps. 
  First,  the function calculates a suitable metrics or semi--metrics   \eqn{m_1+\cdots+m_p}{m1,...,mp} 
  for each curve of the  \code{lfdata} w.r.t. each curve in the \code{lfdataref} object using
   the \code{metric} and \code{par.metric} arguments,  see \code{\link{metric.lp}} or \code{\link{semimetric.NPFDA}}
    for more details.  Second, the function  uses the \eqn{p}--dimensional metrics to construct a new metric, 
    specified in \code{method} argument, by default if \code{method="euclidean"}, i.e.
     \eqn{m:=\sqrt{m_1^2+\cdots+m_p^2}}{m:=\sqrt(m1^2+...+mp^2)}. Finally, the empirical \emph{h}--depth is computed as:  
     \deqn{\hat{f}_h(x_0)=N^{-1}\sum_{i=1}^{N}{K(m/h)}}{\hat{f}(x)=1/N\sum K(m/h)}
where \eqn{x} is dataset with p observed fucntional data,  \eqn{m} is a suitable metric or semi--metric, 
\eqn{K(t)} is an asymmetric kernel function and \code{h} is the bandwidth  parameter.
 }
}
 
\value{
  \item{lmed}{ Index deepest element \code{median}.}
  \item{ltrim}{ Index of curves with trimmed mean \code{mtrim}.}   
  \item{dep}{ Depth of each curve of fdataobj w.r.t. fdataori.}     
  \item{dfunc}{second depth function used as multivariate depth, see details section.}
  \item{par.dfunc}{list of parameters for the dfunc depth function.}  
  \item{proj}{ The projection value of each point on the curves. }
  \item{dist}{ Distance matrix between curves or functional data.}    
}

\references{
Cuevas, A., Febrero-Bande, M. and Fraiman, R. (2007). \emph{Robust estimation and classification for functional data via projection-based depth notions.} Computational Statistics 22, 3, 481{-}496.
%Fraiman R, Muniz G. 2001.  \emph{Trimmed means for functional data.} Test 10: 419{-}440.
%Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
}

\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}

\seealso{See Also as \code{\link{Descriptive}}.}
\examples{
\dontrun{
data(tecator)
xx<-tecator$absorp
xx1<-fdata.deriv(xx,1)
xx1<-fdata.deriv(xx,2)
proj<-rproc2fdata(111,c(1,2))$data

 # Fraiman-Muniz Depth
out.FM=depth.FM(xx,trim=0.1,dfunc="FM1",draw=TRUE)
out.FM1p=depth.FMp(list(xx),trim=0.1,dfunc="FM1",draw=TRUE)
out.FM$dep-out.FM1p$dep
out.FM2p=depth.FMp(list(xx,xx1),trim=0.1,dfunc="SD",
draw=TRUE,par.dfunc=list(random=FALSE))


#Modal Depth
out.mode=depth.mode(xx,trim=0.1,draw=TRUEmetric=metric.lp)
out.mode1p=depth.modep(list(xx),trim=0.1,metric=metric.lp)
out.mode$dep-out.mode1p$dep
out.mode2p=depth.modep(list(xx),trim=0.1,metric=metric.lp,
method="manhattan'')

 # Fraiman-Muniz Depth
out.RP=depth.RP(xx,trim=0.1,dfunc="FM1",draw=TRUE)
out.RP1p=depth.RPp(list(xx),trim=0.1,dfunc="FM1",draw=TRUE)
out.RP$dep-out.RP1p$dep
## Double Random Projections
out.RPD=depth.RPD(xx,deriv=c(0,1),dfunc2=depth.FM,trim=0.1,draw=TRUE)
out.RP2p=depth.RPp(list(xx,xx1),trim=0.1,dfunc="FM",draw=TRUE,
par.dfunc=list(random=FALSE),proj=out.RPD$proj)
out.RPD$dep-out.RP2p$dep
}
}
\keyword{descriptive} 