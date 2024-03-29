% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/depth.fdata.r
\name{depth.fdata}
\alias{depth.fdata}
\alias{depth.mode}
\alias{Depth}
\alias{depth.FM}
\alias{depth.RP}
\alias{depth.RPD}
\alias{depth.RT}
\alias{depth.FSD}
\alias{depth.KFSD}
\title{Computation of depth measures for functional data}
\usage{
depth.mode(
  fdataobj,
  fdataori = fdataobj,
  trim = 0.25,
  metric = metric.lp,
  h = NULL,
  scale = FALSE,
  draw = FALSE,
  ...
)

depth.RP(
  fdataobj,
  fdataori = fdataobj,
  trim = 0.25,
  nproj = 50,
  proj = "vexponential",
  dfunc = "TD1",
  par.dfunc = list(),
  scale = FALSE,
  draw = FALSE,
  ...
)

depth.RPD(
  fdataobj,
  fdataori = fdataobj,
  nproj = 20,
  proj = 1,
  deriv = c(0, 1),
  trim = 0.25,
  dfunc2 = mdepth.LD,
  method = "fmm",
  draw = FALSE,
  ...
)

depth.RT(
  fdataobj,
  fdataori = fdataobj,
  trim = 0.25,
  nproj = 10,
  proj = 1,
  xeps = 1e-07,
  draw = FALSE,
  ...
)

depth.KFSD(
  fdataobj,
  fdataori = fdataobj,
  trim = 0.25,
  h = NULL,
  scale = FALSE,
  draw = FALSE
)

depth.FSD(
  fdataobj,
  fdataori = fdataobj,
  trim = 0.25,
  scale = FALSE,
  draw = FALSE
)

depth.FM(
  fdataobj,
  fdataori = fdataobj,
  trim = 0.25,
  scale = FALSE,
  dfunc = "FM1",
  par.dfunc = list(scale = TRUE),
  draw = FALSE
)
}
\arguments{
\item{fdataobj}{The set of new curves to evaluate the depth.
\code{\link{fdata}} class object.}

\item{fdataori}{The set of reference curves respect to which the depth is
computed.  \code{\link{fdata}} class object.}

\item{trim}{The alpha of the trimming.}

\item{metric}{Metric function, by default \code{\link{metric.lp}}. Distance
matrix between \code{fdataobj} and \code{fdataori}.}

\item{h}{Bandwidth parameter. 
 \itemize{ 
\item If \code{h} is a numerical value, the procedure considers the argument
value as the bandwidth.  
\item If is \code{NULL} (by default) the bandwidth is provided as the 
 15\%--quantile of the distance among curves of
\code{fdataori}.
\item If \code{h} is a character string (like \code{"0.15"}), the procedure
reads the numeric value and consider it as the quantile of the distance in
\code{fdataori} (as in the second case).  
}}

\item{scale}{=TRUE, the depth is scaled respect to depths in
\code{fdataori}.}

\item{draw}{=TRUE, draw the curves, the sample median and trimmed mean.}

\item{\dots}{Further arguments passed to or from other methods. For
\code{depth.mode} parameters for \code{metric}. For random projection
depths, parameters to be included in \code{rproc2fdata} not included before.}

\item{nproj}{The number of projections. Ignored if a \code{fdata} class
object is provided in \code{proj}}

\item{proj}{if a \code{fdata} class, projections provided by the user.
Otherwise, it is the \code{sigma} parameter of \code{\link{rproc2fdata}}
function.}

\item{dfunc}{type of univariate depth function used inside depth function:
"FM1" refers to the original Fraiman and Muniz univariate depth (default),
"TD1" Tukey (Halfspace),"Liu1" for simplical depth, "LD1" for Likelihood
depth and "MhD1" for Mahalanobis 1D depth. Also, any user function
fulfilling the following pattern \code{FUN.USER(x,xx,...)} and returning a
\code{dep} component can be included.f}

\item{par.dfunc}{List of parameters for \emph{dfunc}.}

\item{deriv}{Number of derivatives described in integer vector \code{deriv}.
\code{=0} means no derivative.}

\item{dfunc2}{Multivariate depth function (second step depth function) in
RPD depth, by default \code{\link{mdepth.LD}}. Any user function with the
pattern \code{FUN.USER(x,xx,...)} can be employed.}

\item{method}{Type of derivative method. See \code{\link{fdata.deriv}} for
more details.}

\item{xeps}{Accuracy. The left limit of the empirical distribution function.}
}
\value{
Return a list with:
\itemize{
\item {median}{ Deepest curve.} 
\item {lmed}{ Index deepest element \code{median}.}
\item {mtrim}{ \code{fdata} class object with the average from the \code{(1-trim)\%} deepest curves. }
\item {ltrim}{ Indexes of curves that conform the trimmed mean \code{mtrim}. }
\item {dep}{ Depth of each curve of fdataobj w.r.t. fdataori.}
\item {dep.ori}{ Depth of each curve of fdataori w.r.t. fdataori.}
\item {proj}{ The projection value of each point on the curves. } 
\item {dist}{ Distance matrix between curves or functional data.}
}
}
\description{
Several depth measures can be computed for functional data for descriptive
or classification purposes.
}
\details{
Type of depth functions: Fraiman and Muniz (FM)
depth, modal depth, random Tukey (RT), random projection (RP) depth and
double random projection depth (RPD).
\itemize{ 
\item \code{\link{depth.FM}} computes the integration of an univariate depth
along the axis x (see Fraiman and Muniz 2001). It is also known as
Integrated Depth.

\item \code{\link{depth.mode}} implements the modal depth (see Cuevas et al
2007).

\item \code{\link{depth.RT}} implements the Random Tukey depth (see
Cuesta--Albertos and Nieto--Reyes 2008).

\item \code{\link{depth.RP}} computes the Random Projection depth (see
Cuevas et al. 2007).

\item \code{\link{depth.RPD}} implements a depth measure based on random
projections possibly using several derivatives (see Cuevas et al. 2007).

\item \code{\link{depth.FSD}} computes the Functional Spatial Depth (see
Sguera et al. 2014).

\item \code{\link{depth.KFSD}} implements the Kernelized Functional Spatial
Depth (see Sguera et al. 2014).  
} 
\itemize{ 
\item The \code{\link{depth.mode}} function calculates the depth of a datum
accounting the number of curves in its neighbourhood. By default, the
distance is calculated using \code{\link{metric.lp}} function although any
other distance could be employed through argument \code{metric} (with the
general pattern \code{USER.DIST(fdataobj,fdataori)}).
\item The \code{\link{depth.RP}} function summarizes the random projections
through averages whereas the \code{\link{depth.RT}} function uses the
minimum of all projections.
\item The \code{\link{depth.RPD}} function involves the original
trajectories and the derivatives of each curve in two steps. It builds
random projections for the function and their derivatives (indicated in the
parameter \code{deriv}) and then applies a depth function (by default
\code{\link{depth.mode}}) to this set of random projections (by default the
Tukey one).
\item The \code{\link{depth.FSD}} and \code{\link{depth.KFSD}} are the
implementations of the default versions of the functional spatial depths
proposed in Sguera et al 2014. At this moment, it is not possible to change
the kernel in the second one.#' 
}
}
\examples{
\dontrun{
#Ex: CanadianWeather data
tt=1:365
fdataobj<-fdata(t(CanadianWeather$dailyAv[,,1]),tt)
# Fraiman-Muniz Depth
out.FM=depth.FM(fdataobj,trim=0.1,draw=TRUE)
#Modal Depth
out.mode=depth.mode(fdataobj,trim=0.1,draw=TRUE)
out.RP=depth.RP(fdataobj,trim=0.1,draw=TRUE)
out.RT=depth.RT(fdataobj,trim=0.1,draw=TRUE)
out.FSD=depth.FSD(fdataobj,trim=0.1,draw=TRUE)
out.KFSD=depth.KFSD(fdataobj,trim=0.1,draw=TRUE)
## Double Random Projections
out.RPD=depth.RPD(fdataobj,deriv=c(0,1),dfunc2=mdepth.LD,
trim=0.1,draw=TRUE)
out<-c(out.FM$mtrim,out.mode$mtrim,out.RP$mtrim,out.RPD$mtrim)
plot(fdataobj,col="grey")
lines(out)
cdep<-cbind(out.FM$dep,out.mode$dep,out.RP$dep,out.RT$dep,out.FSD$dep,out.KFSD$dep)
colnames(cdep)<-c("FM","mode","RP","RT","FSD","KFSD")
pairs(cdep)
round(cor(cdep),2)
}

}
\references{
Cuevas, A., Febrero-Bande, M., Fraiman, R. (2007). Robust
estimation and classification for functional data via projection-based depth
notions. \emph{Computational Statistics} 22, 3, 481-496.

Fraiman R, Muniz G. (2001). Trimmed means for functional data. \emph{Test}
10: 419-440.

Cuesta--Albertos, JA, Nieto--Reyes, A. (2008) The Random Tukey Depth.
\emph{Computational Statistics and Data Analysis} Vol. 52, Issue 11,
4979-4988.

Febrero-Bande, M, Oviedo de la Fuente, M. (2012).  Statistical Computing in
Functional Data Analysis: The R Package fda.usc. \emph{Journal of
Statistical Software}, 51(4), 1-28. \url{https://www.jstatsoft.org/v51/i04/}

Sguera C, Galeano P, Lillo R (2014). Spatial depth based classification for
functional data. \emph{TEST} 23(4):725--750.
}
\seealso{
See Also as \code{\link{Descriptive}}.
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente
\email{manuel.oviedo@udc.es}
}
\keyword{descriptive}
