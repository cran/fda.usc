\name{Depth for univariate fdata}
\Rdversion{1.1}
\alias{Depth}
\alias{depth.FM}
\alias{depth.mode}
\alias{depth.RP}
\alias{depth.RPD}
\alias{depth.RT}
\title{ Provides the depth measure for functional data}
\description{
 Compute measure of centrality of the functional data. Type of depth function: 
 Fraiman and Muniz (FM) depth, modal depth, random tukey (RT), random project (RP) depth and double random project depth (RPD).\cr
\itemize{
  \item  The \code{\link{depth.FM}} function is based on integrated an univariate depth measure along the axis x. 
  \item The \code{\link{depth.mode}} function implements the modal depth that selects the curve most densely surrounded by others as the deepest one. By default, the distance is calculated  using \code{\link{metric.lp}} function. 
  \item The \code{\link{depth.RT}} function implements a depth measure based on random projections using a half-space Tukey method.
  \item The \code{\link{depth.RP}} function implements a depth measure based on random projections.  
  \item The \code{\link{depth.RPD}} function implements a depth measure based on random projections and using several derivatives. 
 }
}
\usage{
depth.FM(fdataobj,fdataori=fdataobj,trim=0.25,dfunc="FM1",scale=FALSE,
xeps=1e-8,draw=FALSE,...)

depth.mode(fdataobj,fdataori=fdataobj,trim=0.25,
metric=metric.lp,h=NULL,scale=FALSE,draw=FALSE,...)

depth.RT(fdataobj,fdataori=fdataobj, trim = 0.25, nproj = 10, 
proj = 1, xeps = 1e-07, draw = FALSE, ...)

depth.RP(fdataobj,fdataori=fdataobj,trim=0.25,nproj=50,proj=1,scale=FALSE,
xeps=1e-8,draw=FALSE,...)

depth.RPD(fdataobj,fdataori=fdataobj,nproj=50,proj=1,deriv=c(0,1),trim=0.25,
dfunc2=depth.mode,method="fmm",draw=FALSE,...)  
}

\arguments{
 \item{fdataobj}{  A set of new curves to evaluate the depth. \code{\link{fdata}} class object. }    
 \item{fdataori}{  A set of original curves where the depth is computed.  \code{\link{fdata}} class object.}  
  \item{trim}{ The alpha of the trimming.}
  \item{nproj}{ The number of projections.}
    %\item{proj}{ Random projections.  If is a character: create the random projection using a covariance matrix by process indicated in the argument (by default, proj=1, sigma=diag(ncol(fdataobj))), else if is a  matrix or fdata, it is a random projection provided by the user.} 
    
        \item{proj}{   matrix or fdata class object with a random projection provided by the user. Otherwise, it is a \code{sigma} parameter of \code{\link{rproc2fdata}} function.} 
        
        
        
  \item{dfunc}{type of univariate depth function used in Framiman and Muniz depth, \code{depth.FM}. "TD1" for Tukey (Halfspace) depth, "SD" for simplical (Liu) depth and "FM1" by default for Fraiman and Muniz depth.}
  \item{dfunc2}{second depth function in  RPD depth, by default \code{\link{depth.mode}}.}
  \item{deriv}{  Number of derivatives described in integer vector \code{deriv}. \code{=0} means no derivative.}
  \item{method}{  Type of derivative method. See   \code{\link{fdata.deriv}} for more details.}
  \item{h}{ Bandwidth, \code{h>0}. Default argument values are provided as the 15\%--quantile of the distance between \code{fdataobj} and \code{fdataori}.}
  \item{metric}{ Metric function, by default \code{\link{metric.lp}}. Distance matrix between \code{fdataobj} and \code{fdataori}.}
  \item{scale}{ =TRUE, scale the depth.}   
  \item{xeps}{ Accuracy. The left limit  of the empirical distribution function.}
  \item{draw}{ =TRUE, draw the curves, the sample median and trimmed mean.}
  \item{\dots}{ Further arguments passed to or from other methods.}
}

\details{
 \itemize{
  \item FM depth measures the conditional quantile on all points and then integrated. The \code{\link{depth.FM}} method sorts the depths along the axis x, ie, performs a ranking of depths.  
  \item The modal depth \code{\link{depth.mode}} function calculates the depth of a datum accounting the number of curves in the neighbourhood. The mode of a functional data will be the deepest curve.
  \item The \code{\link{depth.RP}} function calculates depth using univariate depth tools through random projection method (RP). 
  \item The  \code{\link{depth.RPD}} function also involves the derivatives of each curve. This function calculates the depth in two steps. It builds random projections for the funcion and their derivatives (indicated in the parameter \code{deriv}) without looking at the functional data. Then it applyes a depth function  (by default \code{\link{depth.mode}}) to the set of random projections  created for calculate the functional depth by the Tukey method.
  \item The \code{\link{depth.RT}} function builds random projections  and calculates the functional depth by the Tukey method combining the information of all projections.
 }
}

\value{
    \item{median}{ Deepest curve.}
    \item{lmed}{ Index deepest element \code{median}.}
    \item{mtrim}{ \code{fdata} class object with the average from the \code{(1-trim)\%} deepest curves. }
    \item{ltrim}{ Index of curves with trimmed mean \code{mtrim}. }    
    \item{dep}{ Depth of each curve of fdataobj w.r.t. fdataori.}
    \item{dep.ori}{ Depth of each curve of fdataori w.r.t. fdataori.}        
    \item{proj}{ The projection value of each point on the curves. }
    \item{dist}{ Distance matrix between curves or functional data.}    
}

\references{
Cuevas, A., Febrero-Bande, M. and Fraiman, R. (2007). \emph{Robust estimation and classification for functional data via projection-based depth notions.} Computational Statistics 22, 3, 481{-}496.

%Febrero-Bande, M., Galeano, P., and Gonzalez-Manteiga, W. (2008).  \emph{Outlier detection in functional data by depth measures with application to identify abnormal NOx levels}. Environmetrics 19, 4, 331{-}345.

Fraiman R, Muniz G. 2001.  \emph{Trimmed means for functional data.} Test 10: 419{-}440.

Cuesta Albertos, J. A. and Nieto Reyes, A.  \emph{The Random Tukey Depth.}
Computational Statistics and Data Analysis (2008), Vol. 52, Issue 11, 4979{-}4988.

Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
}


\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{See Also as \code{\link{Descriptive}}.}
\examples{
#Ex: CanadianWeather data
tt=1:365
fdataobj<-fdata(t(CanadianWeather$dailyAv[,,1]),tt)
# Fraiman-Muniz Depth
out.FM=depth.FM(fdataobj,trim=0.1,draw=TRUE)
#Modal Depth
out.mode=depth.mode(fdataobj,trim=0.1,draw=TRUE)
out.RP=depth.RP(fdataobj,trim=0.1,draw=TRUE)
out.RT=depth.RT(fdataobj,trim=0.1,draw=TRUE)
\dontrun{
## Double Random Projections
out.RPD=depth.RPD(fdataobj,deriv=c(0,1),dfunc2=depth.FM,trim=0.1,draw=TRUE)
out<-c(out.FM$mtrim,out.mode$mtrim,out.RP$mtrim,out.RPD$mtrim)
 plot(out)
}
}
\keyword{descriptive} 