\name{Depth for univariate fdata}
\Rdversion{1.1}
\alias{Depth}
\alias{depth.FM}
\alias{depth.mode}
\alias{depth.RP}
\alias{depth.RPD}
\alias{depth.RT}
\alias{depth.FSD}
\alias{depth.KFSD}
\title{ Provides the depth measure for functional data}
\description{
 Compute measure of centrality of the functional data. Type of depth function: 
 Fraiman and Muniz (FM) depth, modal depth, random tukey (RT), random project (RP) depth and double random project depth (RPD).\cr
\itemize{

\item  \code{\link{depth.FM}} computes the integration of an univariate depth along the axis x (see Fraiman and Muniz 2001). It is also known as Integrated Depth. 

\item \code{\link{depth.mode}} implements the modal depth (see Cuevas et al 2007).    

\item \code{\link{depth.RT}} implements the Random Tukey depth (see Cuesta--Albertos and Nieto--Reyes 2008).

\item \code{\link{depth.RP}} computes the Random Projection depth (see Cuevas et al. 2007).   

\item \code{\link{depth.RPD}} implements a depth measure based on random projections possibly using several derivatives (see Cuevas et al. 2007).   

\item \code{\link{depth.FSD}} computes the Functional Spatial Depth (see Sguera et al. 2014).   

\item \code{\link{depth.KFSD}} implements the Kernelized Functional Spatial Depth (see Sguera et al. 2014).
 }
}
\usage{
depth.FM(fdataobj,fdataori=fdataobj,trim=0.25,scale=FALSE,
dfunc="FM1",par.dfunc = list(scale = TRUE), draw = FALSE)

depth.mode(fdataobj,fdataori=fdataobj,trim=0.25,
metric=metric.lp,h=NULL,scale=FALSE,draw=FALSE,...)

depth.RT(fdataobj,fdataori=fdataobj, trim = 0.25, nproj = 10, 
proj = 1, xeps = 1e-07, draw = FALSE, ...)

depth.RP(fdataobj,fdataori=fdataobj,trim=0.25,nproj=50,
proj="vexponential",dfunc="TD1",par.dfunc=list(),scale=FALSE,
draw=FALSE,...)

depth.RPD(fdataobj,fdataori=fdataobj,nproj=50,proj=1,deriv=c(0,1),
trim=0.25,dfunc2=depth.mode,method="fmm",draw=FALSE,...)  

depth.FSD(fdataobj, fdataori = fdataobj, 
          trim = 0.25, scale = FALSE, draw = FALSE)
          
depth.KFSD(fdataobj, fdataori = fdataobj, trim = 0.25,
           h=NULL,scale = FALSE, draw = FALSE)
}

\arguments{
 \item{fdataobj}{  A set of new curves to evaluate the depth. \code{\link{fdata}} class object. }    
 \item{fdataori}{  A set of original curves where the depth is computed.  \code{\link{fdata}} class object.}  
  \item{trim}{ The alpha of the trimming.}
  \item{nproj}{ The number of projections.}
    %\item{proj}{ Random projections.  If is a character: create the random projection using a covariance matrix by process indicated in the argument (by default, proj=1, sigma=diag(ncol(fdataobj))), else if is a  matrix or fdata, it is a random projection provided by the user.} 
    
\item{proj}{   matrix or fdata class object with a random projection provided by the user. If is a character: create the random projection using a covariance matrix by process indicated in the argument. Otherwise, it is a \code{sigma} parameter of \code{\link{rproc2fdata}} function.}         
        
\item{dfunc}{type of univariate depth function used inside depth function: "FM1" refers to the original Fraiman and Muniz univariate depth (default), "TD1" Tukey (Halfspace),"Liu1" for simplical depth, "LD1" for Likelihood depth and "MhD1" for Mahalanobis 1D depth. Also, any user function fulfilling the following pattern \emph{FUN.USER(x,xx,...)} and returning a \code{dep} component can be included.}
  \item{par.dfunc}{List of parameters for \emph{dfunc}.}
  
  \item{dfunc2}{Second depth function in  RPD depth, by default \code{\link{depth.mode}}.}
  
  \item{deriv}{  Number of derivatives described in integer vector \code{deriv}. \code{=0} means no derivative.}
  
  \item{method}{  Type of derivative method. See   \code{\link{fdata.deriv}} for more details.}
 
 \item{h}{ Bandwidth parameter. 
  \enumerate{
  
  \item If \code{h} is a numerical value, the procedure considers the argument value as bandwidth.
 \item If is \code{NULL} (by default) the bandwidth is provided as the 15\%--quantile of the distance between \code{fdataobj} and \code{fdataori}.
  
  \item If \code{h} is a character string (like \code{"q-0.15"}), the procedure reads the numeric value from the third position of the character to the end and uses it to compute the quantile of the distance between \code{fdataobj} and \code{fdataori} (as in the second case).
  }}
 
 \item{metric}{ Metric function, by default \code{\link{metric.lp}}. Distance matrix between \code{fdataobj} and \code{fdataori}.}
  
  \item{scale}{ =TRUE, scale the depth.}   

  \item{xeps}{ Accuracy. The left limit  of the empirical distribution function.}
  
  \item{draw}{ =TRUE, draw the curves, the sample median and trimmed mean.}
  
  \item{\dots}{ Further arguments passed to or from other methods.}
}

\details{
 \itemize{
 
\item The modal depth \code{\link{depth.mode}} function calculates the depth of a datum accounting the number of curves in its neighbourhood. By default, the distance is calculated using \code{\link{metric.lp}} function although any other distance could be employed through argument \code{metric} (with the general pattern \code{USER.DIST(fdataobj,fdataori)}).

\item The \code{\link{depth.RP}} function summarizes the random projections through averages whereas the \code{\link{depth.RT}} function uses the minimum of all projections. 
 
\item The  \code{\link{depth.RPD}} function involves the original trajectories and the derivatives of each curve in two steps. It builds random projections for the function and their derivatives (indicated in the parameter \code{deriv}) and then applies a depth function  (by default \code{\link{depth.mode}}) to this set of random projections (by default the Tukey one).

\item The \code{\link{depth.FSD}} and \code{\link{depth.KFSD}} are the implementations of the default versions of the functional spatial depths proposed in Sguera et al 2014. At this moment, it is not possible to change the kernel in the second one. 

}
}

\value{
    \item{median}{ Deepest curve.}
    \item{lmed}{ Index deepest element \code{median}.}
    \item{mtrim}{ \code{fdata} class object with the average from the \code{(1-trim)\%} deepest curves. }
    \item{ltrim}{ Indexes of curves that conform the trimmed mean \code{mtrim}. }    
    \item{dep}{ Depth of each curve of fdataobj w.r.t. fdataori.}
    \item{dep.ori}{ Depth of each curve of fdataori w.r.t. fdataori.}        
    \item{proj}{ The projection value of each point on the curves. }
    \item{dist}{ Distance matrix between curves or functional data.}    
}

\references{
Cuevas, A., Febrero-Bande, M., Fraiman, R. (2007). Robust estimation and classification for functional data via projection-based depth notions. \emph{Computational Statistics} 22, 3, 481{-}496.

%Febrero-Bande, M., Galeano, P., and Gonzalez-Manteiga, W. (2008).  \emph{Outlier detection in functional data by depth measures with application to identify abnormal NOx levels}. Environmetrics 19, 4, 331{-}345.

Fraiman R, Muniz G. (2001). Trimmed means for functional data. \emph{Test} 10: 419{-}440.

Cuesta--Albertos, JA, Nieto--Reyes, A. (2008)  The Random Tukey Depth.
\emph{Computational Statistics and Data Analysis} Vol. 52, Issue 11, 4979{-}4988.

Febrero-Bande,  M, Oviedo de la Fuente, M. (2012).  Statistical Computing in Functional Data Analysis: The R Package fda.usc.
\emph{Journal of Statistical Software}, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}

Sguera C, Galeano P, Lillo R (2014). Spatial depth based classification for functional data. \emph{TEST} 23(4):725--750.
}

\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{See Also as \code{\link{Descriptive}}.}
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
out.RPD=depth.RPD(fdataobj,deriv=c(0,1),dfunc2=depth.FM,
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
\keyword{descriptive} 