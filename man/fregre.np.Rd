\name{fregre.np}
\Rdversion{1.1}
\alias{fregre.np}
\title{ Functional regression with scalar response using non-parametric kernel estimation}
\description{
 	Computes functional regression between functional explanatory variables and scalar response using kernel estimation.
}
\usage{
fregre.np(fdataobj,y,h=NULL,Ker=AKer.norm,metric=metric.lp,
type.S=S.NW,...)
}
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
	\item{y}{ Scalar response with length \code{n}.}
	\item{h}{ Smoothing parameter or bandwidth.}	
	\item{Ker}{ Type of asymmetric kernel used, by default asymmetric normal kernel.}	
  \item{metric}{ Metric function, by default \code{\link{metric.lp}}.}
  \item{type.S}{ Type of smothing matrix \code{S}. By default \code{S} is calculated by Nadaraya-Watson kernel estimator (\code{S.NW}).}   
  \item{\dots}{ Arguments to be passed  for \code{\link{metric.lp}} o other metric function.}
}
\details{ 
The distance between curves is calculated using the \code{\link{metric.lp}} although any other semimetric could be used (see \code{\link{semimetric.basis}} or \code{\link{semimetric.NPFDA}} functions). 
Different asymmetric kernels can be used, see   \code{\link{Kernel.asymmetric}}.\cr
}
\value{
Return:\cr
\item{call}{ The matched call.}
\item{fitted.values}{  Estimated scalar response.}
\item{H}{  Hat matrix.}
\item{residuals}{ \code{y} minus \code{fitted values}.}
\item{df}{ The residual degrees of freedom.}
\item{r2}{ Coefficient of determination.}
\item{sr2}{ Residual variance.}
\item{y}{ Response.}
\item{fdataobj}{ Functional explanatory data.}
\item{mdist}{ Distance matrix between \code{x} and \code{newx}.}
\item{Ker}{ Asymmetric kernel used.}
\item{h.opt}{ smoothing parameter or bandwidth.}
}
\references{
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional data analysis.}
Springer Series in Statistics, New York. \cr

Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University Press, 1994.
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:   \code{\link{fregre.np.cv}}, \code{\link{summary.fregre.fd}} and \code{\link{predict.fregre.fd}}  .\cr
Alternative method: \code{\link{fregre.basis}},cand  \code{\link{fregre.pc}}.

}
\examples{

data(tecator)
absorp=tecator$absorp.fdata
ind=1:129
x=absorp[ind,]
y=tecator$y[ind,2]
res.np=fregre.np(x,y,Ker=AKer.epa)
summary.fregre.fd(res.np)
res.np2=fregre.np(x,y,Ker=AKer.tri)
summary.fregre.fd(res.np2)

# with other semimetrics.
res.pca1=fregre.np(x,y,Ker=AKer.tri,metri=semimetric.pca,q=1)
summary.fregre.fd(res.pca1)
res.deriv=fregre.np(x,y,metri=semimetric.deriv)
summary.fregre.fd(res.deriv)
x.d2=fdata.deriv(x,nderiv=1,method="fmm",class.out='fdata')
res.deriv2=fregre.np(x.d2,y)
summary.fregre.fd(res.deriv2)
x.d3=fdata.deriv(x,nderiv=1,method="bspline",class.out='fdata')
res.deriv3=fregre.np(x.d3,y)
summary.fregre.fd(res.deriv3)
}
\keyword{regression}