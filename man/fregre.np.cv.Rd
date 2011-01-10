\name{fregre.np.cv}
\Rdversion{1.1}
\alias{fregre.np.cv}
\title{ Cross-validation functional regression with scalar response using kernel estimation.}
\description{
 	Computes functional regression between functional explanatory variables and scalar response using asymetric kernel estimation by cross-validation method.
}
\usage{
fregre.np.cv(fdataobj,y,h=NULL,Ker=AKer.norm,metric=metric.lp,
type.CV = GCV.S,type.S=S.NW,par.CV=list(trim=0),...)
}
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
  \item{y}{ Scalar response with length \code{n}.}
  \item{h}{ Smoothing parameter or bandwidth.}
  \item{Ker}{ Type of asymmetric kernel used, by default asymetric normal kernel.}	  	
  \item{metric}{ Metric function, by default \code{\link{metric.lp}}.}
  \item{type.CV}{ Type of cross-validation. By default generalized cross-validation \code{\link{GCV.S}} method.} 
  \item{type.S}{ Type of smothing matrix \code{S}. By default \code{S} is calculated by Nadaraya-Watson kernel estimator (\code{S.NW}).}
  \item{par.CV}{ List of parameters for \code{type.CV}: \code{trim}, the alpha of the trimming and \code{draw=TRUE}.}	  
  \item{\dots}{ Arguments to be passed  for \code{\link{metric.lp}} o other metric function.}
}
\details{
The function estimates the value of smoothing parameter or the bandwidth through the cross validation methods: \code{\link{GCV.S}} or \code{\link{CV.S}}. It computes the distance between curves using the \code{\link{metric.lp}}, although any other semimetric could be used (see \code{\link{semimetric.basis}} or \code{\link{semimetric.NPFDA}} functions). 
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
\item{gcv}{  CV or GCV values.}
\item{h.opt}{ smoothing parameter or bandwidth that minimizes CV or GCV method.}
\item{h}{ Vector of smoothing parameter or bandwidth.}
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
See Also as:   \code{\link{fregre.np}}, \code{\link{summary.fregre.fd}} and \code{\link{predict.fregre.fd}}  .\cr
Alternative method: \code{\link{fregre.basis.cv}} and   \code{\link{fregre.np.cv}}.

}
\examples{
data(tecator)
absorp=tecator$absorp.fdata
ind=1:129
x=absorp[ind,]
y=tecator$y[ind,2]
Ker=AKer.tri
res.np=fregre.np.cv(x,y,Ker=Ker)
summary.fregre.fd(res.np)
res.np2=fregre.np.cv(x,y,type.CV=GCV.S,criteria="Shibata")
summary.fregre.fd(res.np2)

# with other semimetrics.
res.pca1=fregre.np.cv(x,y,Ker=Ker,metric=semimetric.pca,q=1)
summary.fregre.fd(res.pca1)
res.deriv=fregre.np.cv(x,y,Ker=Ker,metric=semimetric.deriv)
summary.fregre.fd(res.deriv)

x.d2=fdata.deriv(x,nderiv=1,method="fmm",class.out='fdata')
res.deriv2=fregre.np.cv(x.d2,y,Ker=Ker)
summary.fregre.fd(res.deriv2)
x.d3=fdata.deriv(x,nderiv=1,method="bspline",class.out='fdata')
res.deriv3=fregre.np.cv(x.d3,y,Ker=Ker)
summary.fregre.fd(res.deriv3)
}
\keyword{regression}