\name{min.np}
\Rdversion{1.1}
\alias{min.np}
\title{ Smoothing of functional data using nonparametric kernel estimation }
\description{
Smoothing of functional data using nonparametric kernel estimation with cross-validation (CV) or generalized cross-validation (GCV) methods. 
}
\usage{
%min.np(fdataobj,h=NULL,W=diag(1,nrow=np,ncol=np),Ker=Ker.norm,
%type.CV=GCV.S,type.S=S.NW,par.CV=list(trim=0,draw=FALSE),...)
\method{min}{np}(fdataobj,h=NULL,W=NULL,Ker=Ker.norm,
type.CV=GCV.S,type.S=S.NW,par.CV=list(trim=0,draw=FALSE),verbose=FALSE,\dots)
}

\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
  \item{h}{ Smoothing parameter or bandwidth.}
  \item{W}{ Matrix of weights.}
  \item{Ker}{ Type of kernel used, by default normal kernel.}
  \item{type.CV}{ Type of cross-validation. By default generalized cross-validation (GCV) method. Possible values are \emph{GCV.S} and \emph{CV.S}}
  \item{type.S}{ Type of smothing matrix \code{S}. By default \code{S} is calculated by Nadaraya-Watson kernel estimator (\code{S.NW}). Possible values are \emph{S.NW} and \emph{S.KNN}}
  \item{par.CV}{ List of parameters for type.CV: \code{trim}, the alpha of the trimming and \code{draw=TRUE}.}	    
  \item{verbose}{ If \code{TRUE} information about  GCV values and input parameters is printed. Default is \code{FALSE}.}
  \item{\dots}{ Further arguments passed to or from other methods. Arguments to be passed  for kernel method.}
}

\details{  

Calculate the minimum GCV for a vector of values of the smoothing parameter \code{h}.  \cr
Nonparametric smoothing is performed by the kernel function. The type of kernel to use with the parameter \code{Ker} and the type of smothing matrix \code{S} to use with the parameter \code{type.S} can be selected by the user, see function \code{\link{Kernel}}.\cr
W is the matrix of weights of the discretization points.
\cr

}
\value{
Returns GCV or CV values calculated for input parameters.\cr
    \item{gcv}{ GCV  or CV for a vector of values of the smoothing parameter \code{h}  }
    \item{fdataobj}{ \code{\link{fdata}} class object.}  
    \item{fdata.est}{ Estimated \code{fdata} class object.}
    \item{h.opt}{ \code{h} value that minimizes CV or GCV method.}
    \item{S.opt}{ Smoothing matrix for the minimum  CV or GCV method.}
    \item{gcv.opt}{ Minimum  of CV or GCV method.}
	  \item{h}{ Smoothing parameter or bandwidth.}    
}
\references{
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional data analysis.}
Springer Series in Statistics, New York. 

Wasserman, L. \emph{All of Nonparametric Statistics}. Springer Texts in Statistics, 2006. 

Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University Press, 1994.

Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}

}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{
See Also as  \code{\link{S.NW}}. \cr
Alternative method:  \code{\link{min.basis}}
}
\examples{

# Exemple, phoneme DATA
data(phoneme)
mlearn<-phoneme$learn[1:100]

out1<-min.np(mlearn,type.CV=CV.S,type.S=S.NW)
np<-ncol(mlearn)
# variance calculations
y<-mlearn
out<-out1
i<-1
z=qnorm(0.025/np)
fdata.est<-out$fdata.est
tt<-y[["argvals"]]
var.e<-Var.e(y,out$S.opt)
var.y<-Var.y(y,out$S.opt)
var.y2<-Var.y(y,out$S.opt,var.e)

# plot estimated fdata and point confidence interval
upper.var.e<-fdata.est[i,]-z*sqrt(diag(var.e))
lower.var.e<-fdata.est[i,]+z*sqrt(diag(var.e))
dev.new()
plot(y[i,],lwd=1,
ylim=c(min(lower.var.e$data),max(upper.var.e$data)),xlab="t")
lines(fdata.est[i,],col=gray(.1),lwd=1)
lines(fdata.est[i,]+z*sqrt(diag(var.y)),col=gray(0.7),lwd=2)
lines(fdata.est[i,]-z*sqrt(diag(var.y)),col=gray(0.7),lwd=2)
lines(upper.var.e,col=gray(.3),lwd=2,lty=2)
lines(lower.var.e,col=gray(.3),lwd=2,lty=2)
legend("bottom",legend=c("Var.y","Var.error"),
col = c(gray(0.7),gray(0.3)),lty=c(1,2))

}

\keyword{nonparametric} 