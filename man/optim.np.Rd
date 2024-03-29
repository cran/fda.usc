% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optim.np.r
\name{optim.np}
\alias{optim.np}
\alias{min.np}
\title{Smoothing of functional data using nonparametric kernel estimation}
\usage{
optim.np(
  fdataobj,
  h = NULL,
  W = NULL,
  Ker = Ker.norm,
  type.CV = GCV.S,
  type.S = S.NW,
  par.CV = list(trim = 0, draw = FALSE),
  par.S = list(),
  correl = TRUE,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{fdataobj}{\code{\link{fdata}} class object.}

\item{h}{Smoothing parameter or bandwidth.}

\item{W}{Matrix of weights.}

\item{Ker}{Type of kernel used, by default normal kernel.}

\item{type.CV}{Type of cross-validation. By default generalized
cross-validation (GCV) method. Possible values are \emph{GCV.S} and
\emph{CV.S}}

\item{type.S}{Type of smothing matrix \code{S}. By default \code{S} is
calculated by Nadaraya-Watson kernel estimator (\code{S.NW}). Possible
values are \code{\link{S.KNN}},  \code{\link{S.LLR}}, 
 \code{\link{S.LPR}} and  \code{\link{S.LCR}}.}

\item{par.CV}{List of parameters  for type.CV: \code{trim}, the alpha of the
trimming and \code{draw=TRUE}}

\item{par.S}{List of parameters for \code{type.S}:
\code{tt} for argvals,  \code{h} for bandwidth,  
 \code{Ker} for kernel, etc.}

\item{correl}{logical. If \code{TRUE} the bandwidth parameter \code{h} is computed following the 
procedure described for De  Brabanter et al. (2018). (option avalaible since v1.6.0 version)}

\item{verbose}{If \code{TRUE} information about GCV values and input
parameters is printed. Default is \code{FALSE}.}

\item{\dots}{Further arguments passed to or from other methods. Arguments to
be passed for kernel method.}
}
\value{
Returns GCV or CV values calculated for input parameters.
\itemize{
\item \code{gcv}{ GCV or CV for a vector of values of the smoothing parameter
\code{h} } 
\item \code{fdataobj}{ \code{\link{fdata}} class object.}
\item \code{fdata.est}{ Estimated \code{fdata} class object.} 
\item \code{h.opt}{ \code{h} value that minimizes CV or GCV method.} 
\item \code{S.opt}{ Smoothing matrix for the minimum CV or GCV method.} 
\item \code{gcv.opt}{ Minimum of CV or GCV method.} 
\item \code{h}{ Smoothing parameter or bandwidth.}
}
}
\description{
Smoothing of functional data using nonparametric kernel
estimation with cross-validation (CV) or generalized cross-validation 
(GCV) methods.
}
\details{
Calculate the minimum GCV for a vector of values of the smoothing parameter
\code{h}.
Nonparametric smoothing is performed by the kernel function.
The type of kernel to use with the parameter \code{Ker} and the type of
smothing matrix \code{S} to use with the parameter \code{type.S} can be
selected by the user, see function \code{\link{Kernel}}.
W is the matrix of weights of the discretization points.
}
\note{
min.np deprecated.
}
\examples{
\dontrun{
# Exemple, phoneme DATA
data(phoneme)
mlearn<-phoneme$learn[1:100]

out1<-optim.np(mlearn,type.CV=CV.S,type.S=S.NW)
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

}
\references{
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional
data analysis.} Springer Series in Statistics, New York.

Wasserman, L. \emph{All of Nonparametric Statistics}. Springer Texts in
Statistics, 2006.

Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University
Press, 1994. 

De Brabanter, K., Cao, F., Gijbels, I., Opsomer, J. (2018). Local polynomial regression with correlated errors in random design and
unknown correlation structure.  \emph{Biometrika}, 105(3), 681-69.

Febrero-Bande, M., Oviedo de la Fuente, M. (2012). Statistical Computing in
Functional Data Analysis: The R Package fda.usc. \emph{Journal of
Statistical Software}, 51(4), 1-28. \url{https://www.jstatsoft.org/v51/i04/}
}
\seealso{
Alternative method:  \code{\link{optim.basis}}
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente
\email{manuel.oviedo@udc.es}
}
\keyword{nonparametric}
