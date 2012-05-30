\name{fregre.plm}
\Rdversion{1.1}
\alias{fregre.plm}
\title{ Semi-functional partially linear model with scalar response.}
\description{
	Computes functional regression between functional (and non functional) explanatory variables and scalar response using asymmetric kernel estimation.\cr
	
An extension of the non-parametric functional regression models is the semi-functional partial linear model proposed in Aneiros-Perez and Vieu (2005). This model uses a non-parametric kernel procedure as that described in \code{\link{fregre.np}}. The output \eqn{y}
is scalar. A functional covariate \eqn{X} and a multivariate non functional covariate \eqn{Z} are considered.

\deqn{y =\emph{r(X)}+\sum_{j=1}^{p}{Z_j\beta_j}+\epsilon}{y = r(X)+\sum_(j=1:p) Z_j \beta_j+\epsilon}

The unknown smooth real function \eqn{r} is  estimated by means of \deqn{\hat{r}_{h}(X)=\sum_{i=1}^{n}{w_{n,h}(X,X_{i})(Y_{i}-Z_{i}^{T}\hat{\beta}_{h})}}{\hat{r}_{h}(X)=\sum_(i=1:n) w_{n,h}(X,X_i) (Y_i - Z_i^{T} \beta.est_h)} where \eqn{W_h} is the weight function:

 \eqn{w_{n,h}(X,X_{i})=\frac{K(d(X,X_i)/h)}{\sum_{j=1}^{n}K(d(X,X_j)/h)}}{w_{n,h}(X,X_i)= K( d( X , X_i )/h ) / \sum_(j=1:n) K( d( X, X_j )/h )} with smoothing parameter \eqn{h}, an asymmetric kernel \eqn{K} and a metric or semi-metric \eqn{d}.  In \code{fregre.plm()} by default \eqn{W_h} is a functional version of the Nadaraya-Watson-type weights (\code{type.S=S.NW}) with asymmetric normal kernel (\code{Ker=AKer.norm}) in \eqn{L_2}\cr (\code{metric=metric.lp} with \code{p=2}). The unknown parameters \eqn{\beta_j}  for  the multivariate non functional covariates are  estimated by means of \eqn{\hat{\beta}_j=(\tilde{Z}_{h}^{T}\tilde{Z}_{h})^{-1}\tilde{Z}_{h}^{T}\tilde{Z}_{h}}{\beta.est_j=(Z_h'Z_h)^{-1} Z_h^{T}Z_h} where \eqn{\tilde{Z}_{h}=(I-W_{h})Z}{Z_h=(I-W_h)Z} with the smoothing parameter \eqn{h}. The errors \eqn{\epsilon} are independent, with zero mean, finite variance \eqn{\sigma^2} and \eqn{E[\epsilon|Z_1,\ldots,Z_p,X(t)]=0}{E[\epsilon|Z_1,...,Z_p,X(t)]=0}.\cr
	
}
\usage{
fregre.plm(formula,data,h=NULL,Ker=AKer.norm,metric=metric.lp,
type.CV = GCV.S,type.S=S.NW,par.CV=list(trim=0,draw=FALSE),
par.S=list(w=1),...)}
\arguments{
 \item{formula}{ an object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under \code{Details}.}
  \item{data}{ List that containing the variables in the model. }
	\item{Ker}{ Type of asymmetric kernel used, by default asymmetric normal kernel.}	
	\item{h}{ Smoothing parameter or bandwidth.}
  \item{metric}{ Metric function, by default \code{\link{metric.lp}}.}
  \item{type.CV}{ Type of cross-validation. By default generalized cross-validation \code{\link{GCV.S}} method.} 
  \item{type.S}{ Type of smothing matrix \code{S}. By default \code{S} is calculated by Nadaraya-Watson kernel estimator (\code{S.NW}).}   
  \item{par.CV}{ List of parameters for \code{type.CV}: \code{trim}, the alpha of the trimming and \code{draw=TRUE}.}	  
   \item{par.S}{ List of parameters for \code{type.S}: \code{w}, the weights.}	      
   \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{ 
The first item in the \code{data} list  is called \emph{"df"} and  is a data frame  with the response and non functional explanatory variables, as \code{link{lm}}. If non  functional data into the formula then \code{\link{lm}} regression is performed.\cr
Functional variable (\code{fdata} or \code{fd} class) is introduced in the second item in the \code{data} list.  If only  functional variable into the formula then \code{\link{fregre.np.cv}} is performed.\cr
 
The function estimates the value of smoothing parameter or the bandwidth \code{h} through  Generalized Cross-validation \code{GCV} criteria.
It computes the distance between curves using the \code{\link{metric.lp}}, although you can also use other metric function. \cr
Different asymmetric kernels can be used, see   \code{\link{Kernel.asymmetric}}.\cr
}
\value{
\item{call}{ The matched call.}
\item{fitted.values}{  Estimated scalar response.}
\item{residuals}{ \code{y} minus \code{fitted values}.}
\item{df}{ The residual degrees of freedom.}
\item{H}{  Hat matrix.}
\item{r2}{ Coefficient of determination.}
\item{sr2}{ Residual variance.}
\item{y}{ Scalar response.}
\item{fdataobj}{ Functional explanatory data.}
\item{XX}{ Non functional explanatory data.}
\item{mdist}{ Distance matrix between curves.}
\item{betah}{ beta coefficient estimated}
\item{data}{  List that containing the variables in the model.}
\item{Ker}{ Asymmetric kernel used.}
\item{h.opt}{ Value that minimizes CV or GCV method.}
\item{h}{ Smoothing parameter or bandwidth.}\item{data}{ List that containing the variables in the model.}
\item{gcv}{  GCV values.}
\item{formula}{ formula.} 
}
\references{
Aneiros-Perez G. and Vieu P. (2005). \emph{Semi-functional partial linear regression}.  Statistics {&} Probability Letters, 76:1102-1110.\cr
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional data analysis.} Springer Series in Statistics, New York. \cr
Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University Press, 1994.
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:  \code{\link{predict.fregre.plm}} and  \code{\link{summary.fregre.fd}}\cr
Alternative methods: \code{\link{fregre.lm}}, \code{\link{fregre.np}} and   \code{\link{fregre.np.cv}}
}
\examples{

data(tecator)
x=tecator$absorp.fdata[1:129]
dataf=tecator$y[1:129,]

f=Fat~Water+x
ldata=list("df"=dataf,"x"=x)
res.plm=fregre.plm(f,ldata)
summary(res.plm)

# with 2nd derivative of functional data
x.fd=fdata.deriv(x,nderiv=2)
f2=Fat~Water+x.fd
ldata2=list("df"=dataf,"x.fd"=x.fd)
res.plm2=fregre.plm(f2,ldata2)
summary(res.plm2)

}
\keyword{regression}