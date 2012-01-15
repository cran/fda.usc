\name{fregre.ppc.cv}   
\Rdversion{1.1}
\alias{fregre.ppc.cv}
\alias{fregre.ppls.cv}
\title{ Functional penalized PC (or PLS) regression with scalar response using selection of number of PC (or PLS) components}
\description{																						
Functional Regression with scalar response using selection of number of  penalized principal componentes PPC(or partial least squares components PPLS)  through cross-validation. The algorithm selects the PPLS components with best estimates the response. The selection is performed by cross-validation (CV) or  Model Selection Criteria (MSC). After is computing  functional regression using the best selection of PPC (or PPLS) components.
	}
\usage{
fregre.ppc.cv(fdataobj, y, kmax=8,  lambda = 0, P = c(0, 0, 1), 
    criteria = "SIC", ...) 

fregre.ppls.cv(fdataobj, y, kmax=8,  lambda = 0, P = c(0, 0, 1), 
    criteria = "SIC", ...) 
}
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
   \item{y}{ Scalar response with length \code{n}.}
  \item{kmax}{The number of components to include in the model.} 
  \item{lambda}{ Vector with the amounts of penalization. Default value is 0, i.e. no penalization is used.
  If \code{lambda=TRUE} the algorithm computes a sequence of lambda values.}
 \item{P}{ If \code{P} is a vector:  \code{P} are coefficients to define the penalty matrix object. By default \code{P=c(0,0,1)} penalize the second derivative (curvature) or acceleration.
  If \code{P} is a matrix: P is the penalty matrix object.} 
 \item{criteria}{ Type of cross-validation (CV) or  Model Selection Criteria (MSC) applied. Possible values are \emph{"CV"}, \emph{"AIC"}, \emph{"AICc"}, \emph{"SIC"}.} 	
 \item{\dots}{ Further arguments passed to \code{\link{fregre.ppc}} or \code{\link{fregre.ppls}}}  
}
\details{
The algorithm is as follows: 
	\itemize{
		\item 	Select the  bests components (\code{pc.opt} or \code{pls.opt}) with minimum MSC criteria by stepwise regression using  \code{\link{fregre.ppc}}  or \code{\link{fregre.ppls}} in each step. 
 \item Fit the functional PPLS regression between  \eqn{\tilde{X}(t)} and  \eqn{Y} using the best selection of FPLS components \code{pls.opt}.  
 }
For more details in estimation process see \code{\link{fregre.ppc}} or \code{\link{fregre.ppls}}.\cr
 
The criteria selection is done by cross-validation (CV) or Model Selection Criteria (MSC). 
\itemize{
	\item Predictive Cross-Validation:  \eqn{PCV(k_n)=\frac{1}{n}\sum_{i=1}^{n}{\Big(y_i -\hat{y}_{(-i,k_n)}	 \Big)^2}}{PCV(k_n)=1/n \sum_(i=1:n)  (y_i - \hat{y}_{-i})^2},\cr 
	\code{criteria}=``CV''
	
	\item Model Selection Criteria:  \eqn{MSC(k_n)=log \left[ \frac{1}{n}\sum_{i=1}^{n}{\Big(y_i-\hat{y}_i\Big)^2} \right] +p_n\frac{k_n}{n} }{MSC(k_n)=log [ 1/n \sum_(i=1:n){ (y_i- \hat{y}_i )^2} ] +p_n k_n/n } \cr
	   \eqn{p_n=\frac{log(n)}{n}}{p_n=log(n)/n}, \code{criteria}=``SIC'' (by default)\cr
	   \eqn{p_n=\frac{log(n)}{n-k_n-2}}{p_n=log(n)/(n-k_n-2)}, \code{criteria}=``SICc''\cr
	   \eqn{p_n=2}, \code{criteria}=``AIC''\cr
	   \eqn{p_n=\frac{2n}{n-k_n-2}}{p_n=2n/(n-k_n-2)}, \code{criteria}=``AICc''\cr
	   \eqn{p_n=\frac{2log(log(n))}{n}}{p_n=2log(log(n))/(n)}, \code{criteria}=``HQIC''\cr 

\item The generalized minimum description length (gmdl) criteria: \cr

 \eqn{gmdl(k_n)=log \left[ \frac{1}{n-k_n}\sum_{i=1}^{n}{\Big(y_i-\hat{y}_i\Big)^2} \right] +K_n log \left(\frac{(n-k_n)\sum_{i=1}^{n}\hat{y}_i^2}{{\sum_{i=1}^{n}\Big(y_i-\hat{y}_i\Big)^2} }\right)+log(n) }
 
 %{MSC(k_n)=log [ 1/(n-k_n) \sum_(i=1:n){ (y_i- < X_i , \beta_(i,k_n) > )^2} ] +p_n k_n/n }
\item The rho criteria:  \eqn{rho(k_n)=log \left[ \frac{1}{n-k_n}\sum_{i=1}^{n}\left(\frac{y_i-\hat{y}_i}{1-H_{ii}} \right)^2\right]}	
}  where \code{criteria} is an argument that controls the type of validation used in the selection of the smoothing parameter \code{kmax}\eqn{=k_n} and penalized parameter \code{lambda}\eqn{=\lambda}.\cr

\code{criteria=``CV''} is not recommended: time-consuming.
}
\value{
Return:\cr
\item{pls.opt}{ Index of PC or PLS components selected.}
\item{MSC.min}{ Minimum Model Selection Criteria (MSC) value for the (\code{pc.opt} or \code{pls.opt}) components.}
\item{MSC}{ Minimum  Model Selection Criteria (MSC) value for \code{kmax} components.}
\item{fregre.ppls}{ Fitted regression object by the best (\code{pc.opt} or \code{pls.opt}) components.}
}
\author{ Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es} }
\seealso{
See also as:   \code{\link{fregre.ppls}} and  \code{\link{fregre.ppc}} .
}
\references{
Preda C. and Saporta G. \emph{PLS regression on a stochastic process}. Comput. Statist. Data Anal. 48 (2005): 149{-}158.

Kraemer, N., Sugiyama M. (2011). \emph{The Degrees of Freedom of Partial Least Squares Regression}. Journal of the American Statistical Association. Volume 106, 697-705.

Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
}
\examples{

data(tecator)
x<-tecator$absorp.fdata[1:129]
y<-tecator$y$Fat[1:129]
# fregre.ppls.cv(x,y,8)
# fregre.ppls.cv(x,y,8,lambda=TRUE)

}
\keyword{regression}
