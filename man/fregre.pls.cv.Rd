\name{fregre.pls.cv}
\Rdversion{1.1}
\alias{fregre.pls.cv}
\title{ Functional PLS regression with scalar response using selection of number of  PLS components}
\description{																						
Functional Regression with scalar response using selection of number of partial leas squares (PLS) components through cross-validation. The algorithm selects the PLS components with best estimates the response. The selection is performed by cross-validation (CV) or  Model Selection Criteria (MSC). After is computing  functional regression using the best selection of PLS components.
	}
\usage{fregre.pls.cv(fdataobj, y, kmax=8, criteria = "SICc",...)}
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
  \item{kmax}{The number of PLS components to include in the model.}
  \item{y}{ Scalar response with length \code{n}.}
 \item{criteria}{ Type of cross-validation (CV) or  Model Selection Criteria (MSC) applied. Possible values are \emph{"CV"}, \emph{"AIC"}, \emph{"AICc"}, \emph{"SIC"}.} 	
 \item{\dots}{ Further arguments passed to \code{\link{fregre.pls}}}  
}
\details{
The algorithm is as follows: 
	\itemize{
		\item 	Select the  bests PLS components \code{pls.opt} with minimum MSC criteria by stepwise regression using  \code{\link{fregre.pls}} in each step. 
 \item Fit the functional PLS regression between  \eqn{X(t)} and  \eqn{Y} using the best selection of FPLS components \code{pls.opt}.  
 }
For more details in estimation process see \code{\link{fregre.pls}}.\cr
 
The criteria selection is done by cross-validation (CV) or Model Selection Criteria (MSC). 
\itemize{
	\item Predictive Cross-Validation:  \eqn{PCV(k_n)=\frac{1}{n}\sum_{i=1}^{n}{\Big(y_i -\big< X_i,\hat{\beta}_{(-i,k_n)} \big>\Big)^2}}{PCV(k_n)=1/n \sum_(i=1:n)  (y_i - < X_i , \beta.est_(-i,k_n) >)^2},\cr 
	\code{criteria}=``CV''
	
	\item Model Selection Criteria:  \eqn{MSC(k_n)=log \left[ \frac{1}{n}\sum_{i=1}^{n}{\Big(y_i-\big< X_i,\hat{\beta}_{(i,k_n)} \big>\Big)^2} \right] +p_n\frac{k_n}{n} }{MSC(k_n)=log [ 1/n \sum_(i=1:n){ (y_i- < X_i , \beta_(i,k_n) > )^2} ] +p_n k_n/n } \cr
	   \eqn{p_n=\frac{log(n)}{n}}{p_n=log(n)/n}, \code{criteria}=``SIC'' (by default)\cr
	   \eqn{p_n=\frac{log(n)}{n-k_n-2}}{p_n=log(n)/(n-k_n-2)}, \code{criteria}=``SICc''\cr
	   \eqn{p_n=2}, \code{criteria}=``AIC''\cr
	   \eqn{p_n=\frac{2n}{n-k_n-2}}{p_n=2n/(n-k_n-2)}, \code{criteria}=``AICc''
}
% falta poner CV criteria
where \code{criteria} is an argument that controls the type of validation used in the selection of the smoothing parameter \code{kmax}\eqn{=k_n}.\cr

\code{criteria=``CV''} is not recommended: time-consuming.
}
\value{
Return:\cr
\item{pls.opt}{ Index of PLS components selected.}
\item{MSC.min}{ Minimum Model Selection Criteria (MSC) value for the code{pls.opt} components.}
\item{MSC}{ Minimum  Model Selection Criteria (MSC) value for \code{kmax} components.}
\item{fregre.pls}{ Fitted regression object by the best \code{pls.opt} components.}
}
\author{ Febrero-Bande, M. and Oviedo de la Fuente, M. }
\seealso{
See Also as:   \code{\link{fregre.pls}}, \code{\link{summary.fregre.fd}} and \code{\link{predict.fregre.fd}}.\cr
Alternative method: \code{\link{fregre.pc}}, \code{\link{fregre.basis}} and   \code{\link{fregre.np}}.
}
\references{
Preda C. and Saporta G. \emph{PLS regression on a stochastic process}. Comput. Statist. Data Anal. 48 (2005): 149{-}158.\cr
Kraemer, N., Sugiyama M. (2011). \emph{The Degrees of Freedom of Partial Least Squares Regression}. Journal of the American Statistical Association. Volume 106, 697-705.
}
\examples{

data(tecator)
x<-tecator$absorp.fdata[1:129]
y<-tecator$y$Fat[1:129]
fregre.pls.cv(x,y,8)

}
\keyword{regression}
