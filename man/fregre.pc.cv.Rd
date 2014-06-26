\name{fregre.pc.cv}   
\Rdversion{1.1}
\alias{fregre.pc.cv}

\title{ Functional penalized PC regression with scalar response using selection of number of PC components}
\description{																						
Functional Regression with scalar response using selection of number of  (penalized) principal components PC through cross-validation. The algorithm selects the PC with best estimates the response. The selection is performed by cross-validation (CV) or  Model Selection Criteria (MSC). After is computing  functional regression using the best selection of principal components.
	}
\usage{
fregre.pc.cv(fdataobj, y, kmax=8,  lambda = 0, P = c(1, 0, 0), 
    criteria = "SIC",weights=rep(1,len=n),...) 
}
\arguments{
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
   \item{y}{ Scalar response with length \code{n}.}
  \item{kmax}{The number of components to include in the model.} 
  \item{lambda}{ Vector with the amounts of penalization. Default value is 0, i.e. no penalization is used.
  If \code{lambda=TRUE} the algorithm computes a sequence of lambda values.}
 \item{P}{ The vector of coefficients to define the penalty matrix object. For example,  if \code{P=c(1,0,0)}, ridge regresion is computed and if \code{P=c(0,0,1)}, penalized regression is computed penalizing the second derivative (curvature). }
 \item{criteria}{ Type of cross-validation (CV) or  Model Selection Criteria (MSC) applied. Possible values are \emph{"CV"}, \emph{"AIC"}, \emph{"AICc"}, \emph{"SIC"},\emph{"SICc"}.} 	
   \item{weights}{ weights}
 \item{\dots}{ Further arguments passed to \code{\link{fregre.pc}} or \code{\link{fregre.pls}}}  
}
\details{ 
%%%%%%%%%%%
The algorithm selects the best principal components \code{pc.opt} from the first \code{kmax} PC and (optionally) the best penalized parameter \code{lambda.opt} from a sequence of non-negative numbers \code{lambda}. \cr

If \code{kmax} is a integer (by default and recomended) the procedure is as follows  (see example 1):
  \itemize{
		\item Calculate the best principal component (\emph{pc.order[1]}) between \code{kmax} by \code{\link{fregre.pc}}.
 \item Calculate the second-best principal component (\code{pc.order [2]}) between the \code{(kmax-1)} by \code{\link{fregre.pc}} and calculate the criteria value of the two principal  components.
 \item The process (point 1 and 2) is repeated until \code{kmax} principal component (\emph{pc.order[kmax]}).
 \item The process (point 1, 2 and 3) is repeated for each \code{lambda} value.
 \item The method selects the principal components (\code{pc.opt}=\code{pc.order[1:k.min]}) and (optionally) the lambda parameter with minimum MSC criteria.
 }
If \code{kmax} is a sequence of integer the procedure is as follows (see example 2):
	\itemize{
		\item 	The method selects the  best principal components  with minimum MSC criteria by stepwise regression using  \code{\link{fregre.pc}} in each step.
 \item The process (point 1) is repeated for each \code{lambda} value.
 \item The method selects the principal components (\code{pc.opt}=\code{pc.order[1:k.min]}) and (optionally) the lambda parameter with minimum MSC criteria.\cr
 }
Finally, is computing  functional PC regression between functional explanatory variable \eqn{X(t)} and scalar response \eqn{Y} using the best selection of PC \code{pc.opt} and ridge parameter \code{rn.opt}.  \cr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
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

%\item The generalized minimum description length (gmdl) criteria: \cr
%
% \eqn{gmdl(k_n)=log \left[ \frac{1}{n-k_n}\sum_{i=1}^{n}{\Big(y_i-\hat{y}_i\Big)^2} \right] +K_n log \left(\frac{(n-k_n)\sum_{i=1}^{n}\hat{y}_i^2}{{\sum_{i=1}^{n}\Big(y_i-\hat{y}_i\Big)^2} }\right)+log(n) }
 
 %{MSC(k_n)=log [ 1/(n-k_n) \sum_(i=1:n){ (y_i- < X_i , \beta_(i,k_n) > )^2} ] +p_n k_n/n }
%\item The rho criteria:  \eqn{rho(k_n)=log \left[ \frac{1}{n-k_n}\sum_{i=1}^{n}\left(\frac{y_i-\hat{y}_i}{1-H_{ii}} \right)^2\right]}	
} 
where \code{criteria} is an argument that controls the type of validation used in the selection of the smoothing parameter \code{kmax}\eqn{=k_n} and penalized parameter \code{lambda}\eqn{=\lambda}.\cr
}
\value{
Return:\cr
\item{fregre.pc}{ Fitted regression object by the best (\code{pc.opt}) components.}
\item{pc.opt}{ Index of PC components selected.}
\item{MSC.min}{ Minimum Model Selection Criteria (MSC) value for the (\code{pc.opt} components.}
\item{MSC}{ Minimum  Model Selection Criteria (MSC) value for \code{kmax} components.}
}
\note{\code{criteria=``CV''} is not recommended: time-consuming.}

\author{ Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es} }
\seealso{
See also as:\code{\link{fregre.pc}} .
}
\references{
Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
}
\examples{

\dontrun{
data(tecator)
x<-tecator$absorp.fdata[1:129]
y<-tecator$y$Fat[1:129]
# no penalization
 res.pc1=fregre.pc.cv(x,y,8)
# 2nd derivative penalization
 res.pc2=fregre.pc.cv(x,y,8,lambda=TRUE,P=c(0,0,1))
#Ridge regression
res.pc3=fregre.pc.cv(x,y,1:8,lambda=TRUE,P=1) 
}
}
\keyword{regression}
