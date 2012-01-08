\name{fregre.pc.cv}
\Rdversion{1.1}
\alias{fregre.pc.cv}
%\alias{fregre.combn}
\title{Vaidation criteria for Functional Principal Component (and Ridge)  Regression using selection of number of Principal Components}
\description{																						
Functional regression with scalar response using selection of number of principal components through cross-validation. The algorithm selects the principal components with best estimates the response. The	select is performed by cross-validation (CV) or  model selection criteria (MSC).  After is computing  functional regression using the best selection of principal components (PC).
	}
\usage{
fregre.pc.cv(fdataobj, y, kmax = 8, rn = 0, criteria = "SIC",weights= rep(1,n), ...) 
}
\arguments{
   \item{fdataobj}{ \code{\link{fdata}} class object or \code{fdata.comp} class object created by \code{\link{create.pc.basis}} function.} 
	\item{y}{ Scalar response with length \code{n}.}
	\item{kmax}{ Maximum number of principal comoponents.}
	\item{criteria}{ Type of cross-validation (CV) or  Model Selecttion Criteria (MSC) applied. Possible values are \emph{"CV"}, \emph{"AIC"}, \emph{"AICc"}, \emph{"SIC"} and \emph{"SICc"}.} 	
		\item{rn}{Ridge parameter. Default value is \code{rn=0}, i.e. no penalization is used. If =TRUE, a sequence of values is provided}
		\item{weights}{ weights}
 \item{\dots}{ Further arguments passed to \code{\link{fregre.pc}}}  
}
\details{
The algorithm selects the best principal components \code{pc.opt} from the first \code{kmax} PC and (optionally) the best ridge parameter \code{rn.opt} from a sequence of non-negative numbers \code{rn}. \cr

If \code{kmax} is a integer (by default and recomended) the procedure is as follows  (see example 1):
	\itemize{
		\item Calculate the best principal component (\emph{pc.order[1]}) between \code{kmax} by \code{\link{fregre.pc}}.
 \item Calculate the second-best principal component (\code{pc.order [2]}) between the \code{(kmax-1)} by \code{\link{fregre.pc}} and calculate the criteria value of the two principal  components.
 \item The process (point 1 and 2) is repeated until \code{kmax} principal component (\emph{pc.order[kmax]}).
 \item The process (point 1, 2 and 3) is repeated for each \code{rn} value.
 \item The method selects the principal components (\code{pc.opt}=\code{pc.order[1:k.min]}) and (optionally) the ridge parameter with minimum MSC criteria.
 }
If \code{kmax} is a sequence of integer the procedure is as follows (see example 2):
	\itemize{
		\item 	The method selects the  best principal components  with minimum MSC criteria by stepwise regression using  \code{\link{fregre.pc}} in each step.
 \item The process (point 1) is repeated for each \code{rn} value.\cr 
 \item The method selects the principal components (\code{pc.opt}=\code{pc.order[1:k.min]}) and (optionally) the ridge parameter with minimum MSC criteria.\cr
 }
Finally, is computing  functional PC regression between functional explanatory variable \eqn{X(t)} and scalar response \eqn{Y} using the best selection of PC \code{pc.opt} and ridge parameter \code{rn.opt}.  \cr

%For more details in estimation proces see \code{\link{fregre.pc}}.\cr
 
The criteria selection is done by cross-validation (CV) or Model Selection Criteria (MSC). 
\itemize{
	\item Predictive Cross-Validation:  \eqn{PCV(k_n)=\frac{1}{n}\sum_{i=1}^{n}{\Big(y_i -\big< X_i,\hat{\beta}_{(-i,k_n)} \big>\Big)^2}}{PCV(k_n)=1/n \sum_(i=1:n)  (y_i - < X_i , \beta.est_(-i,k_n) >)^2},\cr 
	\code{criteria}=``CV''
	
	\item Model Selection Criteria:  \eqn{MSC(k_n)=log \left[ \frac{1}{n}\sum_{i=1}^{n}{\Big(y_i-\big< X_i,\hat{\beta}_{(i,k_n)} \big>\Big)^2} \right] +p_n\frac{k_n}{n} }{MSC(k_n)=log [ 1/n \sum_(i=1:n){ (y_i- < X_i , \beta_(i,k_n) > )^2} ] +p_n k_n/n } \cr
	   \eqn{p_n=\frac{log(n)}{n}}{p_n=log(n)/n}, \code{criteria}=``SIC'' (by default)\cr
	   \eqn{p_n=\frac{log(n)}{n-k_n-2}}{p_n=log(n)/(n-k_n-2)}, \code{criteria}=``SICc''\cr
	   \eqn{p_n=2}, \code{criteria}=``AIC''\cr
	   \eqn{p_n=\frac{2n}{n-k_n-2}}{p_n=2n/(n-k_n-2)}, \code{criteria}=``AICc''\cr
	   \eqn{p_n=\frac{2log(log(n))}{n}}{p_n=2log(log(n))/(n)}, \code{criteria}=``HQIC''\cr 

\item The generalized minimum description length (gmdl) criteria: \cr

 \eqn{gmdl(k_n)=log \left[ \frac{1}{n-k_n}\sum_{i=1}^{n}{\Big(y_i-\big< X_i,\hat{\beta}_{(i,k_n)} \big>\Big)^2} \right] +K_n log \left(\frac{\sum_{i=1}^{n}\big< X_i,\hat{\beta}_{(i,k_n)} \big>^2}{\frac{\sum_{i=1}^{n}\Big(y_i-\big< X_i,\hat{\beta}_{(i,k_n)} \big>\Big)^2}{n-k_n} }\right)+log(n) }
 
%{MSC(k_n)=log [ 1/(n-k_n) \sum_(i=1:n){ (y_i- < X_i , \beta_(i,k_n) > )^2} ] +p_n k_n/n }
\item The rho criteria:  \eqn{rho(k_n)=log \left[ \frac{1}{n-k_n}\sum_{i=1}^{n}\left(\frac{y_i-\big< X_i,\hat{\beta}_{(i,k_n)} \big>}{1-H_{ii}}\right)^2 \right]}	    
	    
%	    "SICc"=log(s2) + log(n) * ck/(n-ck-2),
%	     "HQIC"=log(s2) + 2*log(log(n)) * ck/n,
%     RSS<-sum(out$residuals^2)
             %s2 <- RSS/n
             %S0<-sqrt(RSS/(n-ck)) #sigmahat
             %FF<-sum(out$fitted.values^2)/(ck*RSS/(n-ck))
             % "gmdl"=(n*log(RSS/(n-ck))+ck*log(FF)+log(n))/2,
             % "rho"={A<-out$residuals;B<-1-diag(out$H)/n; D1<-(A/B)^2;sum(D1)})
}

  
where \code{criteria} is an argument that controls the type of validation used in the selection of the smoothing parameter \code{kmax}\eqn{=k_n}.\cr

\code{criteria=``CV''} is not recommended: time-consuming.
}
% incluir rn 
\value{
Return:\cr
\item{fregre.pc}{regression model fitted  with best selection of FPC \code{(pc.opt)} and ridge parameter \code{(rn.opt)}.}
\item{pc.opt}{ Index of principal components selected.}
\item{rn.opt}{ Ridge parameter selected.}
\item{PC.order}{Matrix whose columns correspond to the selected index of each principal component with the lowest MSC value by each ridge parameter value (in rows)}
\item{MSC.order}{ Model Selection Criteria (MSC) value for complete sequence principal components selected  by ridge parameter value in each iteration.}
}
\references{
Cai T., Hall P. 2006. \emph{Prediction in functional linear regression}. Annals of Statistics 34: 2159{-}2179.

Cardot H, Ferraty F, Sarda P. 1999. \emph{Functional linear model}. Statistics and Probability Letters 45: 11{-}22.

Hall P, Hosseini{-}Nasab M. 2006. \emph{On properties of functional principal components analysis}. Journal of the Royal Statistical Society B 68: 109{-}126.

Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{
See Also as:   \code{\link{fregre.pc}}, \code{\link{summary.fregre.fd}} and \code{\link{predict.fregre.fd}}.\cr
Alternative method:  \code{\link{fregre.basis.cv}} and   \code{\link{fregre.np.cv}}.
}
\examples{
#################################################################
data(aemet)
names(aemet)
x=aemet$temp
y=apply(aemet$logprec$data,1,sum)

# Example 1: criteria selects the best FPC (recommended)
# kmax is a integer (maximum FPC to consider)
out1=fregre.pc.cv(x,y,8,criteria="AICc")

 
## Example 2: Three equivalent outputs
## using the functional data
# out2=fregre.pc.cv(x,y,8,criteria="SICc") 

## using the FPC components (recommended)
# pc<-fdata2pc(x,8) 											
# out111=fregre.pc.cv(pc,y,criteria="SICc")

## using the FPC components
# out11=fregre.pc.cv(out1$fregre.pc$fdata.comp,y,criteria="SICc") 


## Example 3: stepwise selection criteria (unrecommended)
## kmax is a sequence of integer
# out3=fregre.pc.cv(x,y,1:8,criteria="AICc")
# out4=fregre.pc.cv(x,y,1:8,criteria="SICc")

}
\keyword{regression}
