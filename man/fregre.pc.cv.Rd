\name{fregre.pc.cv}
\Rdversion{1.1}
\alias{fregre.pc.cv}
\alias{fregre.combn}
\title{Vaidation criteria for Functional Principal Component (and Ridge)  Regression using selection of number of principal components}
\description{																						
Functional Regression with scalar response using selection of number of principal components through cross-validation. The algorithm selects the  principal components with best estimates the response. The	select is performed by cross-validation (CV) or  Model Selecttion Criteria (MSC).  After is computing  functional regression using the best selection of Principal Components.
	}
\usage{
fregre.pc.cv(fdataobj, y, kmax = 8, rn = 0, criteria = "SIC", ...) 

}
\arguments{
   \item{fdataobj}{ \code{\link{fdata}} class object or \code{fdata.comp} class object created by \code{\link{create.pc.basis}} function.} 
	\item{y}{ Scalar response with length \code{n}.}
	\item{kmax}{ Maximum number of principal comoponents.}
	\item{criteria}{ Type of cross-validation (CV) or  Model Selecttion Criteria (MSC) applied. Possible values are \emph{"CV"}, \emph{"AIC"}, \emph{"AICc"}, \emph{"SIC"} and \emph{"SICc"}.} 	
		\item{rn}{Ridge parameter. Default value is \code{rn=0}, i.e. no penalization is used. If =TRUE, a sequence of values is provided}
 \item{\dots}{ Further arguments passed to \code{\link{fregre.pc}}}  
}
\details{
The algorithm selects the best principal components (PC) from the first \code{kmax} PC and (optionally) the best ridge parameter from a sequence of non-negative numbers. The procedure is as follows
 1.-Calculate the best principal component (\emph{pc.order[1]}) between \code{kmax} first.\cr
 2.-Calculate the second-best principal component (\emph{pc.order [2]}) between the \code{(kmax-1)} and calculate the criteria value of the two principal  components.\cr
 3.-The process (point 1 and 2) is repeated until \code{kmax} principal component (\emph{pc.order[kmax]}).\cr
 4.-The process (point 1, 2 and 3) is repeated for each \code{rn} value.\cr 
 The method selects the principal components (\code{pc.opt}=\code{pc.order[1:k.min]}) and (optionally) the ridge parameter with minimum MSC criteria.\cr
 5.-Finally, is computing  functional PC regression between functional explanatory variable and scalar response using the best selection of PC and ridge parameter.  For more details see \code{\link{fregre.pc}}.
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
Cai T., Hall P. 2006. \emph{Prediction in functional linear regression}. Annals of Statistics 34: 2159{-}2179.\cr 

Cardot H, Ferraty F, Sarda P. 1999. \emph{Functional linear model}. Statistics and Probability Letters 45: 11{-}22.\cr 

Hall P, Hosseini{-}Nasab M. 2006. \emph{On properties of functional principal components analysis}. Journal of the Royal Statistical
Society B 68: 109{-}126.
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:   \code{\link{fregre.pc}}, \code{\link{summary.fregre.fd}} and \code{\link{predict.fregre.fd}}.\cr
Alternative method:  \code{\link{fregre.basis.cv}} and   \code{\link{fregre.np.cv}}.
}
\examples{
data(tecator)
absorp=tecator$absorp.fdata
ind=1:129
x=absorp[ind,]
y=tecator$y$Fat[ind]
out=fregre.pc.cv(x,y,7,criteria="SICc")
out$pc.opt
summary(out$fregre.pc)
pc<-fdata2pc(x,8)
out2=fregre.pc.cv(pc,y,criteria="SICc")
out2$pc.opt
summary(out2$fregre.pc)
out3=fregre.pc.cv(pc,y,criteria="SICc",rn=TRUE)
summary(out3$fregre.pc)

}
\keyword{regression}