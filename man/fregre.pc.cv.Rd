\name{fregre.pc.cv}
\Rdversion{1.1}
\alias{fregre.pc.cv}
\alias{fregre.combn}
\title{Functional Regression using selection of number of principal components}
\description{																						
Functional Regression with scalar response using selection of number of principal components through cross-validation. The algoritm selects the  principal components with best estimates the response. The	select is performed by cross-validation (CV) or  Model Selecttion Criteria (MSC).  After is computing  functional regression using the best selection of Principal Components.
	}
\usage{
fregre.pc.cv(fdataobj,y,kmax=8,criteria="SIC",...)
}
\arguments{
 \item{fdataobj}{ \code{\link{fdata}} class object.}  
	\item{y}{ Scalar response with length \code{n}.}
	\item{kmax}{ Maximum number of principal comoponents.}
	\item{criteria}{ Type of cross-validation (CV) or  Model Selecttion Criteria (MSC) applied. Possible values are \emph{"CV"}, \emph{"AIC"}, \emph{"AICc"}, \emph{"SIC"}.} 	
 \item{\dots}{ Further arguments passed to \code{\link{fregre.pc}}}  
}
\details{
The algorithm for selecting the best principal components from the first \code{kmax} is:

 1.-Calculate the best principal component (\emph{pc.order[1]}) between \code{kmax} first.\cr
 2.-Calculate the second-best principal component (\emph{pc.order [2]}) between the \code{(kmax-1)} and calculate the criteria value of the two principal  components.\cr
 3.-The process (point 1 and 2) is repeated until \code{kmax} principal component (\emph{pc.order[kmax]}).\cr
 The method select the  principal components (\code{pc.opt}=\code{pc.order[1:k.min]}) with minimum MSC criteria.
 
 4.-Finally, is computing  functional regression between functional explanatory variable and scalar response using the best selection of Principal Components, for more details see \code{\link{fregre.pc}}.
}
\value{
Return:\cr
\item{pc.opt}{ Index of principal comoponents selected.}
\item{MSC.min}{ Minimum Model Selecttion Criteria (MSC) value for \code{pc.opt} principal components.}
\item{pc.order}{ Complete sequence of principal components selected in each iteration.}
\item{MSC}{ Complete sequence of minimum  Model Selecttion Criteria (MSC) value for \code{pc.opt} principal components in each iteration.}
\item{fregre.pc}{ object for \code{pc.opt} principal components.}
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
Alternative method: \code{\link{fregre.basis.cv}} and   \code{\link{fregre.np.cv}}.
}
\examples{
data(tecator)
absorp=tecator$absorp.fdata
ind=1:129
x=absorp[ind,]
y=tecator$y$Fat[ind]
out=fregre.pc.cv(x,y,7,"SICc")
out$pc.opt
summary(out$fregre.pc)
}
\keyword{regression}