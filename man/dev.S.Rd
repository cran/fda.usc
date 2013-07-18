\name{dev.S}
\Rdversion{1.1}
\alias{dev.S}
\title{ The deviance score .}
\description{
 Returns the deviance of a fitted model object by GCV score.
}
\usage{
dev.S(y, S, obs,family = gaussian(),off,offdf,criteria="GCV",
W = diag(1, ncol = ncol(S), nrow = nrow(S)), trim = 0, draw = FALSE,...)
}

\arguments{
  \item{y}{ Matrix  of  set cases with dimension (\code{n} x \code{m}), where \code{n} is the number of curves and \code{m} are the points observed in each curve.}
    \item{obs}{ observed response.}
  \item{S}{ Smoothing matrix.}
  \item{family}{ a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.) }
    \item{off}{ off}
    \item{offdf}{ off, degrees of freedom}
  \item{criteria}{ The penalizing function. By default \emph{"Rice"} criteria. Possible values are \emph{"GCV"}, \emph{"AIC"}, \emph{"FPE"}, \emph{"Shibata"}, \emph{"Rice"}. }
  \item{W }{ Matrix of weights.}    
	\item{trim}{ The alpha of the trimming.}
  \item{draw}{ =TRUE, draw the curves, the sample median and trimmed mean.}	  
  \item{\dots}{ Further arguments passed to or from other methods.}  
}

\details{  
up to a constant, minus twice the maximized log-likelihood. Where sensible, the constant is chosen so that a saturated model has deviance zero.\cr

 \deqn{GCV(h)=p(h) \Xi(n^{-1}h^{-1})}{} \cr

 Where  \cr
 \deqn{p(h)=\frac{1}{n} \sum_{i=1}^{n}{\Big(y_i-r_{i}(x_i)\Big)^{2}w(x_i)}}{} \cr
 and penalty function  \deqn{\Xi()}{} can be selected from the following criteria: \cr \cr
Generalized Cross-validation (GCV):\cr
	  \deqn{\Xi_{GCV}(n^{-1}h^{-1})=(1-n^{-1}S_{ii})^{-2}}{} \cr
Akaike's Information Criterion (AIC): \cr
     \deqn{\Xi_{AIC}(n^{-1}h^{-1})=exp(2n^{-1}S_{ii})}{} \cr
Finite Prediction Error (FPE) \cr
     \deqn{\Xi_{FPE}(n^{-1}h^{-1})=\frac{(1+n^{-1}S_{ii})}{(1-n^{-1}S_{ii})}}{} \cr
Shibata's model selector (Shibata): \cr
     \deqn{\Xi_{Shibata}(n^{-1}h^{-1})=(1+2n^{-1}S_{ii})}{}\cr
Rice's bandwidth selector (Rice):\cr
     \deqn{\Xi_{Rice}(n^{-1}h^{-1})=(1-2n^{-1}S_{ii})^{-1}}{} \cr
}
     
\value{
    \item{res}{ Returns GCV score calculated for input parameters.  } 
}
\references{
Wasserman, L. \emph{All of Nonparametric Statistics}. Springer Texts in Statistics, 2006. 

Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University Press, 1994.

Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{
See Also as  \code{\link{GCV.S}}. \cr
Alternative method:  \code{\link{CV.S}}
}
\examples{

data(phoneme)
mlearn<-phoneme$learn
np<-ncol(mlearn)
tt<-mlearn[["argvals"]]
S1 <- S.NW(tt,2.5)
gcv1 <- dev.S(mlearn$data[1,],obs=(sample(150)), S1,off=rep(1,150),offdf=3)
gcv2 <- dev.S(mlearn$data[1,],obs=sort(sample(150)), S1,off=rep(1,150),offdf=3)
 }

\keyword{utilities} 