\name{CV.S}
\Rdversion{1.1}
\alias{CV.S}

\title{ The cross-validation (CV) score}
\description{
The cross-validation (CV)  score.
}
\usage{
CV.S(y,S,W=diag(ncol(S)),trim=0,draw=FALSE,metric=metric.lp,...)
}

\arguments{
  \item{y}{ Matrix  of  set cases with dimension (\code{n} x \code{m}), where \code{n} is the number of curves and \code{m} are the points observed in each curve.}
  \item{S}{ Smoothing matrix, see \code{\link{S.NW}}, \code{\link{S.LLR}} or \eqn{S.KNN}.}
  \item{W }{ Matrix of weights.}    
	\item{trim}{ The alpha of the trimming.}
  \item{draw}{ =TRUE, draw the curves, the sample median and trimmed mean.}	
  \item{metric}{ Metric function, by default \code{\link{metric.lp}}.}
  \item{\dots}{ Further arguments passed to or from other methods.}  
  }

\details{
Compute the leave-one-out cross-validation score.\cr
A.-If \code{trim=0}:\cr
\deqn{CV(h)=\frac{1}{n} \sum_{i=1}^{n}{\Bigg(\frac{y_i-r_{i}(x_i)}{(1-S_{ii})}\Bigg)^{2}w(x_{i})}}{CV(h)=1/n\, \sum_i ((y_i\, -\, r_{i}(x_i))\, /\, (1\, -\, S_ii))^2\,  w(x_i),\, i=1,...,n} 
\eqn{S_{ii}}{Sii} is the ith diagonal element of the smoothing matrix \eqn{S}{S}.

B.-If \code{trim>0}:\cr
\deqn{CV(h)=\frac{1}{l} \sum_{i=1}^{l}{\Bigg(\frac{y_i-r_{i}(x_i)}{(1-S_{ii})}\Bigg)^{2}w(x_{i})}}{CV(h)=1/n\ \sum_i ((y_i-r_{i}(x_i))/(1-S_ii))^2 w(x_i),\, i=1,...,l}
\eqn{S_{ii}}{Sii} is the ith diagonal element of the smoothing matrix \eqn{S}{S} and l the index of  \code{(1-trim)} curves with less error.
}
\value{
    \item{res}{ Returns CV score calculated for input parameters.  } 
}
\references{
Wasserman, L. \emph{All of Nonparametric Statistics}. Springer Texts in Statistics, 2006.
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as  \code{\link{min.np}} \cr
Alternative method:  \code{\link{GCV.S}}
}
\examples{

data(tecator)
x<-tecator$absorp.fdata
np<-ncol(x)
tt<-1:np
 S1 <- S.NW(tt,3,Ker.epa)
 S2 <- S.LLR(tt,3,Ker.epa)
 S3 <- S.KNN(tt,3,Ker.epa)
 S4 <- S.NW(tt,5,Ker.epa)
 S5 <- S.LLR(tt,5,Ker.epa)
 S6 <- S.KNN(tt,5,Ker.epa)
 cv1 <- CV.S(x, S1)
 cv2 <- CV.S(x, S2)
 cv3 <- CV.S(x, S3)
 cv4 <- CV.S(x, S4)
 cv5 <- CV.S(x, S5)
 cv6 <- CV.S(x, S6)
 cv7 <- CV.S(x, S4,trim=0.1,draw=TRUE)
 cv1;cv2;cv3;cv4;cv5;cv6;cv7
 }

\keyword{utilities} 