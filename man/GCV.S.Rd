\name{GCV.S}
\Rdversion{1.1}
\alias{GCV.S}
\title{ The generalized cross-validation (GCV) score.}
\description{
 The generalized cross-validation (GCV) score. 
}
\usage{
GCV.S(y,S,criteria="GCV",W=diag(1,ncol=ncol(S),nrow=nrow(S)),
trim=0,draw=FALSE,metric=metric.lp,...)
}

\arguments{
  \item{y}{ Matrix  of  set cases with dimension (\code{n} x \code{m}), where \code{n} is the number of curves and \code{m} are the points observed in each curve.}
  \item{S}{ Smoothing matrix, see \code{\link{S.NW}}, \code{\link{S.LLR}} or \eqn{S.KNN}. }
  \item{criteria}{ The penalizing function. By default \emph{"Rice"} criteria. Possible values are \emph{"GCV"}, \emph{"AIC"}, \emph{"FPE"}, \emph{"Shibata"}, \emph{"Rice"}. }
  \item{W }{ Matrix of weights.}    
	\item{trim}{ The alpha of the trimming.}
  \item{draw}{ =TRUE, draw the curves, the sample median and trimmed mean.}	  
  \item{metric}{ Metric function, by default \code{\link{metric.lp}}.}
  \item{\dots}{ Further arguments passed to or from other methods.}  
}

\details{  
 \deqn{GCV(h)=p(h) \Xi(n^{-1}h^{-1})}{GCV(h)=p(h) \Xi(1/(n h))} \cr

 Where  \cr
A.-If \code{trim=0}:\cr 
\deqn{p(h)=\frac{1}{n} \sum_{i=1}^{n}{\Big(y_i-r_{i}(x_i)\Big)^{2}w(x_i)}}{p(h)=1/n \sum_{i=1:n} ( y_i - r_i(x_i) )^2 w(x_i)} \cr
B.-If \code{trim>0}:\cr
\deqn{p(h)=\frac{1}{l} \sum_{i=1}^{l}{\Big(y_i-r_{i}(x_i)\Big)^{2}w(x_i)}}{p(h)=1/l \sum_{i=1:l} ( y_i - r_i(x_i) )^2 w(x_i)} \cr
\code{l}: index of \code{(1-trim)} curves with less error.
where \eqn{h} is the bandwidth parameter, \eqn{w} the weights and the penalty function \eqn{\Xi} can be selected from the following criteria: \cr 
\itemize{
\item{Generalized Cross-validation (GCV):\cr
	  \deqn{\Xi_{GCV}(n^{-1}h^{-1})=(1-n^{-1}S_{ii})^{-2}}{\Xi_{GCV}(1/(n h))=(1-1/n S_{ii})^{-2} }}
\item{Akaike's Information Criterion (AIC): \cr
     \deqn{\Xi_{AIC}(n^{-1}h^{-1})=exp(2n^{-1}S_{ii})}{\Xi_{AIC}(1/(n h))=exp(2n^{-1} S_{ii})}}
\item{Finite Prediction Error (FPE): \cr
     \deqn{\Xi_{FPE}(n^{-1}h^{-1})=\frac{(1+n^{-1}S_{ii})}{(1-n^{-1}S_{ii})}}{\Xi_{FPE}(1/(n h))= (1+1/n S_{ii}) / (1-1/n S_{ii})} }
\item{Shibata's model selector (Shibata): \cr
     \deqn{\Xi_{Shibata}(n^{-1}h^{-1})=(1+2n^{-1}S_{ii})}{\Xi_{Shibata}(1/(n h))=(1+2n^{-1} S_{ii})}}
\item{Rice's bandwidth selector (Rice):\cr
     \deqn{\Xi_{Rice}(n^{-1}h^{-1})=(1-2n^{-1}S_{ii})^{-1}}{\Xi_{Rice}(1/(n h))=(1-2 n^{-1}S_{ii})^{-1}}}
}
where \eqn{S_{ii}} the \eqn{i}th diagonal element of the smoothing matrix \eqn{S}, in see \code{\link{S.NW}}, \code{\link{S.LLR}} or \eqn{S.KNN}. 
}     
\value{
    \item{res}{ Returns GCV score calculated for input parameters.  } 
}
\references{
Wasserman, L. \emph{All of Nonparametric Statistics}. Springer Texts in Statistics, 2006. \cr

Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University Press, 1994.
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as  \code{\link{min.np}}. \cr
Alternative method:  \code{\link{CV.S}}
}
\examples{

data(phoneme)
mlearn<-phoneme$learn
tt<-1:ncol(mlearn)
S1 <- S.NW(tt,2.5)
S2 <- S.LLR(tt,2.5)
S3 <- S.KNN(tt,2.5)
gcv1 <- GCV.S(mlearn, S1)
gcv2 <- GCV.S(mlearn, S2)
gcv3 <- GCV.S(mlearn, S3)
gcv4 <- GCV.S(mlearn, S1,criteria="AIC")
gcv5 <- GCV.S(mlearn, S1,trim=0.01,draw=TRUE)
gcv1; gcv2; gcv3; gcv4; gcv5
 }

\keyword{utilities} 