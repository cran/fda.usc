\name{metric.hausdorff}
\Rdversion{1.1}
\alias{metric.hausdorff}
\title{
Compute the Hausdorff distances between two curves.
}
\description{	
Hausdorff distance is the greatest of all the distances from a point in one curve to the closest point in the other curve (been closest the euclidean distance).
}
\usage{
metric.hausdorff(fdata1, fdata2 = fdata1) 
}
\arguments{  
  \item{fdata1}{Curves 1 of \code{fdata} class. The dimension of \code{fdata1} object is (\code{n1} x \code{m}), where \code{n1} is the number of points observed in \code{t} coordinates with lenght \code{m}.}
  \item{fdata2}{Curves 2 of \code{fdata} class. The dimension of \code{fdata2} object is (\code{n2} x \code{m}), where \code{n2} is the number of points observed in \code{t} coordinates with lenght \code{m}.}
}
\details{
Let \eqn{G(X)=\left\{ (t,X(t))\in R^2 \right\}}{G(X)={(t,X(t)) \in R^2}} and 
\eqn{G(Y)=\left\{(t,Y(t))\in R^2\right\}}{G(Y)={(s,Y(s)) \in R^2}}  be two graphs of the considered curves \eqn{X} and \eqn{Y} respectively, the Hausdorff distance \eqn{d_H(X, Y)} is defined as,

\deqn{ d_H(X,Y)=max\left\{ sup_{x\in G(X)} inf_{y\in G(Y)} d_2(x,y), sup_{y\in G(Y)} inf_{x\in G(X)}d_2(x,y)\right\},}{d_H(X,Y)=max{sup_{x\in G(X)} inf_{y\in G(Y)}d_2(x,y),sup_{y\in G(Y)} inf_{x\in G(X)}d_2(x,y)},}
where \eqn{d_2(x,y)} is the euclidean distance, see \code{\link{metric.lp}.}
}
%\references{
%http://en.wikipedia.org/wiki/Hausdorff_distance
%Huttenlocher, D. P., Klanderman, G. A., & Rucklidge, W. J. (1993). Comparing images using the Hausdorff distance. Pattern Analysis and Machine Intelligence, IEEE Transactions on, 15(9), 850-863.
%Genovese, C. R., Perone-Pacifico, M., Verdinelli, I. and Wasserman, L. (2012). Manifold estimation and singular deconvolution under hausdorff loss. \emph{The Annals of Statistics} 40, 941-963.
%}
\author{Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}}
\examples{

\dontrun{   
data(poblenou)
nox<-poblenou$nox[1:6]
# Hausdorff vs maximum distance
out1<-metric.hausdorff(nox)       
out2<-metric.lp(nox,lp=0) 
out1
out2
par(mfrow=c(1,3))
plot(nox)
plot(hclust(as.dist(out1)))
plot(hclust(as.dist(out2)))
}   
}
\keyword{ cluster }