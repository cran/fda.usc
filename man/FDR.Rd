\name{FDR}
\Rdversion{1.1}
\alias{FDR}
\alias{pvalue.FDR}
\title{	False Discorvery Rate (FDR)}
\description{ 	Compute the False Discovery Rate for a vector of p-values and alpha	value.}
\usage{
FDR(pvalues,alpha=0.95,dep=1)
pvalue.FDR(pvalues,dep=1)
}
\arguments{
  \item{pvalues}{ Vector of p-values}
	\item{alpha}{  Alpha value (level of significance).}
  \item{dep}{ Parameter dependence test. By default \code{dep = 1}, direct dependence between tests.}
}
\value{
Return:\cr
	\item{out.FDR}{ \code{=TRUE}. If  there are significative differences.}
  \item{pv.FDR}{  p-value for False Discovery Rate test.}
  }
\details{
\code{FDR} method is used for multiple hypothesis testing to correct problems of multiple contrasts.\cr
If \code{dep = 1}, the tests are positively correlated, for example when many tests are the same contrast. \cr
If  \code{dep < 1} the tests are negatively correlated.
}  
\references{
Benjamini, Y., Yekutieli, D. (2001). \emph{The control of the false discovery rate in multiple testing under dependency}. Annals of Statistics. 29 (4): 1165{-}1188. DOI{:}10.1214/aos/1013699998.
}
\author{
Febrero-Bande, M.  and Oviedo de la Fuente, M.
}
\seealso{
Function used in \code{\link{anova.RPm}} 
}
\examples{
 p=seq(1:50)/1000
 FDR(p)
 pvalue.FDR(p)
 FDR(p,alpha=0.9999)
 FDR(p,alpha=0.9)
 FDR(p,alpha=0.9,dep=-1)
}