\name{anova.onefactor}
\Rdversion{1.1}
\alias{anova.onefactor}
\title{One--way anova model for functional data}
\description{
%The function tests the null hypothesis for functional data the equality of respective mean functions. 
%The function tests ANOVA models for functional data  perform contrasts for the factors in the formula.

%The function \code{anova.Ftest} tests the  null hypothesis of a equality of mean functions of functional data based on the an asymptotic version of the anova F--test,
%  \deqn{H_0:\, m_1=\ldots=m_k}{H_0: m_1=,\ldots=m_k,}
%  versus an alternative 
%    \deqn{H_a: m_i\neq=m_j,\ for\ any\ i,j,}{H_a: m_i!=m_j, for any i,j,}
             
One--way anova model for k independent samples of functional data. The function contrasts the null hypothesis of equality of mean functions of functional data based on the an asymptotic version of the anova F--test. 
 \deqn{H_0:\, m_1=\ldots=m_k}{H_0: m_1=,\ldots=m_k,}
%Given k independent samples of functional data the function contrasts the null hypothesis of equality of their respective mean functions. The procedure is based on the an asymptotic version of the well-known anova F--test. 
% \deqn{H_0:\, m_1=\ldots=m_k}{H_0: m_1=,\ldots=m_k,}
%the null hypothesis of equality of the average curves in both groups.
   }
\usage{
%anovaCSDA=function(fdata,group,nboot=100,plot=TRUE,verbose=TRUE,...)
\method{anova}{onefactor}(object,group,nboot=100,plot=FALSE,verbose=FALSE,\dots)
}
\arguments{
  \item{object}{ functional response data.  fdata class object with \code{n} curves.}
	\item{group}{ a factor specifying the class for each curve.}
  \item{nboot}{ number of bootstrap samples. }   
  \item{plot}{ if \code{TRUE}, plot the mean of each factor level and the results of test. } 
  \item{verbose}{ if \code{TRUE}, print intermediate results. }   
  \item{\dots}{ further arguments passed to or from other methods.}
}


\details{

The function returns the p--value of  test using one--way anova model over \code{nboot} runs.
}
\value{
returns:\cr
\item{pvalue}{ p-value.}          %   probability of rejecting the null hypothesis H0 at a significance level
\item{stat}{ statistic value of test.}
\item{wm}{ statistic values of bootstrap resamples.}
}
\references{
Cuevas, A., Febrero, M., & Fraiman, R. (2004). \emph{An anova test for functional data.} Computational statistics & data analysis, \bold{47}(1), 111-122.
}
\author{
 Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{
See Also as: \code{\link{anova.RPm}} 
}
\examples{
\dontrun{
data(MCO)
grupo<-MCO$classintact
datos<-MCO$intact
res=anova.onefactor(datos,grupo,nboot=50,plot=TRUE)

grupo<-MCO$classpermea
datos<-MCO$permea
res=anova.onefactor(datos,grupo,nboot=50,plot=TRUE)
}  
}
\keyword{anova}