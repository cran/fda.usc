\name{anova.RPm}
\Rdversion{1.1}
\alias{anova.RPm}
\alias{anova.RPm.boot}
\alias{summary.anova}
%\alias{intercambio.l}
%\alias{intercambio}
\title{ Functional ANOVA with Random Project.}
\description{
  The procedure is based on the analysis of randomly chosen one-dimensional projections. The function tests ANOVA models for functional data with continuous covariates and perform special contrasts for the factors in the formula.
 	}
\usage{
%anova.RPm(object,formula,data.fac,RP=min(30,ncol(object)),
%alpha=0.95,hetero=TRUE,pr=FALSE,w=rep(1,ncol(object)),
%nboot=0,contrast=NULL,...)
\method{anova}{RPm}(object,formula,data.fac,RP=min(30,ncol(object)),
alpha=0.95,hetero=TRUE,pr=FALSE,w=rep(1,ncol(object)),
nboot=0,contrast=NULL,\dots)

}
\arguments{
  \item{object}{ Functional response data.  Object with class fdata with \code{n} curves discretizated in \code{m} points.}
	\item{formula}{ as \link[stats]{formula}.}
	\item{data.fac}{ Explanatory variables.  Data frame  with dimension (\code{n} x \code{p}), where \code{p} are the number of factors or covariate considered.}
  \item{RP}{ Vector of random projections. }	
	\item{alpha}{ Alpha value, by defalult \code{alpha}=0.95.}
  \item{hetero}{ =TRUE (by default) heteroskedastic ANOVA.}
  \item{pr}{ If TRUE, print intermediate results. } 
  \item{w}{ Vector of weights.} 
  \item{nboot}{ Number of bootstrap samples, by defult no bootstrap sample, \code{nboot}=0. } 
  \item{contrast}{ List of special contrast  to be used ; by default no special contrasts are used (\code{contrast}=\code{NULL}).} 
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{
The function allows user-defined contrasts.  The list of contrast  to be used for some of the factors in the formula.  Each  contrast matrix in the list has  \code{r} rows, where \code{r} is the number of factor levels. The user can also request special predetermined contrasts, for example using the \code{\link{contr.helmert}}, \code{\link{contr.sum}} or \code{\link{contr.treatment}} functions.\cr


The function returns (by default) the significance of the variables using the Bonferroni test and the False Discovery Rate test.
Bootstrap procedure provides more precision 
}
\value{
Shows:\cr
\item{proj}{ The projection value of each point on the curves. Matrix with dimension (\code{RP} x \code{m}), where \code{RP} is the number of projection and \code{m} are the points observed in each projection curve.}
\item{mins}{ minimum number for each random projection.}
\item{result}{ p-value for each random projection.}
\item{test.Bonf}{ significance (TRUE or FALSE)  for vector of random projections \code{RP} in columns and factor (and special contrast) by rows.}
\item{p.Bonf}{ p-value   for vector of random projections \code{RP} in columns and factor (and special contrast) by rows.}
\item{test.fdr}{ False Discovery Rate (TRUE or FALSE)   for vector of random projections \code{RP} in columns and factor (and special contrast) by rows.}
\item{p.fdr}{  p-value of False Discovery Rate  for vector of random projections \code{RP} in columns and factor (and special contrast) by rows.}
\item{test.Boot}{ False Discovery Rate (TRUE or FALSE)   for vector of random projections \code{RP} in columns and factor (and special contrast) by rows.}
\item{p.Boot}{  p-value of Bootstrap sambple for vector of random projections \code{RP} in columns and factor (and special contrast) by rows.}
}
\note{If \code{hetero=TRUE} then all factors must be categorical.}
\references{
Cuesta-Albertos, J.A., Febrero-Bande, M. \emph{A simple multiway ANOVA for functional data.}
Test 2010, DOI \bold{10.1007/s11749-010-0185-3}.
}
\author{
Cuesta-Albertos, J.A., Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as: \code{\link{anova.hetero}} 
}
\examples{

# ex anova.hetero
data(phoneme)
names(phoneme)
data=as.data.frame(phoneme$learn[["data"]])
group=phoneme$classlearn
n=nrow(data)
group.rand=as.factor(sample(rep(1:3,len=n),n))
RP=c(2,5,15,30)

#ex 1: real factor and random factor
m03=data.frame(group,group.rand)
resul1=anova.RPm(data,~group+group.rand,m03,RP=c(5,30))
summary.anova(resul1)

#ex 2: real factor with special contrast
m0=data.frame(group)
cr5=contr.sum(5)   #each level vs last level
resul03c1=anova.RPm(data,~group,m0,contrast=list(group=cr5))
summary.anova(resul03c1)

#ex 3: random factor with special contrast
m0=data.frame(group.rand)
cr3=contr.sum(3)   #each level vs last level
resul03c1=anova.RPm(data,~group.rand,m0,contrast=list(group.rand=cr3))
summary.anova(resul03c1)
}
\keyword{anova}