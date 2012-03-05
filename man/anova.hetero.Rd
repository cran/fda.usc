\name{anova.hetero}
\Rdversion{1.1}
\alias{anova.hetero}
\title{ ANOVA for heteroscedastic data}
\description{
  Univariate ANOVA for heteroscedastic data.
}
\usage{
%anova.hetero(object=NULL,formula,pr=FALSE,contrast=NULL,...)
\method{anova}{hetero}(object=NULL,formula,pr=FALSE,contrast=NULL,...)
}
\arguments{
  \item{object}{ A data frame with  dimension (\code{n} x \code{p+1}).  In the first column contains the \code{n} response values and on the following \code{p} columns the explanatory variables specified in the formula.}
  \item{formula}{ as \link[stats]{formula}.}
  \item{pr}{ If TRUE, print intermediate results. } 
  \item{contrast}{ List of special contrast  to be used, by default no special contrasts are used (\code{contrast}=\code{NULL}).}
  \item{\dots}{ Further arguments passed to or from other methods.} 
}
\details{ 
This function fits a univariate analysis of variance model and  allows calculate special contrasts defined by the user.  The list of special contrast  to be used for some of the factors in the formula.  Each  matrix of the list has  \code{r} rows and \code{r-1} columns.\cr

The user can also request special predetermined contrasts, for example using \code{\link{contr.helmert}}, \code{\link{contr.sum}} or \code{\link{contr.treatment}} functions.\cr
}

\value{
Return: \cr
\item{ans}{ A list with components including: the  Beta estimation \code{Est}, the factor degrees of freedom \code{df1}, the residual degrees of freedom \code{df2} and \code{p-value} for each factor.}
\item{contrast}{ List of  special contrasts. \cr}
}
\note{It only works with categorical variables.}
\references{

Brunner, E., Dette, H., Munk, A. \emph{Box-Type Approximations in Nonparametric Factorial Designs.}
Journal of the American Statistical Association, Vol. 92, No. 440 (Dec., 1997), pp. 1494{-}1502.
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:  \code{\link{anova.RPm}} 
}
\examples{

data(phoneme)
ind=1 # beetwen 1:150
fdataobj=data.frame(phoneme$learn[["data"]][,ind])
n=dim(fdataobj)[1]
group<-factor(phoneme$classlearn)

#ex 1: real factor and random factor
group.rand=as.factor(sample(rep(1:3,n),n))
f=data.frame(group,group.rand)
mm=data.frame(fdataobj,f)
colnames(mm)=c("value","group","group.rand")
out1=anova.hetero(object=mm[,-2],value~group.rand,pr=FALSE)
out2=anova.hetero(object=mm[,-3],value~group,pr=FALSE)
out1
out2

#ex 2: real factor, random factor and  special contrasts
cr5=contr.sum(5)  #each level vs last level
cr3=c(1,0,-1)			#first level vs last level
out.contrast=anova.hetero(object=mm[,-3],value~group,pr=FALSE,
contrast=list(group=cr5))
out.contrast
    }
\keyword{anova}
