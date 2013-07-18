\name{flm.Ftest}
\alias{flm.Ftest}
\alias{Ftest.statistic}

\title{F-test for the Functional Linear Model with scalar response}
\description{
   The function \code{flm.Ftest} tests the null hypothesis of no interaction between a functional covariate and a scalar response inside the Functional Linear Model (FLM): \eqn{Y=\big<X,\beta\big>+\epsilon}{Y=<X,\beta>+\epsilon}. The null hypothesis is \eqn{H_0:\,\beta=0}{H_0: \beta=0} and the alternative is \eqn{H_1:\,\beta\neq 0}{H_1: \beta\neq 0}. 
   The way of testing the null hypothesis is via a functional extension of the classical F-test (see Details). 
 }

\usage{
Ftest.statistic (X.fdata, Y)
flm.Ftest (X.fdata, Y, B=5000, verbose=TRUE)
}

\arguments{
  \item{X.fdata}{Functional covariate for the FLM. The object must be in the class \code{\link{fdata}}.}
  \item{Y}{Scalar response for the FLM. Must be a vector with the same number of elements as functions are in \code{X.fdata}.}
  \item{B}{Number of bootstrap replicates to calibrate the distribution of the test statistic. \code{B=5000} replicates are the recommended for carry out the test, although for exploratory analysis (\bold{not inferential}), an acceptable less time-consuming option is \code{B=500}.}
  \item{verbose}{Either to show or not information about computing progress.}
 }

\details{
 The Functional Linear Model with scalar response (FLM), is defined as \eqn{Y=\big<X,\beta\big>+\epsilon}{Y=<X,\beta>+\epsilon},
 for a functional process \eqn{X}{X} such that \eqn{E[X(t)]=0}{E[X(t)]=0}, \eqn{E[X(t)\epsilon]=0}{E[X(t)\epsilon]=0} for all \eqn{t}{t} and for a scalar variable \eqn{Y}{Y} such that \eqn{E[Y]=0}{E[Y]=0}.
 The \emph{functional F-test} is defined as \deqn{T_n=\bigg\|\frac{1}{n}\sum_{i=1}^n (X_i-\bar X)(Y_i-\bar Y)\bigg\|,}{||\frac{1}{n} \sum_{i=1}^n (X_i-\bar X)(Y_i-\bar Y)||,} where \eqn{\bar X}{\bar X} is the functional mean of \eqn{X}{X}, \eqn{\bar Y}{\bar Y} is the ordinary mean of \eqn{Y}{Y} and \eqn{\|\cdot\|}{||.||} is the \eqn{L^2}{L^2} functional norm.
 The statistic is computed with the function \code{Ftest.statistic}. The distribution of the test statistic is approximated by a wild bootstrap on the residuals, using the \emph{golden section bootstrap}.

 For further details of the test and its implementation see Gonzalez-Manteiga \emph{et al.} (2012) and Garcia-Portugues \emph{et al.} (2012), respectively.

}

\value{
  The value for \code{Ftest.statistic} is simply the F-test statistic. The value for \code{flm.Ftest} is an object with class \code{"htest"} whose underlying structure is a list containing the following components:
  \item{statistic}{The value of the F-test statistic.}
  \item{boot.statistics}{A vector of length \code{B} with the values of the bootstrap F-test statistics.}
  \item{p.value}{The p-value of the test.}
  \item{method}{The character string "Functional Linear Model F-test".}
  \item{B}{The number of bootstrap replicates used.}
  \item{data.name}{The character string "Y=<X,0>+e"}
}

\references{

Garcia-Portugues, E., Gonzalez-Manteiga, W. and Febrero-Bande, M. (2012). A goodness--of--fit test for the functional linear model with scalar
response. \url{http://arxiv.org/abs/1205.6167}

Gonzalez-Manteiga, W., Gonzalez-Rodriguez, G., Martinez-Calvo, A. and Garcia-Portugues, E. Bootstrap independence test for functional linear models. 

}

\note{No NA's are allowed neither in the functional covariate nor in the scalar response.}

\author{Eduardo Garcia-Portugues. Please, report bugs and suggestions to \if{latex}{\cr}\email{eduardo.garcia@usc.es}}

\seealso{\code{\link{rwild}}, \code{\link{flm.test}}, \code{\link{dfv.test}}}

\examples{

## Simulated example ##

X=rproc2fdata(n=50,t=seq(0,1,l=101),sigma="OU")

beta0=fdata(mdata=rep(0,length=101)+rnorm(101,sd=0.05),
argvals=seq(0,1,l=101),rangeval=c(0,1))
beta1=fdata(mdata=cos(2*pi*seq(0,1,l=101))-(seq(0,1,l=101)-0.5)^2+
rnorm(101,sd=0.05),argvals=seq(0,1,l=101),rangeval=c(0,1))

# Null hypothesis holds
Y0=drop(inprod.fdata(X,beta0)+rnorm(50,sd=0.1))

# Null hypothesis does not hold
Y1=drop(inprod.fdata(X,beta1)+rnorm(50,sd=0.1))

# Do not reject H0
#flm.Ftest(X,Y0,B=100)
# flm.Ftest(X,Y0,B=5000)

# Reject H0
#flm.Ftest(X,Y1,B=100)
# flm.Ftest(X,Y1,B=5000)

}

\keyword{htest}
\keyword{models}
\keyword{regression}
