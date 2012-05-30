\name{rber.gold}
\alias{rber.gold}

\title{Gold section bootstrap sampling}

\description{
Sampling from a binomial variable with values \eqn{\frac{1-\sqrt{5}}{2}}{(1-\sqrt 5)/2 }, \eqn{\frac{1+\sqrt{5}}{2}}{(1+\sqrt 5)/2} and probabilities \eqn{\frac{5+\sqrt{5}}{10}}{(5+\sqrt 5)/10}, \eqn{\frac{5-\sqrt{5}}{10}}{(5-\sqrt 5)/10}, respectively.
}

\usage{
rber.gold (n)
}

\arguments{
  \item{n}{Number of observations.}
}

\details{
For the construction of wild bootstrap residuals, sampling from a random variable \eqn{V}{V} such that \eqn{E[V^2]=0}{E[V^2]=0} and \eqn{E[V]=0}{E[V]=0} is needed. 
A simple and suitable \eqn{V}{V} is obtained with a binomial variable of the form:
\deqn{P\big\{ V=\frac{1-\sqrt{5}}{2} \big\} = \frac{5+\sqrt{5}}{10} \, and \, P\big\{ V=\frac{1+\sqrt{5}}{2} \big\} = \frac{5-\sqrt{5}}{10},}{P\{ V=(1-\sqrt 5)/2 \} = (5+\sqrt 5)/10 and P\{ V=(1+\sqrt 5)/2 \} = (5-\sqrt 5)/10,}
which leads to the \emph{golden section bootstrap}. If \code{e} denotes a vector of \code{n} residuals, the wild bootstrap residuals would be computed as \code{e*rber.gold(n)}.
}

\value{
A sample of length \code{n} of the random variable \eqn{V}{V}.
}


\author{Eduardo Garcia-Portugues. Please, report bugs and suggestions to \email{eduardo.garcia@usc.es}}

\seealso{\code{\link{flm.test}}, \code{\link{rbinom}}}

\examples{

# Sampling
samp=rber.gold(100)
mean(samp)
sd(samp)
samp

# Construction of wild bootstrap residuals
e=rnorm(200)
e.boot=e*rber.gold(200)
summary(e)
summary(e.boot)

}

\keyword{distribution}

