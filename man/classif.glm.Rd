% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classif.glm.R
\name{classif.glm}
\alias{classif.glm}
\title{Classification Fitting Functional Generalized Linear Models}
\usage{
classif.glm(
  formula,
  data,
  family = binomial(),
  weights = "equal",
  basis.x = NULL,
  basis.b = NULL,
  type = "1vsall",
  prob = 0.5,
  CV = FALSE,
  ...
)
}
\arguments{
\item{formula}{an object of class \code{formula} (or one that can be coerced
to that class): a symbolic description of the model to be fitted. The
details of model specification are given under \code{Details}.}

\item{data}{List that containing the variables in the model.}

\item{family}{a description of the error distribution and link function to
be used in the model. This can be a character string naming a family
function, a family function or the result of a call to a family function.
(See \code{\link{family}} for details of family functions).}

\item{weights}{Weights:
\itemize{
\item if \code{character} string \code{='equal'} same weights for each observation (by default) and
\code{='inverse'} for inverse-probability of weighting.   
\item if \code{numeric} vector of length \code{n}, Weight values of each observation.
}}

\item{basis.x}{List of basis for functional explanatory data estimation.}

\item{basis.b}{List of basis for functional beta parameter estimation.}

\item{type}{If type is\code{"1vsall"}  (by default) 
a maximum probability scheme is applied: requires G binary classifiers.
If type is \code{"majority"}  (only for multicalss classification G > 2) 
a voting scheme is applied: requires  G (G - 1) / 2 binary classifiers.}

\item{prob}{probability value used for binari discriminant.}

\item{CV}{=TRUE, Cross-validation (CV) is done.}

\item{\dots}{Further arguments passed to or from other methods.}
}
\value{
Return \code{glm} object plus:
\itemize{
\item \code{formula}{ formula.}
\item \code{data}{ List that containing the variables in the model.} 
\item \code{group}{ Factor of length \emph{n}} 
\item \code{group.est}{ Estimated vector groups}
\item \code{prob.classification}{ Probability of correct classification by group.}
\item \code{prob.group}{ Matrix of predicted class probabilities. For each
functional point shows the probability of each possible group membership.}
\item \code{max.prob}{ Highest probability of correct classification.}
}
}
\description{
Computes functional classification using functional (and non functional)
explanatory variables by basis representation.

The first item in the \code{data} list is called \emph{"df"} and is a data
frame with the response and non functional explanatory variables, as
\code{\link{glm}}.\cr

Functional covariates of class \code{fdata} or \code{fd} are introduced in
the following items in the \code{data} list.\cr \code{basis.x} is a list of
basis for represent each functional covariate. The basis object can be
created by the function: \code{\link{create.pc.basis}}, \code{\link{pca.fd}}
\code{\link{create.pc.basis}}, \code{\link{create.fdata.basis}} o
\code{\link{create.basis}}.\cr \code{basis.b} is a list of basis for
represent each functional beta parameter. If \code{basis.x} is a list of
functional principal components basis (see \code{\link{create.pc.basis}} or
\code{\link{pca.fd}}) the argument \code{basis.b} is ignored.
}
\note{
If the formula only contains a non functional explanatory variables
(multivariate covariates), the function compute a standard \code{\link{glm}}
procedure.
}
\examples{
\dontrun{
data(phoneme)
ldat <- ldata("df" = data.frame(y = phoneme[["classlearn"]]),
             "x" = phoneme[["learn"]])
a1 <- classif.glm(y ~ x, data = ldat)
summary(a1)
newldat <- ldata("df" = data.frame(y = phoneme[["classtest"]]),
                "x" = phoneme[["test"]])
p1 <- predict(a1,newldat)
table(newldat$df$y,p1)
sum(p1==newldat$df$y)/250
}
}
\references{
Ramsay, James O., and Silverman, Bernard W. (2006), \emph{
Functional Data Analysis}, 2nd ed., Springer, New York.

McCullagh and Nelder (1989), \emph{Generalized Linear Models} 2nd ed.
Chapman and Hall.

Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics
with S}, New York: Springer.  %Wood (2001) mgcv:GAMs and Generalized Ridge
Regression for R. R News 1(2):20-25
}
\seealso{
See Also as: \code{\link{fregre.glm}}.\cr %Alternative method:
\code{\link{classif.gsam}} and \code{\link{classif.gkam}}.
}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente
\email{manuel.oviedo@udc.es}
}
\keyword{classif}
