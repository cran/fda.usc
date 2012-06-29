\name{classif.glm}
\Rdversion{1.1}
\alias{classif.glm}
\title{Classification Fitting Functional Generalized Linear Models}
\description{
	Computes functional classification using functional (and non functional) explanatory variables by basis representation.
}
\usage{
classif.glm(formula,data,family = binomial(),
basis.x=NULL,basis.b=NULL,CV=FALSE,...)
}
\arguments{
 \item{formula}{ an object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under \code{Details}.}
  \item{data}{ List that containing the variables in the model. }
 \item{family}{ a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions). }
	\item{basis.x}{ List of basis for functional explanatory data estimation.}
	\item{basis.b}{ List of basis for functional beta parameter estimation.}
	\item{CV}{ =TRUE, Cross-validation (CV) is done.}		
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{ 
The first item in the \code{data} list  is called \emph{"df"} and  is a data frame  with the response and non functional explanatory variables, as \code{\link{glm}}.\cr

Functional covariates of class  \code{fdata} or \code{fd} are introduced in the following items in the \code{data} list.\cr
\code{basis.x} is a list of basis  for represent each functional covariate. The basis object can be created by the function: \code{\link{create.pc.basis}}, \code{\link{pca.fd}} \code{\link{create.pc.basis}}, \code{\link{create.fdata.basis}} o  \code{\link{create.basis}}.\cr
\code{basis.b} is a list of basis  for represent each functional beta parameter. If \code{basis.x} is a list of functional principal components basis (see \code{\link{create.pc.basis}} or \code{\link{pca.fd}}) the argument \code{basis.b} is ignored. 
}
\value{
Return \code{glm} object plus:\cr
  \item{formula}{ formula.}
  \item{data}{ List that containing the variables in the model.}
  \item{group}{ Factor of length \emph{n}}
  \item{group.est}{ Estimated vector groups}
  \item{prob.classification}{ Probability of correct classification by group.}  
  \item{prob.group}{ Matrix of predicted class probabilities. For each functional point shows the probability of each possible group membership.} 
  }
\references{
Ramsay, James O., and Silverman, Bernard W. (2006), \emph{ Functional Data Analysis}, 2nd ed., Springer, New York. \cr
McCullagh and Nelder (1989), \emph{Generalized Linear Models} 2nd ed. Chapman and Hall. \cr
Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S}, New York: Springer. 
%Wood (2001) mgcv:GAMs and Generalized Ridge Regression for R. R News 1(2):20-25 
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:  \code{\link{fregre.glm}}.\cr
%Alternative method: \code{\link{classif.gsam}} and \code{\link{classif.gkam}}.
}
\examples{
data(phoneme)
mlearn<-phoneme[["learn"]]
glearn<-phoneme[["classlearn"]]
mtest<-phoneme[["test"]]
gtest<-phoneme[["classtest"]]
dataf<-data.frame(glearn)
dat=list("df"=dataf,"x"=mlearn)
a1<-classif.glm(glearn~x, data = dat)
newdat<-list("x"=mtest)
p1<-predict.classif(a1,newdat)
table(gtest,p1)
sum(p1==gtest)/250
}
\keyword{classif}