\name{classif.gkam}
\Rdversion{1.1}
\alias{classif.gkam}
\title{ Classification Fitting Functional Generalized Kernel Additive Models}
\description{
	Computes functional classification using functional explanatory variables using backfitting algorithm.
}
\usage{
classif.gkam(formula,family = binomial(),data,weights = rep(1, nobs),
 par.metric = NULL,par.np=NULL, offset=NULL,
 control = list(maxit = 100,epsilon = 0.001, trace = FALSE ,inverse="solve"),...)

}
\arguments{
 \item{formula}{ an object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under \code{Details}.}
  \item{data}{ List that containing the variables in the model. }
 \item{family}{ a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.) }
  \item{weights}{ weights}
  \item{par.metric}{ List of arguments by covariable to pass to the \code{metric} function by covariable.}
  \item{par.np}{  List of arguments to pass to the \code{fregre.np.cv} function} 
  \item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting.}	     
  \item{control}{ a list of parameters for controlling the fitting process, by default: maxit, epsilon, trace and inverse.}
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{ 
The first item in the \code{data} list  is called \emph{"df"} and  is a data frame  with the response, as \code{\link{glm}}.\cr
Functional covariates of class  \code{fdata} are introduced in the following items in the \code{data} list.\cr
}
\value{
Return \code{gam} object plus:\cr
  \item{formula}{ formula.}
  \item{data}{ List that containing the variables in the model.}
  \item{group}{ Factor of length \emph{n}}
  \item{group.est}{ Estimated vector groups}
  \item{prob.classification}{ Probability of correct classification by group.}  
  \item{prob.group}{ Matrix of predicted class probabilities. For each functional point shows the probability of each possible group membership.}    
}
\references{
Febrero-Bande M. and Gonzalez-Maneiga W. (2011). \emph{Generalized Additive Models for Functional Data}. Recent Advances in functional Data Analyisis and Related Topics. Physica-Velag. \cr
McCullagh and Nelder (1989), \emph{Generalized Linear Models} 2nd ed. Chapman and Hall. \cr
Opsomer J.D. and Ruppert D.(1997). \emph{Fitting a bivariate additive model by local polynomial regression}.Annals of Statistics, \code{25}, 186-211.
}

\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:  \code{\link{fregre.gkam}}.\cr
Alternative method: \code{\link{classif.glm}}.
}
\examples{

# Time-consuming: selection of 2 levels 
data(phoneme)
mlearn<-phoneme[["learn"]][1:100]
glearn<-as.numeric(phoneme[["classlearn"]][1:100])
dataf<-data.frame(glearn)
dat=list("df"=dataf,"x"=mlearn)
a1<-classif.gkam(glearn~x,data=dat)
summary(a1)
mtest<-phoneme[["test"]][1:100]
gtest<-as.numeric(phoneme[["classtest"]][1:100])
newdat<-list("x"=mtest)
p1<-predict.classif(a1,newdat)
table(gtest,p1)

}
\keyword{classif}