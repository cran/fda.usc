\name{fregre.glm}
\Rdversion{1.1}
\alias{fregre.glm}
\title{ Fitting Functional Generalized Linear Models}
\description{
	Computes functional regression between functional (and non functional) explanatory variables and scalar response using basis representation.
}
\usage{
fregre.glm(formula,family = gaussian(),data,basis.x=NULL,basis.b=NULL,CV=FALSE,...)
}
\arguments{
  \item{formula}{ an object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under \code{Details}.} 
  \item{family}{ a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.) }
  \item{data}{ List that containing the variables in the model. } 
  %\item{weigths}{ weights}
	\item{basis.x}{ List of basis for functional explanatory data estimation.}
	\item{basis.b}{ List of basis for functional beta parameter estimation.}
	\item{CV}{ =TRUE, Cross-validation (CV) is done .}	
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{ 
The first item in the \code{data} list  is called \emph{"df"} and  is a data frame  with the response and non functional explanatory variables, as \code{\link{glm}}.\cr

Functional covariates of class  \code{fdata} or \code{fd} are introduced in the following items in the \code{data} list.\cr
\code{basis.x} is a list of basis  for represent each functional covariate. The basis object can be created by the function: \code{\link{create.pc.basis}}, \code{\link{pca.fd}} \code{\link{create.pc.basis}}, \code{\link{create.fdata.basis}} o  \code{\link{create.basis}}.\cr
\code{basis.b} is a list of basis  for represent each functional beta parameter. If \code{basis.x} is a list of functional principal components basis (see \code{\link{create.pc.basis}} or \code{\link{pca.fd}}) the argument \code{basis.b} is ignored. 

%When using functional data derived recommend using a number of basis to represent beta lower than the number of basis used to represent the functional data. 
}
\value{
Return \code{glm} object plus:\cr
\item{basis.x}{ Basis used for \code{fdata} or \code{fd} covariates.}
\item{basis.b}{ Basis used for beta parameter estimation.}
\item{beta.l}{ List of estimated beta parameter of functional covariates.}
\item{data}{ List that containing the variables in the model.}
\item{formula}{ formula.}
\item{CV}{ \$y.pred predicted response by cross-validation.}	
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
See Also as:  \code{\link{predict.fregre.glm}} and \code{\link{summary.glm}}.\cr
Alternative method if \code{family}=\emph{gaussian}: \code{\link{fregre.lm}}.
}
\examples{

data(tecator)
x=tecator$absorp.fdata
y=tecator$y$Fat
tt=x[["argvals"]]
dataf=as.data.frame(tecator$y)

nbasis.x=11
nbasis.b=7
basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)
 
f=Fat~Protein+x
basis.x=list("x"=basis1)
basis.b=list("x"=basis2)
ldata=list("df"=dataf,"x"=x)
res=fregre.glm(f,family=gaussian(),ldata,basis.x=basis.x,basis.b=basis.b)
summary(res)
}
\keyword{regression}