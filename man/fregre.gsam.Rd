\name{fregre.gsam}
\Rdversion{1.1}
\alias{fregre.gsam}
\title{ Fitting Functional Generalized Spectral Additive Models}
\description{
Computes functional GAM model between functional covariate \eqn{(X^1(t_1),\cdots,X^{q}(t_q))}{(X(t_1),...,X(t_q))} (and non functional covariate \eqn{(Z^1,...,Z^p)}{(Z1,...,Zp)}) and scalar response \eqn{Y}.
	
This function is an extension of the functional generalized linear regression models: \code{\link{fregre.glm}} where the \eqn{E[Y|X,Z]} is related to the linear prediction \eqn{\eta} via a link function \eqn{g(\cdot)}{g(.)} with integrated smoothness estimation by the smooth functions \eqn{f(\cdot)}{f(.)}. 

\deqn{g(E[Y|X,Z])=\eta=\alpha+\sum_{i=1}^{p}f_{i}(Z^{i})+\sum_{k=1}^{q}\sum_{j=1}^{k_q}{f_{j}^{k}(\xi_j^k)}}{g(E[Y|X,Z])=\eta=\alpha+\sum_i  f_i(Z_{i})+\sum_k^q\sum_{j=1}^{k_q}{f_j^k(\xi_j^k)}}
where \eqn{\xi_j^k}{\xi_j^k} is the coefficient of the basis  function expansion of \eqn{X^k}, (in PCA analysis \eqn{\xi_j^k}{\xi_j^k} is the score of the \eqn{j}-functional PC of \eqn{X^k}.

The smooth functions \eqn{f(\cdot)}{f(.)} can be added to the right hand side of the formula to specify that the linear predictor depends on smooth functions of predictors using smooth terms \code{\link{s}} and \code{\link{te}} as in  \code{\link{gam}} (or linear functionals of these as \eqn{Z\beta} and \eqn{\big<X(t),\beta\big>}{< X(t),\beta(t) >} in \code{\link{fregre.glm}}). 
}
\usage{
fregre.gsam(formula, family = gaussian(),data=list(),
weights=NULL,basis.x=NULL,basis.b=NULL,CV=FALSE,...) 
}
\arguments{
 \item{formula}{ an object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under \code{Details}.}
 \item{family}{ a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.) }
 \item{data}{ List that containing the variables in the model. } 
 \item{weights}{ weights}
 \item{basis.x}{ List of basis for functional explanatory data estimation.}
 \item{basis.b}{ List of basis for functional beta parameter estimation.}
 \item{CV}{ =TRUE, Cross-validation (CV) is done.}	
 \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{ 
The first item in the \code{data} list  is called \emph{"df"} and  is a data frame  with the response and non functional explanatory variables, as \code{\link{gam}}.\cr

Functional covariates of class  \code{fdata} or \code{fd} are introduced in the following items in the \code{data} list.\cr
\code{basis.x} is a list of basis  for represent each functional covariate. The basis object can be created by the function: \code{\link{create.pc.basis}}, \code{\link{pca.fd}} \code{\link{create.pc.basis}}, \code{\link{create.fdata.basis}} o  \code{\link{create.basis}}.\cr
\code{basis.b} is a list of basis  for represent each functional beta parameter. If \code{basis.x} is a list of functional principal components basis (see \code{\link{create.pc.basis}} or \code{\link{pca.fd}}) the argument \code{basis.b} is ignored. 

}
\value{
Return \code{gam} object plus:\cr
\item{basis.x}{ Basis used for \code{fdata} or \code{fd} covariates.}
\item{basis.b}{ Basis used for beta parameter estimation.}
\item{beta.l}{ List of estimated beta parameter of functional covariates.}
\item{data}{ List that containing the variables in the model.}
\item{formula}{ formula.}
\item{CV}{ \$y.pred predicted response by cross-validation.}	
}
\references{
Muller HG and Stadtmuller U. (2005). \emph{Generalized functional linear models.}  Ann. Statist.33 774-805.

Wood (2001) \emph{mgcv:GAMs and Generalized Ridge Regression for R}. R News 1(2):20-25.

Ramsay, James O., and Silverman, Bernard W. (2006), \emph{ Functional Data Analysis}, 2nd ed., Springer, New York. 

Ramsay, James O., and Silverman, Bernard W. (2006), \emph{ Functional Data Analysis}, 2nd ed., Springer, New York. 

McCullagh and Nelder (1989), \emph{Generalized Linear Models} 2nd ed. Chapman and Hall. 

Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S}, New York: Springer.
}
\note{If the formula only contains a non functional explanatory variables (multivariate covariates), the function compute a standard  \code{\link{glm}} procedure.}
\author{
Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}
\seealso{
See Also as:  \code{\link{predict.fregre.gsam}} and \code{\link{summary.gam}}.\cr
Alternative methods: \code{\link{fregre.glm}} and \code{\link{fregre.gkam}}.
}
\examples{

data(tecator)
x=tecator$absorp.fdata
x.d1<-fdata.deriv(x)
tt<-x[["argvals"]]
dataf=as.data.frame(tecator$y)
nbasis.x=11;nbasis.b=5
basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)
f=Fat~s(Protein)+s(x)
basis.x=list("x"=basis1,"x.d1"=basis1)
basis.b=list("x"=basis2,"x.d1"=basis2)
ldata=list("df"=dataf,"x"=x,"x.d1"=x.d1)
res=fregre.gsam(Fat~Water+s(Protein)+x+s(x.d1),ldata,family=gaussian(),
basis.x=basis.x,basis.b=basis.b)
res

\dontrun{
res2=fregre.gsam(Fat~te(Protein,k=3)+x,data=ldata,family=gaussian())
summary(res2)

##  dropind basis pc
basis.pc0=create.pc.basis(x,c(2,4,7))
basis.pc1=create.pc.basis(x.d1,c(1:3))
basis.x=list("x"=basis.pc0,"x.d1"=basis.pc1)
ldata=list("df"=dataf,"x"=x,"x.d1"=x.d1)  
res.pc=fregre.gsam(f,data=ldata,family=gaussian(),
basis.x=basis.x,basis.b=basis.b)
summary(res.pc)

## fd class object
y=tecator$y$Fat
dataf=as.data.frame(tecator$y)
x.fd<-fdata2fd(x)
x.fd1<-fdata.deriv(x,class.out="fd")
x.fd2<-fdata.deriv(x,class.out="fd",nderiv=2)
basis.b=list("x.fd"=basis2,"x.fd1"=basis2)
ldata=list("df"=dataf,"x.fd"=x.fd,"x.fd1"=x.fd1,"x.fd2"=x.fd2)
res.pc=fregre.gsam(f,data=ldata,family=gaussian(),basis.b=basis.b)

## pca.fd
basis.pc0=pca.fd(x.fd,3)
basis.pc1=pca.fd(x.fd1,3)
basis.pc2=pca.fd(x.fd2,3)
basis.x=list("x.fd"=basis.pc0,"x.fd1"=basis.pc1)
f0=Fat~+x.fd+s(x.fd1)+s(Protein)
res.fd=fregre.gsam(f0,data=ldata,basis.x=basis.x)

##  Binomial family
x=tecator$absorp.fdata
tecator$y$Fat<-ifelse(tecator$y$Fat>20,1,0)
x.d1<-fdata.deriv(x)
dataf=as.data.frame(tecator$y)
ldata=list("df"=dataf,"x"=x,"x.d1"=x.d1)
res.bin=fregre.gsam(Fat~Protein+s(x),ldata,family=binomial())
}
}
\keyword{regression}