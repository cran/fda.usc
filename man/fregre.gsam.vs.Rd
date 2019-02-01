\name{fregre.gsam.vs}
\Rdversion{1.1}
\alias{fregre.gsam.vs}
\title{ Variable Selection using Functional Additive Models}
\description{
Computes functional GAM model between functional covariate \eqn{(X^1(t_1),\cdots,X^{q}(t_q))}{(X(t_1),...,X(t_q))} (and non functional covariate \eqn{(Z^1,...,Z^p)}{(Z1,...,Zp)}) and scalar response \eqn{Y}.
	
This function is an extension of the functional generalized linear regression models: \code{\link{fregre.glm}} where the \eqn{E[Y|X,Z]} is related to the linear prediction \eqn{\eta} via a link function \eqn{g(\cdot)}{g(.)} with integrated smoothness estimation by the smooth functions \eqn{f(\cdot)}{f(.)}. 

\deqn{E[Y|X,Z])=\eta=g^{-1}(\alpha+\sum_{i=1}^{p}f_{i}(Z^{i})+\sum_{k=1}^{q}\sum_{j=1}^{k_q}{f_{j}^{k}(\xi_j^k)})}{E[Y|X,Z]=\eta=g^{-1}(\alpha+\sum_i  f_i(Z_{i})+\sum_k^q\sum_{j=1}^{k_q}{f_j^k(\xi_j^k)})}
where \eqn{\xi_j^k}{\xi_j^k} is the coefficient of the basis  function expansion of \eqn{X^k}, (in PCA analysis \eqn{\xi_j^k}{\xi_j^k} is the score of the \eqn{j}-functional PC of \eqn{X^k}.

The smooth functions \eqn{f(\cdot)}{f(.)} can be added to the right hand side of the formula to specify that the linear predictor depends on smooth functions of predictors using smooth terms \code{\link{s}} and \code{\link{te}} as in  \code{\link{gam}} (or linear functionals of these as \eqn{Z\beta} and \eqn{\big<X(t),\beta\big>}{< X(t),\beta(t) >} in \code{\link{fregre.glm}}). 
}
\usage{
fregre.gsam.vs(data, y, x, alpha, type.basis = "pc", ncomp, kbs,
            #criterio="sp",
           dcor.min=.1, par.model, xydist, trace = FALSE,
           CV = TRUE, ncomp.fix = FALSE, smooth = TRUE)
}

\arguments{
 \item{data}{ List that containing the variables in the model.
 "df" element is a data.frame with the response and scalar covariates (numeric and factors variables are allowed). Functional covariates of class  \code{fdata} or \code{fd} are introduced in the following items in the \code{data} list.} 
 \item{y}{ caracter string with the name of the scalar response variable}
 \item{x}{ caracter string vector with the name of the scalar and functional potential covariates.}                     
 \item{alpha}{alpha value uese to thes the null hypothesis for the test of independence among covariate X and residual e. By default is 0.05.}
 \item{type.basis}{ Type of basis used, by default principal component basis "pc".} 
 \item{ncomp}{ maximum number of basis elements (only used in functional covariates).}
 \item{kbs}{ the dimension of the basis used to represent the smooth term. The default depends on the number of variables that the smooth is a function of.}
 \item{dcor.min}{lower threshold for the variable X to be considered. X is discarded if the distance correlation \eqn{R(X,e)< dcor.min} (e is the residual).}
 %\item{family}{ a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.) }
 %\item{data}{ List that containing the variables in the model. } 
 \item{par.model}{ Model parameters.}
 \item{xydist}{list with the matrices of distances of each variable (all potential covariates and the response) with itself.}
 \item{trace}{ Interactive Tracing and Debugging of Call.}
 \item{CV}{ TRUE, Cross-validation (CV) is done.}	
 \item{ncomp.fix}{if TRUE, a number of basis element is fixed in \code{ncomp}. If FALSE,  the funcion selects the number of PC components between the \code{ncomp}.}
 \item{smooth}{if TRUE, a smooth estimate is made for all covariates included in the model (less for factors). The model is adjusted with the estimated variable linearly or smoothly. If the models are equivalent, the model is adjusted with the linearly estimated variable.}
% \item{\dots}{ Further arguments passed to or from other methods.}
}
\value{
%"model"       "gof"         "i.predictor" "xydist"      "dcor" 
Return a list with the follow elements:
\item{model}{object corresponding to the estimated additive mdoel using the selected variables. Same output as the\code{\link{fregre.gsam}} function.}
\item{gof}{The goodness of fit for each step of VS algorithm.}
\item{i.predictor}{Vector  with 1 if the variable is selected, 0 otherwise.}
\item{dcor}{the value of distance correlation for each pontential covariate and the residual of the model in each step.}
}
\references{
Febrero-Bande, M., Gonz\'{a}lez-Manteiga, W. and Oviedo de la Fuente, M. Variable selection in functional additive regression models, (2018).  Computational Statistics, 1-19. DOI: \url{https://doi.org/10.1007/s00180-018-0844-5}
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
\dontrun{ 
data(tecator)
x=tecator$absorp.fdata
x1<-fdata.deriv(x)
x2<-fdata.deriv(x,nderiv=2)
y=tecator$y$Fat
xcat0<-cut(rnorm(length(y)),4)
xcat1<-cut(tecator$y$Protein,4)
xcat2<-cut(tecator$y$Water,4)
ind <- 1:129
dat    <- data.frame("Fat"=y, x1$data, xcat1, xcat2)
ldat <- list("df"=dat[ind,],"x"=x[ind,],"x1"=x1[ind,],"x2"=x2[ind,])
# 3 functionals (x,x1,x2), 3 factors (xcat0, xcat1, xcat2)
# and 100 scalars (impact poitns of x1) 

# Time consuming
 res.gam1<-fregre.gsam.vs("Fat",data=ldat)
 summary(res.gam1$model)
 res.gam1$i.predictor

covar <- c("xcat0","xcat1","xcat2","x","x1","x2")
res.gam2<-fregre.gsam.vs(y="Fat", x=covar, data=ldat)
summary(res.gam2$model)
res.gam2$i.predictor

# Prediction like fregre.gsam() 
newldat <- list("df"=dat[-ind,],"x"=x[-ind,],"x1"=x1[-ind,],
                "x2"=x2[-ind,])
pred.gam1<-predict.fregre.gsam(res.gam1$model,newldat)
pred.gam2<-predict.fregre.gsam(res.gam2$model,newldat)
plot(dat[-ind,"Fat"],pred.gam1)
points(dat[-ind,"Fat"],pred.gam2,col=2)
}
}
\keyword{regression}