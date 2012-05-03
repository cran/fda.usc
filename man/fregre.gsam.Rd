\name{fregre.gsam}
\Rdversion{1.1}
\alias{fregre.gsam}
\title{ Fitting Functional Generalized Spectral Additive Models}
\description{
	Computes functional regression between functional (and non functional) explanatory variables and scalar response using basis representation.
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

%When using functional data derived recommend using a number of basis to represent beta lower than the number of basis used to represent the functional data. 
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
Ramsay, James O., and Silverman, Bernard W. (2006), \emph{ Functional Data Analysis}, 2nd ed., Springer, New York. \cr
Muller HG and Stadtmuller U. (2005). \emph{Generalized functional linear models.}  Ann. Statist.33 774-805.\cr
Ramsay, James O., and Silverman, Bernard W. (2006), \emph{ Functional Data Analysis}, 2nd ed., Springer, New York. \cr
McCullagh and Nelder (1989), \emph{Generalized Linear Models} 2nd ed. Chapman and Hall. \cr
Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S}, New York: Springer. \cr
Wood (2001) \emph{mgcv:GAMs and Generalized Ridge Regression for R}. R News 1(2):20-25

}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
\seealso{
See Also as:  \code{\link{predict.fregre.gsam}} and \code{\link{summary.gam}}.\cr
Alternative methods: \code{\link{fregre.glm}}.
}
\examples{

data(tecator)
x=tecator$absorp.fdata
#tecator$y$Fat=floor(tecator$y$Fat)
x.d1<-fdata.deriv(x)
tt=x[["argvals"]]
dataf=as.data.frame(tecator$y)
nbasis.x=11;nbasis.b=7
basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)
f=Fat~s(Protein)+s(x)
basis.x=list("x"=basis1,"x.d1"=basis1)
basis.b=list("x"=basis2,"x.d1"=basis2)
ldata=list("df"=dataf,"x"=x,"x.d1"=x.d1)
res=fregre.gsam(Fat~Water+s(Protein)+x+s(x.d1),ldata,family=gaussian(),basis.x=basis.x,basis.b=basis.b)
res2=fregre.gsam(Fat~te(Protein,k=3)+x,data=ldata,family=gaussian())

res=fregre.gsam(Fat~s(Water)+x.d1,data=ldata,family=gaussian(),basis.x=basis.x,basis.b=basis.b)

#gam(Fat~s(Protein,k=3),data=dataf)

#res$beta.l$x$coefs
res2=fregre.gsam(Fat~te(Protein)+x,data=ldata,family=gaussian())
#res2$beta.l$x$coefs
a<-predict.fregre.gsam(res)


# dropind basis pc
basis.pc0=create.pc.basis(x,c(1:7))
basis.pc1=create.pc.basis(x.d1,c(1:7))
basis.x=list("x"=basis.pc0,"x.d1"=basis.pc1)
ldata=list("df"=dataf,"x"=x,"x.d1"=x.d1)
#f.pc=Fat~Proteins+s(x) #
#res.pc=fregre.gsam(Fat~Protein+s(x),data=ldata,family=gaussian(),basis.x=basis.x,basis.b=basis.b)
#summary(res.pc)

########################
# clase fd
data(tecator)
ind<-1:129
x=tecator$absorp.fdata[ind]
y=tecator$y$Fat[ind]
tt=x[["argvals"]]
dataf=as.data.frame(tecator$y[ind,])

x.fd<-fdata2fd(x)
x.fd1<-fdata.deriv(x,class.out="fd")
x.fd2<-fdata.deriv(x,class.out="fd",nderiv=2)
basis.b=list("x.fd"=basis2,"x.fd1"=basis2)
ldata=list("df"=dataf,"x.fd"=x.fd,"x.fd1"=x.fd1,"x.fd2"=x.fd2)
f0=Fat~+x.fd+s(x.fd1)
res.pc=fregre.gsam(f0,data=ldata,family=gaussian())
res.pc2=fregre.gsam(f0,data=ldata,family=gaussian(),basis.b=basis.b)

#pca.fd
basis.pc0=pca.fd(x.fd,3)
basis.pc1=pca.fd(x.fd1,3)
basis.pc2=pca.fd(x.fd2,3)
basis.x=list("x.fd"=basis.pc0,"x.fd1"=basis.pc1)
f0=Fat~+x.fd+s(x.fd1)+s(Protein)
res.fd=fregre.gsam(f0,data=ldata,basis.x=basis.x)


data(tecator)
x=tecator$absorp.fdata
tecator$y$Fat<-ifelse(tecator$y$Fat>20,1,0)
x.d1<-fdata.deriv(x)
tt=x[["argvals"]]
dataf=as.data.frame(tecator$y)
nbasis.x=11;nbasis.b=7
basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)
f=Fat~s(Protein)+s(x)
basis.x=list("x"=basis1,"x.d1"=basis1)
basis.b=list("x"=basis2,"x.d1"=basis2)
ldata=list("df"=dataf,"x"=x,"x.d1"=x.d1)
res=fregre.gsam(Fat~Protein+s(x),ldata,family=binomial(),basis.x=basis.x,basis.b=basis.b)
res2=fregre.gsam(Fat~s(Protein)+x,data=ldata,family=binomial())

res=fregre.gsam(f,data=ldata,family=poisson(),basis.x=basis.x,basis.b=basis.b)
#res=fregre.gsam(Fat~Protein+x,data=ldata,family=gaussian()) no functional variables
res=fregre.gsam(Fat~s(Protein)+s(Water)+s(x)+x.d1,data=ldata,family=gaussian())
summary(res)
plot(res)
res=fregre.gsam(Fat~Protein+s(x),data=ldata,family=gaussian())
res=fregre.gsam(Fat~s(Protein)+x,data=ldata,family=gaussian())
res=fregre.gsam(f,ldata,family=gaussian(),basis.x=basis.x,basis.b=basis.b)

}
\keyword{regression}