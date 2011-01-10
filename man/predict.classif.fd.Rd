\name{predict.classif.fd}
\Rdversion{1.1}
\alias{predict.classif.fd}
\title{ Predicts from a fitted classif.fd object.}
\description{
  Classifier of Functional Data by kernel method using Functional Data Object of class \code{classif.fd}. 
}
\usage{
\method{predict}{classif.fd}(object,new.fdataobj,probs.group=FALSE, \dots)
%predict.classif.fd(object, new.fdataobj,probs.group=FALSE,...)
}
\arguments{
  \item{object}{ Object \code{object} estimated by: k nearest neighbors method \code{classif.knn.fd}, kernel method  \code{classif.kernel.fd} or  kernel method  with fucntional basis representation \code{classif.kernel.fb}.  } 
  \item{new.fdataobj}{ New functional explanatory data of \code{fdata} class.}  
  \item{probs.group}{ =TRUE, for each functional data shows the probability of each group membership.}  
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{
Returns the groups or classes predicted using a previously trained model.
  
}
\value{
 \item{group.pred}{Vector of groups or classes predicted}
 \item{prob.group}{For each functional data shows the probability of each group membership.}  
 }
\references{
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametricc functional data analysis.}
Springer Series in Statistics, New York.

Ferraty, F. and Vieu, P. (2006). \emph{NPFDA in practice}. Free access on line at
\url{http://www.lsp.ups-tlse.fr/staph/npfda/}

Ramsay, James O., and Silverman, Bernard W. (2006), \emph{ Functional Data Analysis}, 2nd ed., Springer, New York. 

Ramsay, James O., and Silverman, Bernard W. (2002), \emph{ Applied Functional Data Analysis}, Springer, New York. 

}
\author{Febrero-Bande, M. and Oviedo de la Fuente, M.}
\seealso{ 
See also \code{\link{classif.knn.fd}}, \code{\link{classif.kernel.fd}},  \cr
\code{\link{classif.kernel.fb}} and \code{\link{summary.classif.fd}} 
}
\examples{


data(phoneme)
mlearn<-phoneme[["learn"]]
glearn<-phoneme[["classlearn"]]
mtest<-phoneme[["test"]]
gtest<-phoneme[["classtest"]]

#	ESTIMATION knn
out1=classif.knn.fd(mlearn,glearn,knn=c(3,5,7))
summary.classif.fd(out1)
#	PREDICTION knn
pred1=predict.classif.fd(out1,mtest,TRUE)
table(pred1$group.pred,gtest)
MisclasPredict <- sum(pred1$group.pred != gtest)/length(gtest)
MisclasPredict

#	ESTIMATION kernel
h=2^(0:5)
out2=classif.kernel.fd(mlearn,glearn,h=h)
summary.classif.fd(out2)
#	PREDICTION kernel
pred2=predict.classif.fd(out2,mtest,TRUE)
table(pred2$group.pred,gtest)
MisclasPredict <- sum(pred2$group.pred != gtest)/length(gtest)
MisclasPredict

#	ESTIMATION kernel functional basis
vbasis=seq(7,37,by=10)
out3=classif.kernel.fb(mlearn,glearn,h=h,type.basis="fourier",
par.basis=list(nbasis=vbasis))
summary.classif.fd(out3)
#	PREDICTION kernel functional basis
pred3=predict.classif.fd(out3,mtest,TRUE)
table(pred3$group.pred,gtest)
MisclasPredict <- sum(pred3$group.pred != gtest)/length(gtest)
MisclasPredict
}
\keyword{classif}