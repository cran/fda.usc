\name{predict.classif}
\Rdversion{1.1}
\alias{predict.classif}
\title{ Predicts from a fitted classif object.}
\description{
  Classifier of Functional Data by kernel method using Functional Data Object of class \code{classif}. 
}
\usage{
\method{predict}{classif}(object, new.fdataobj,type="class",\dots)
%predict.classif(object, new.fdataobj,type="class",...)
}
\arguments{
  \item{object}{ Object \code{object} estimated by: k nearest neighbors method \code{classif.knn}, kernel method  \code{classif.kernel}.} 
  \item{new.fdataobj}{ New functional explanatory data of \code{fdata} class.}  
  \item{type}{ !=''class'', for each functional data shows the probability of each group membership.}  
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
See also \code{\link{classif.knn}}, \code{\link{classif.kernel}} and \code{\link{summary.classif}} .
}
\examples{

data(phoneme)
mlearn<-phoneme[["learn"]]
glearn<-phoneme[["classlearn"]]

#	ESTIMATION knn
out1=classif.knn(glearn,mlearn,knn=c(3,7))
summary.classif(out1)
#	PREDICTION knn
mtest<-phoneme[["test"]][c(40:60,140:160,201:210)]
gtest<-phoneme[["classtest"]][c(40:60,140:160,201:210)]
pred1=predict.classif(out1,mtest)
table(pred1,gtest)
MisclasPredict <- sum(pred1 != gtest)/length(gtest)
MisclasPredict

##	ESTIMATION kernel 
# h=2^(0:5)
# using metric distances computed in classif.knn
# out2=classif.kernel(glearn,mlearn,h=h,metric=out1$mdist)
# summary.classif(out2)
#	PREDICTION kernel
# pred2=predict.classif(out2,mtest)
# table(pred2,gtest)
# MisclasPredict <- sum(pred2!= gtest)/length(gtest)
# MisclasPredict

}
\keyword{classif}