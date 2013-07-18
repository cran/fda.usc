\name{classif.depth}
\Rdversion{1.1}
\alias{classif.depth}
\title{ Classifier from Functional Data}
\description{
Classification of functional data using maximum depth.}
\usage{
classif.depth(group,fdataobj,newfdataobj,depth="RP",par.depth=list(),CV="none")

}
\arguments{
 \item{group}{ Factor of length \emph{n}}
 \item{fdataobj}{  \code{fdata}, \code{matrix} or \code{data.frame} class object of training set cases.}  
 \item{newfdataobj}{  \code{fdata}, \code{matrix} or \code{data.frame} class object of training set cases.}  
 	\item{depth}{ Type of depth function from functional data: 
		\itemize{
		\item   \code{FM}: Fraiman and Muniz depth.
		\item  \code{mode}: modal depth.
		\item  \code{RT}: random Tukey depth.
		\item  \code{RP}: random project depth.
		\item  \code{RPD}: double random project depth.
		}
		}  
 	\item{par.depth}{ List of parameters for \code{depth}.}	  
  \item{CV}{ =``none'' \code{group.est=group.pred}, =TRUE \code{group.est} is estimated by cross-validation, =FALSE \code{group.est} is estimated.}		
}
%\details{}
\value{ 
  \item{group.est}{ Vector of classes of train sample data.}
  \item{group.pred}{ Vector of classes of test sample data.}  
  \item{prob.classification}{ Probability of correct classification by group.}  
	\item{max.prob}{ Highest probability of correct classification.}
  \item{fdataobj}{ \code{\link{fdata}} class object.}  
  \item{group}{ Factor of length \emph{n}.}  
  %\item{prob.group}{ Matrix of predicted class probabilities. For each functional point shows the probability of each possible group membership.}   
  }
\references{
Cuevas, A., Febrero-Bande, M. and Fraiman, R. (2007). \emph{Robust estimation and classification for functional data via projection-based depth notions.} Computational Statistics 22, 3, 481{-}496.  
}
\author{
Febrero-Bande, M. and Oviedo de la Fuente, M.
}
%\seealso{See Also as \code{\link{classif.DD}}}
\examples{
 
data(phoneme)
mlearn<-phoneme[["learn"]]
mtest<-phoneme[["test"]]
glearn<-phoneme[["classlearn"]]
gtest<-phoneme[["classtest"]]

a1<-classif.depth(glearn,mlearn,depth="RP")
table(a1$group.est,glearn)
table(a1$group.pred,gtest)

# a2<-classif.depth(glearn,mlearn,depth="RP",CV=TRUE)
# a3<-classif.depth(glearn,mlearn,depth="RP",CV=FALSE)
# a4<-classif.depth(glearn,mlearn,mtest,"RP")
# a5<-classif.depth(glearn,mlearn,mtest,"RP",CV=TRUE)     
# table(a5$group.est,glearn)
# table(a5$group.pred,gtest)
# a6<-classif.depth(glearn,mlearn,mtest,"RP",CV=FALSE)
# table(a6$group.est,glearn)
# table(a6$group.pred,gtest)

}
\keyword{classif}