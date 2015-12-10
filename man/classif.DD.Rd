\name{classif.DD}
\Rdversion{1.1}
\alias{classif.DD}
\title{ DD-Classifier Based on DD-plot for G levels and h dephts }
\description{
Fits Nonparametric Classification Procedure Based on DD--plot (depth-versus-depth plot) for G groups.}
\usage{
classif.DD(group,fdataobj,depth="FM",classif="glm",w,
par.classif=list(),par.depth=list(),
control=list(verbose=FALSE,draw=TRUE,col=NULL,alpha=.25))
}



\arguments{
  \item{group}{ Factor of length \emph{n} with \emph{G} levels.}
  \item{fdataobj}{ \code{\link{data.frame}}, \code{\link{fdata}} or list with the multivariate, functional or both covariates respectively.}     
\item{depth}{ Character vector specifying the type of depth functions to use, see \code{Details}. }     
\item{classif}{ Character vector specifying the type of classifier method to use, see \code{Details}.}
  \item{w}{ Optional case weights, weights for each value of \code{depth} argument, see \code{Details}.}
	\item{par.depth}{ List of parameters for \code{depth} function.}	  
  \item{par.classif}{ List of parameters for \code{classif} procedure.}	      
  \item{control}{  List of parameters for controlling the process. If \code{draw=TRUE} print DD-plot of two samples based on data depth.
  \code{col}, The colors for points in DD--plot.  \code{alpha}, the alpha transparency used in the background of DD--plot, a number in [0,1].}  
  }
  \details{
Make the group classification of a training dataset using DD-classifier estimation in the following steps.\cr  
\enumerate{
		\item  The function computes the selected \code{depth} measure of the points in \code{fdataobj} w.r.t. a subsample of each G level
   group.  The user can be specify the parameters for depth function in \code{par.depth}.\cr

(i) Type of depth function from functional data, see \code{\link{Depth}}: 
   \itemize{
		\item  \code{"FM"}: Fraiman and Muniz depth.
		\item  \code{"mode"}: h--modal depth.
%    \item  \code{"MB"}: Modified Band depth.		
		\item  \code{"RT"}: random Tukey depth.
		\item  \code{"RP"}: random project depth.
		\item  \code{"RPD"}: double random project depth.
		}
(ii) Type of depth function from multivariate functional data, see \code{\link{Depth.pfdata}}:
\itemize{
  	\item  \code{"FMp"}: Fraiman and Muniz depth with common support. Suppose that all p--fdata objects have the same support (same rangevals),  see \code{\link{depth.FMp}}.
		\item  \code{"modep"}: h--modal depth using a p--dimensional metric, see \code{\link{depth.modep}}.
%\code{\link{metric.dist}} function is used to compute the distances between the rows of a data matrix (as \code{\link{dist}} function.  
		\item  \code{"RPp"}: random project depth using a p--variate depth with the projections,\cr
     see \code{\link{depth.RPp}}.
			}
If the procedure requires to compute a distance such as in  \code{"knn"} or \code{"np"} classifier or \code{"mode"} depth, the user must use a proper distance function:  \code{\link{metric.lp}} for functional data and \code{\link{metric.dist}} for multivariate data.\cr

(iii) Type of depth function from multivariate data, see \code{\link{Depth.Multivariate}}:
 \itemize{
   \item \code{"SD"}: Simplicial depth (for bivariate data).  
   \item \code{"HS"}: Half-space depth.
   \item \code{"MhD"}:  Mahalanobis dept.
   \item \code{"RD"}: random projections depth.
   \item \code{"LD"}: Likelihood depth.  
      }
      
      
\item The function calculates the misclassification rate based on data depth computed in step (1) using the following classifiers.

  \itemize{
  	\item   \code{"MaxD"}: Maximum depth.
    \item   \code{"DD1"}: Search the best separating polynomial of degree 1.
    \item   \code{"DD2"}: Search the best separating polynomial of degree 2.%Polynomial classifier with 2 degrees.
    \item   \code{"DD3"}: Search the best separating polynomial of degree 3.%Polynomial classifier with 3 degrees.
		\item   \code{"glm"}: Logistic regression is computed using Generalized Linear Models\cr  \code{\link{classif.glm}}.
    \item   \code{"gam"}: Logistic regression is computed using Generalized Additive Models\cr \code{\link{classif.gsam}}. 
  	\item  \code{"lda"}:  Linear Discriminant Analysis is computed using \code{\link{lda}}.
\item  \code{"qda"}:  Quadratic Discriminant Analysis is computed using \code{\link{qda}}.
		\item   \code{"knn"}:  k-Nearest Neighbour classification is computed using \code{\link{classif.knn}}. 
		\item   \code{"np"}:  Non-parametric Kernel classifier is computed using \code{\link{classif.np}}. 
	}
The user can be specify the parameters for classifier function in \code{par.classif} such as the smoothing parameter \code{par.classif[[``h'']]}, if \code{classif="np"} or the k-Nearest Neighbour \code{par.classif[[``knn'']]}, if \code{classif="knn"}.	

In the case of polynomial classifier (\code{"DD1"}, \code{"DD2"} and \code{"DD3"}) uses the original procedure proposed by Li et al. (2012), rotating the DD-plot (to exchange abscise and ordinate).  The procedure extend to multi-class problems by incorporating the method of \emph{majority voting} in the case of polynomial classifier and the method \emph{One vs the Rest} in the logistic case (\code{"glm"} and \code{"gam"}).

  }
 }
 
 
 
\value{
  \item{group.est}{ Estimated vector groups by classified method selected.}
%  \item{group.estMD}{ Estimated vector groups by majority depth classified method.}  
  \item{misclassification}{ Probability of misclassification.}  
  \item{prob.classification}{ Probability of correct classification by group level.} 
  \item{dep}{ Data frame with the depth of the curves for functional data (or points for multivariate data) in \code{fdataobj} w.r.t. each \code{group} level.}	  
  \item{depth}{ Character vector specifying the type of depth functions used.}
  \item{par.depth}{ List of parameters for \code{depth} function.}	  
  \item{classif}{ Type of classifier used.}   
  \item{par.classif}{ List of parameters for \code{classif} procedure.}
  \item{w}{ Optional case weights.}
  \item{fit}{ Fitted object by \code{classif} method using the depth as covariate.}  
}

 

\references{
Li, J., P.C., Cuesta-Albertos, J.A. and Liu, R. \emph{DD--Classifier: Nonparametric Classification Procedure Based on DD-plot}. 
Journal of the American Statistical Association (2012), Vol. 107, 737--753. 

Cuesta-Albertos, J.A., Febrero-Bande, M. and Oviedo de la Fuente, M. \emph{The DDG-classifier in the functional setting}. Submitted. 
}

\author{
This version was created by  Manuel Oviedo de la Fuente and Manuel Febrero Bande and includes the original version for polynomial classifier created by Jun Li, Juan A. Cuesta-Albertos and Regina Y. Liu.
}

\seealso{
See Also as \code{\link{predict.classif.DD}}
}
\examples{

\dontrun{

# DD-classif for functional data
data(tecator)
ab=tecator$absorp.fdata
ab1=fdata.deriv(ab,nderiv=1)
ab2=fdata.deriv(ab,nderiv=2)
gfat=factor(as.numeric(tecator$y$Fat>=15))

# DD-classif for p=1 functional  data set
out01=classif.DD(gfat,ab,depth="mode",classif="np")
out02=classif.DD(gfat,ab2,depth="mode",classif="np")
# DD-plot in gray scale
ctrl<-list(draw=T,col=gray(c(0,.5)),alpha=.2)
out02bis=classif.DD(gfat,ab2,depth="mode",classif="np",control=ctrl)

# 2 depth functions (same curves) 
out03=classif.DD(gfat,list(ab2,ab2),depth=c("RP","mode"),classif="np")
# DD-classif for p=2 functional data set
ldata<-list("ab"=ab2,"ab2"=ab2)
# Weighted version 
out04=classif.DD(gfat,ldata,depth="mode",classif="np",w=c(0.5,0.5))
# Model version
out05=classif.DD(gfat,ldata,depth="mode",classif="np")
# Integrated version (for multivariate functional data)
out06=classif.DD(gfat,ldata,depth="modep",classif="np")

# DD-classif for multivariate data
data(iris)
group<-iris[,5]
x<-iris[,1:4]
out10=classif.DD(group,x,depth="RP",classif="lda")
summary.classif(out10)
out11=classif.DD(group,list(x,x),depth=c("MhD","RP"),classif="lda")
summary.classif(out11)

# DD-classif for functional data: G levels 
data(phoneme)
mlearn<-phoneme[["learn"]]
glearn<-as.numeric(phoneme[["classlearn"]])-1
out20=classif.DD(glearn,mlearn,depth="FM",classif="glm")
out21=classif.DD(glearn,list(mlearn,mlearn),depth=c("FM","RP"),classif="glm")
summary.classif(out20)
summary.classif(out21)
}

}
\keyword{classif}