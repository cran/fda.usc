      
\name{fda.usc-package}
\Rdversion{1.0.0}
\alias{fda.usc-package}
\alias{fda.usc}
\docType{package}
\title{Functional Data Analysis and Utilities for Statistical Computing (fda.usc)}
\description{
This package carries out exploratory and descriptive analysis of functional data exploring its most important features: such as depth measurements or functional outliers detection, among others. \cr
It also helps to explain and model the relationship between a dependent variable and independent (regression models) and make predictions. Methods for supervised or unsupervised classification of a set of functional data regarding a feature of the data are also included. Finally,  it can perform analysis of variance model (ANOVA) for functional data.
 
Sections of fda.usc-package: \cr
 \cr
\tabular{ll}{
\tab A.- Functional Data Representation  \cr
\tab B.- Functional Outlier Detection \cr
\tab C.- Functional Regression Model	\cr
\tab D.- Functional Supervised Classification	\cr
\tab E.- Functional Non-Supervised Classification	\cr
\tab F.- Functional ANOVA	\cr
\tab G.- Auxiliary functions:\cr
}


  A.- Functional Data Representation	\cr
The functions included in this section allow to define, transform, manipulate and represent a functional dataset in many ways including derivatives, non-parametric kernel methods or basis representation.\cr

\tabular{ll}{
\tab \code{\link{fdata}}  \cr
\tab \code{\link{plot.fdata}} \cr
\tab \code{\link{fdata.deriv}}  \cr
\tab \code{\link{CV.S}} \cr
\tab \code{\link{GCV.S}} \cr
\tab \code{\link{min.np}} \cr
\tab \code{\link{min.basis}} \cr
\tab \code{\link{S.NW}} \cr
\tab \code{\link{S.LLR}} \cr
\tab \code{\link{S.basis}} \cr
\tab \code{\link{Var.e}} \cr
\tab \code{\link{Var.y}} \cr
}

  B.- Functional Depth and Functional Outlier Detection 	\cr

The functional data depth calculated by the different depth  functions implemented that could be use as a measure of centrality or outlyingness.\cr
\cr
B.1-Depth methods \code{\link{Depth}}:\cr
\tabular{ll}{
\tab \code{\link{depth.FM}} \cr
\tab \code{\link{depth.mode}} \cr
\tab \code{\link{depth.RP}} \cr
\tab \code{\link{depth.RT}} \cr
\tab \code{\link{depth.RPD}} \cr
\tab \code{\link{Descriptive}} \cr
}
B.2-Functional Outliers detection methods:\cr
\tabular{ll}{
\tab \code{\link{outliers.depth.trim}} \cr
\tab \code{\link{outliers.depth.pond}} \cr
\tab \code{\link{outliers.thres.lrt}} \cr
\tab \code{\link{outliers.lrt}} \cr
}
  C.- Functional Regression	Models\cr

C.1. Functional explanatory covariate  and scalar response\cr
The functions included in this section allow the estimation of different functional regression models with a scalar response and a single functional explicative covariate.\cr
\tabular{ll}{
\tab \code{\link{fregre.pc}}  \cr
\tab \code{\link{fregre.pc.cv}}  \cr
\tab \code{\link{fregre.pls}}  \cr
\tab \code{\link{fregre.pls.cv}}  \cr
\tab \code{\link{fregre.basis}}  \cr
\tab \code{\link{fregre.basis.cv}}  \cr
\tab \code{\link{fregre.np}}  \cr
\tab \code{\link{fregre.np.cv}}  \cr
}

C.2. Test for the functional linear model (FLM) with scalar response.\cr	
\tabular{ll}{
\tab \code{\link{flm.Ftest}}, F-test for the FLM with scalar response  \cr
\tab \code{\link{flm.test}}, Goodness-of-fit test for the FLM with scalar response  \cr
\tab \code{\link{PCvM.statistic}},  PCvM statistic for the FLM with scalar response  \cr
}

C.3. Functional and non functional explanatory covariates.\cr
The functions in this section extends those regression models in previous section in several ways.  \cr
Semifunctional partial linear regression \code{\link{fregre.plm}}  is an extension of functional nonparameric regression  \code{\link{fregre.np}}  allowing  include non-functional variables.\cr
Functional  linear regression \code{\link{fregre.lm}}, functional generalized linear regression \code{\link{fregre.glm}} and 
functional generalized spectral additive model \code{\link{fregre.gsam}}  are an extensions of  \code{\link{fregre.basis}} and \code{\link{fregre.pc}} allowing  include  more than one functional variable and other non-functional variables, as \code{\link{lm}} or \code{\link{glm}} functions.\cr

\tabular{ll}{
\tab \code{\link{fregre.plm}}  \cr
\tab \code{\link{fregre.lm}}  \cr
\tab \code{\link{fregre.glm}}  \cr
\tab \code{\link{fregre.gsam}}  \cr
\tab \code{\link{fregre.gkam}}  \cr
}


C.4. Functional explanatory covariate and  response\cr
The functions included in this section allow the estimation of  functional regression models with a functional response and a single functional explicative covariate.\cr
\tabular{ll}{
\tab \code{\link{fregre.basis.fr}}  \cr
}

  D.- Functional Supervised Classification	\cr
This section allows the estimation of the groups in a training set of functional data \code{fdata} class by different nonparametric methods of supervised classification. Once these classifiers have been trained, they can be used to predict on new functional data.\cr	
\cr
Package allows the estimation of the groups in a training set of functional data by different methods of supervised classification. \cr
\cr

\tabular{ll}{
\tab \code{\link{classif.knn}}  \cr
\tab \code{\link{classif.kernel}} \cr
\tab \code{\link{classif.glm}} \cr
\tab \code{\link{classif.gsam}} \cr
\tab \code{\link{classif.gkam}} \cr
\tab \code{\link{classif.DD}} \cr
\tab \code{\link{classif.depth}} \cr
}

	E.- Functional Non-Supervised Classification	\cr
This section allows the estimation of the groups in a  functional data set \code{fdata} class by kmeans method. \cr
\tabular{ll}{
\tab	\code{\link{kmeans.fd}} \cr
}

	F.- Functional ANOVA	\cr
\tabular{ll}{
\tab	\code{\link{anova.onefactor}} \cr
\tab	\code{\link{anova.RPm}} \cr
\tab	\code{\link{anova.hetero}} \cr
}

	G.- Utilities and auxiliary functions:\cr
	\tabular{ll}{
\tab \code{\link{fdata.bootstrap}} \cr
\tab \code{\link{fdata2fd}}  \cr
\tab \code{\link{fdata2pc}} \cr
\tab \code{\link{fdata2pls}} \cr
\tab \code{\link{summary.fdata.comp}} \cr
\tab \code{\link{cond.F}} \cr
\tab \code{\link{cond.quantile}} \cr
\tab \code{\link{cond.mode}} \cr
\tab \code{\link{FDR}} \cr
\tab \code{\link{Kernel}} \cr
\tab \code{\link{Kernel.asymmetric}} \cr
\tab \code{\link{Kernel.integrate}} \cr
\tab \code{\link{metric.lp}} \cr 
\tab \code{\link{metric.kl}} \cr 
\tab \code{\link{metric.hausdorff}} \cr 
\tab \code{\link{metric.dist}} \cr 
\tab \code{\link{semimetric.NPFDA}} \cr 
\tab \code{\link{semimetric.basis}} \cr 
}
	}
\details{ 
\tabular{ll}{
Package: \tab fda.usc\cr
Type: \tab Package\cr
Version: \tab 1.2.0\cr
Date: \tab 2014-05-05\cr
License: \tab GPL-2 \cr
LazyLoad: \tab yes\cr
Url: \url{ http://eio.usc.es/pub/MAESFE/}, \url{http://eio.usc.es/pub/gi1914/}
}
}
\author{
\emph{Authors:}
Manuel Febrero Bande \email{manuel.febrero@usc.es} and Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}

\emph{Contributors:}
 Pedro Galeano, Alicia Nieto-Reyes, Eduardo Garcia-Portugues \email{eduardo.garcia@usc.es} and STAPH group \url{http://www.lsp.ups-tlse.fr/staph/}
 %pedro.galeano@uc3m.es,  alicia.nieto@unican.es STAPH group maintaining the page http://www.lsp.ups-tlse.fr/staph/
 

\emph{Maintainer:} Manuel Oviedo de la Fuente \email{manuel.oviedo@usc.es}
}


\references{
Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
}
\keyword{ package }