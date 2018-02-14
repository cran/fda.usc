# #############################################
#  Author: Manuel Oviedo de la Fuente
#  November, 2017
#  R code for manuscrit: "Determining optimum wavelengths for leaf water 
#  content estimation from reflectance: a distance correlation approach"


##################################
# Find local maxima of x
# Returns the set of locations of the local maxima
localMaxima <- function(x) {
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}


##################################
##################################
# Arguments
# y: name of the response variable
# covar: vector with the names of the predictor variables
# data: data frame 
# We compute the Distance correlation (and  t-test) of multivariate and functional independence 
# (wrapper functions of energy package).
# tol: Tolerance values for the selection criteria. The variables whose distance correlation 
#      values are greater than the tolerance are selected.
# pvalue: p-value of the t-test. 
# plot: logical value, if TRUE plots the distance correlation curve for  each covariate in multivariate case
# and in each discretization points (argvals) in the functional case.
# smo: logical value, If true the distance correlation function is smoothed.
#     using b-spline representation with a suitable number of basis elements.
# local.pc: not used in the manuscript
# ceros.pc: not used in the manuscript
# verbose: print iterative and relevant steps of the procedure.
# Return a list of two elements:
# a) cor: the value of distance correlation for eahc covariate
# b) maxLocal: index or locations of local maxima distance correlations
LMDC.select <- function(y, covar, data, tol = .06, pvalue = .05,
                        plot = FALSE, local.dc = TRUE, 
                        smo = FALSE, verbose = FALSE){
  
  yy <- data[,y]
  xx<-data[,covar,drop=F]
  nn<-ncol(xx)  
  dc<-numeric(nn)
  for (i in 1:nn){
    a<-dcor.xy(xx[,i,drop=F],yy)
    dc[i]<-a$estimate
    if (verbose)   cat(i,a$estimate,a$p.value,a$p.value<pvalue & a$estimate>tol,"\n")    
  }  
  if (smo) {
    nbase<-ifelse(nn<50,floor(nn/2),floor(nn^(4/5)))
    #print(nbase)
    dc1<-fdata2fd(fdata(dc),nbasis=nbase)
    dc2<-fdata(dc1,1:nn)$data[1,]
  }
  else dc2<-dc
  regre<-TRUE #REGRESION O CLASIFICACION
  if (is.factor(yy)) regre<-FALSE
  #if(local.pc)   pc<-fregre.pc.cv(fdata(xx),yy)[[1]]
  max.pc1<-max.pc2<-max.pc3<-NULL
  if (local.dc)  maxLocal<-localMaxima(dc2) 
  #if (local.pc)  max.pc1<-localMaxima(pc$beta.est$data[1,])
  #if (ceros.pc)  max.pc3<- c(1,which(diff(pc$beta.est$data[1,]>0)!=0))
  maxLocal<-unique(c(maxLocal,max.pc1,max.pc3))
  maxLocal2<-intersect(which(dc2>tol),maxLocal)
  xorder<-order(dc2[maxLocal2],decreasing =TRUE)
  dc2<-dc2[maxLocal2][xorder]
  maxLocal2<-maxLocal2[xorder]
  if (plot){
    par(mfrow=c(1,1))
    plot(dc2)
    lines(dc,col=4)
    abline(v=maxLocal,col=2)
    abline(v=maxLocal2,col=3,lwd=2)    
  }
  nvar<-length(maxLocal2)
  return(list(dcor=dc,maxLocal=maxLocal2))
}
#######################################################

#######################################################
# Argurments
# y: name of the response variable
# covar: vector with the names of the predictor variables
# data: data frame for model estimation  (data train)
# newdata: (optional) new data frame for model prediction (data test)
# We compute the Distance correlation (and  t-test) of multivariate and functional independence 
# (wrapper functions of energy package).
# tol: Tolerance values for the selection criteria. The variables whose distance correlation 
#      values are greater than the tolerance are selected.
# pvalue: p-value of covariate significance (used in "lm" and "gam" methods).

# method: character, name of the regression method called. Model options
#   "lm":  step-wise lm regression model (lm function, stats package)
#   "gam": step-wise gam regression model (gam function, mgcv package).
#   "svm": Support vector machine (svm function, e1071 package).
#   "knn": k-nearest neighbor regression (knnn.reg function,F NN package).
#   "lars": Least Angle Regression using Lasso (lars function, lars package).
#   "glmnet": Lasso and Elastic-Net Regularized Generalized Linear Models (glmnet and cv.glmnet function, glmnet package).
#   "rpart": Recursive partitioning for regression a (rpart function, rpart package).

# par.method: list of arguments related with the regression method called.
# plot: logical value, not implemented
# verbose: print iterative and relevant steps of the procedure.

# Return a list of two elements:
# a) model: final estimated model
# b) xvar: not implemented
# c) pred: (optional) vector with the predicted values
# d) nvar: number of covariates

#######################################################

#######################################################
LMDC.regre <- function(y,covar,data,newdata,pvalue=.05,
                       method="lm", par.method = NULL,
                       plot=FALSE,verbose=FALSE){
  edf<-Inf
  nvar <- length(covar)
  if (missing(newdata)) pred <- FALSE
  else pred <- TRUE
  if (is.null(covar)){
    return(list(model=lm(data[,y]~ 1,data=data), xvar=xvar, pred=rep(mean(data[,y]),len=nrow(newdata))) )   
  }
  pred0 <- NULL
  nvar <- length(covar)
  #xnames<-colnames(xx)[maxLocal2]
  xnames <- covar
  xvar<-NULL
  
  ff<-paste(y,"~1")
  if (method == "lm"){
    #model0<-lm(formula(ff),data=data)
    par.method$formula <- formula(ff)
    par.method$data <- data
    model0<-do.call(method,par.method)
    for (i in 1:nvar) {
      if (verbose)   print(i) 
      xentra<-xnames[i]
      xvar2<- c(xvar,xentra)
      ff<-paste(y,"~",paste(xvar2,collapse="+"),collapse="")
      if (verbose)  {print("lm");    print(1);print(ff)}
      par.method$formula <- formula(ff)
      model <- do.call(method,par.method)
      if (rev(summary(model)$coefficients[,"Pr(>|t|)"])[1]<pvalue){
        if (verbose)  { print("entra");print(xentra);print(summary(model))}
        model0 <- model
        xvar <- xvar2  
      }     
    } # fin for
    edf <- summary(model0)$df[1]
    nvar<-edf-1
  } #fin lm
  
  if (method == "gam"){
    if (!is.null(par.method$k)) {
      ik<-which(names(par.method)=="k")
      print(par.method)
      k<-par.method$k
      par.method<-par.method[-ik]
      }
    else k <- 4
    ff<-as.formula(ff)
    model0 <- gam(ff,data=data)
    par.method2 <- list("formula"=ff,"data"= data)
    par.method<-c(par.method2,par.method)
   for (i in 1:nvar) {
      if (verbose)   print(i) 
      xentra<-xnames[i]
      xvar2<- c(xvar,xentra)   
      ff<-as.formula(paste(y,"~",paste("s(",xvar2,",k=",k,")",collapse="+"),collapse=""))
      par.method$formula<-ff
      if (verbose)  {print("gam");    print(1);print(ff)}
      model <- do.call(method,par.method)
        if (rev(summary(model)$s.table[,"p-value"])[1]<pvalue){
        if (verbose)  { print("entra");print(xentra);print(summary(model))}
        model0 <- model
        xvar <- xvar2  
      }     
    } # fin for
    edf<- sum(model0$edf)
    nvar <- ncol(model0$model)-1
  } #fin Gam
  ######################################################  
  if (method == "svm"){   
    if (is.null(par.method)) par.method=list("cost"=100,"gamma"=1,"kernel"="linear")
    par.method$x <- data[,covar,drop=F]
    par.method$y <- data[,"y"]
    model0 <- do.call(method,par.method)
    
  }
  if (method == "rpart"){
    ff<-as.formula(paste(y,"~",paste(covar,collapse="+"),collapse=""))
    #model0 <- rpart(ff, data = data)
    par.method$formula <- ff
    par.method$data <- data
    model0 <- do.call(method,par.method)
  }
  if (method == "knn"){
     par.method$train <- data[,covar,drop=F]
    par.method$test <- newdata[,covar,drop=F]
    par.method$y <- data[,"y"]
    model0 <- do.call("knn.reg",par.method)
    pred0<-model0$pred
  }
  if ( method =="lars") {
    if (is.null(par.method)) 
      par.method= list(type="lasso",normalize=FALSE,intercept = TRUE,use.Gram=FALSE)    
    x0<-as.matrix(data[,covar])
    par.method$x <- x0
    par.method$y <- data[,"y"]
    model0 <- do.call(method, par.method)
    templam <- range(model0$lambda)
    lambda<-seq(templam[1], templam[2], length=200)
    cv <- do.call("cvlars", 
          list(x=x0, y=data[,"y"],K=10,lambda=lambda, trace = FALSE,
               intercept = TRUE,normalize=FALSE, type="lasso",use.Gram=F))
    #cv <- cvlars(x=x0, y=data[,"y"],K=10,lambda=lambda, trace = FALSE,
    #               intercept = TRUE,normalize=FALSE, type="lasso",use.Gram=F)
    minl<-lambda[which.min(cv)]
    pred0 <- do.call("predict.lars",
        list("object" = model0, "newx" = as.matrix(newdata[,covar]), 
             "s" = minl,"type"= "fit", "mode"= "lambda"))$fit
    #pred0 <- predict.lars( model0 ,newx=as.matrix(newdata[,covar]), s=minl,type="fit",mode="lambda")$fit
    cv <- do.call("lars::cv.lars",list("x" = x0, "y" = data[,"y"],"K" = 10))
    #cv <- lars::cv.lars(x0, data[,"y"],K=10)
#    pred0 <- predict.lars( cv  ,newx=as.matrix(newdata[,covar]), s=minl,type="fit",mode="lambda")$fit
    ideal_l1_ratio <- cv$index[which.max(cv$cv - cv$cv.error <= min(cv$cv))]
    #obj <- lars::lars(x0, data[,"y"])
    obj <- do.call("lars::lars",list("x"=x0, "y"=data[,"y"]))
    scaled_coefs <- scale(obj$beta, FALSE, 1 / obj$normx)
    l1 <- apply(X = scaled_coefs, MARGIN = 1, FUN = function(x) sum(abs(x)))
    coef(obj)[which.max(l1 / do.call("tail",list("x"=l1, "n"=1)) > ideal_l1_ratio),]
    pred0 <- do.call("predict.lars",list("object"= model0,
          "newx"=as.matrix(newdata[,covar]), 
          "s"=minl,"type"="fit","mode"="lambda"))$fit
    #pred0 <- predict.lars( model0 ,newx=as.matrix(newdata[,covar]), s=minl,type="fit",mode="lambda")$fit
    nvar <- sum( coef(obj)[which.max(l1 / do.call("tail",list("x"=l1, "n"=1)) > ideal_l1_ratio),]>0)
    }
  
  if ( method =="glmnet") {
    #x should be a matrix with 2 or more columns 
    x0 <- as.matrix(data[,covar,drop=F])
    newx0 <- as.matrix(newdata[,covar,drop=F])
    if (ncol(x0)==1){
      x0<-cbind(x0,1:nrow(x0))
      newx0<-cbind(newx0,1:nrow(newx0))
    }
    if (is.null(par.method)) 
      par.method= list(family="gaussian", standardize=TRUE, nfolds=10)
    
    par.method$x <- x0
    par.method$y <- data[,"y"]
    model0 <- do.call("cv.glmnet",par.method)    
    pred0<-predict(model0,newx=newx0, s=model0$lambda.min)
    
    #min.scale<-which.min(model0$cvm)
    #min.lambda<-model0$lambda.min
    #temp<-coef( model0$glmnet.fit, s=min.lambda)
    #temp<-temp[-1]
    #ntheta1<-rep(0,length(coef(model0)))
    #ntheta1[which(abs(temp) > 0)]<-temp[which(abs(temp) > 0)]
    #edf <- length(which(abs(ntheta1)>0))
    c<-coef(model0 ,s='lambda.min',exact=TRUE)
    inds<-which(c!=0)
    variables<-row.names(c)[inds]
    variables<-setdiff(variables,'(Intercept)')
    nvar<- length(variables)
    edf <-nvar +1
  }
  if ( method =="nnet") {
    if (is.null(par.method)) 
      par.method<-list( size = 5, rang = .1,decay = 5e-6, maxit =1000,linout=T)
    par.method$x <- data[,covar,drop=F]
    par.method$y <- data[,"y"]
    model0 <- do.call("nnet",par.method)    
     }
  if ( method =="mars") {
    #if (is.null(par.method)) 
     # par.method<-list( size = 5, rang = .1,decay = 5e-6, maxit =1000,linout=T)
    par.method$x <- data[,covar,drop=F]
    par.method$y <- data[,"y"]
    model0 <- do.call("mars",par.method)    
  }
  if (method == "npreg"){
    #ff<-as.formula(paste(y,"~",paste(covar,collapse="+"),collapse=""))
    #model0 <- rpart(ff, data = data)
    #par.method$formula <- ff
    #par.method$data <- data
    par.method$txdat<-data[,covar,drop=F]
    par.method$tydat<-data[,"y"]
    model0 <- do.call(method,par.method)
  }
  if (method == "flam"){
    #if (is.null(par.method)) 
   # par.method=list("alpha" = 0.75, n.fold = 2)
    par.method$x <- data[,covar,drop=F]
    par.method$y <- data[,"y"]
    model0 <- do.call("flamCV",par.method)
    alpha <- model0$alpha
    lambda <- model0$lambda.cv
    pred0<- predict(model0$flam.out, new.x =newdata[,covar,drop=F],
                     lambda = lambda, alpha = alpha)
  }
  if (method == "novas"){
    #if (is.null(par.method)) 
    # par.method=list("alpha" = 0.75, n.fold = 2)
    par.method$COVARIATES<- data[,covar,drop=F]
    par.method$Responses  <- data[,"y"]
    model0 <- do.call("novas",par.method)
    pred0<- predict(model0, newdata[,covar,drop=F])
    nvar <- edf <- length(model0$model)
  }
  if ( method =="cosso") {
    #x should be a matrix with 2 or more columns 
    x0 <- as.matrix(data[,covar,drop=F])
    newx0 <- as.matrix(newdata[,covar,drop=F])
    if (ncol(x0)==1){
      x0<-cbind(x0,1:nrow(x0))
      newx0<-cbind(newx0,1:nrow(newx0))
    }
    model0<-do.call("cosso",list("x"=x0,"y"=data[,y],"family"="Gaussian"))
    xvar<-do.call("predict.cosso",
          list("object"=model0,"M"=2,"type"="nonzero"))
    pred0<-do.call("predict.cosso",list("object"=model0,"xnew"=newx0,
                                        "M"=2,"type"="fit"))
    nvar <- length(xvar)
    edf <-nvar+1
  }
  if (method != "knn" & method != "cosso" & method != "lars" & method != "novas" & method != "glmnet" & method != "flam"){
    #print(method);    print(pred)
    if (pred) pred0 <- predict(model0, newdata=newdata[,covar,drop=F])
    else pred0 <- NULL
  }
  return(list(model=model0, xvar=xvar, pred=pred0,edf=edf,nvar=nvar))     
}  
################################################################################
################################################################################
