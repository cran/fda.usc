################################################################################
predict.fregre.fr<-function(object,new.fdataobj=NULL,...){
if (is.null(object)) stop("No fregre.fd object entered")
if (is.null(new.fdataobj)) return(object$fitted.values)
if (object$call[[1]]=="fregre.basis.fr" || object$call[[1]]=="fregre.basis.fr.cv"){
beta.est<-object$coefficients
isfdx<-is.fd(new.fdataobj)
if (isfdx) {
  xcoef<-new.fdataobj$coef
  ncurves<-ncol(xcoef)
  }
else {
 xfdobj<-Data2fd(argvals =new.fdataobj$argvals, y = t(new.fdataobj$data), basisobj = object$basis.s)
 xcoef<-xfdobj$coef
 ncurves<-ncol(xcoef)
 if (any(new.fdataobj$argvals!=object$x$argvals)) stop("Incorrect argvals")
}
H = t(xcoef) %*% object$H 
beta.xest = beta.est %*% t(H)
beta.xfd   = fd(beta.xest, object$basis.t)
yhat = eval.fd(object$argvals.y,object$alpha.est) %*% matrix(1,1,ncurves) + eval.fd(object$argvals.y, beta.xfd)
if (isfdx) {
 fitted.values  <- fd(coef=yhat, basisobj=object$y$basis, fdnames=object$y$fdnames)
}
else {
 fitted.values<-fdata(t(yhat),new.fdataobj$argvals,new.fdataobj$rangeval,new.fdataobj$names)
}
return(fitted.values)
}
else stop("predict is only implemented for fregre.basis.fr output object")
}
################################################################################
