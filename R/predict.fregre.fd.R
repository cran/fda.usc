predict.fregre.fd<-function(object,new.fdataobj=NULL,...){
if (is.null(object)) stop("No fregre.fd object entered")
if (is.null(new.fdataobj)) stop("No newx entered")
if (!is.fdata(new.fdataobj)) new.fdataobj=fdata(new.fdataobj,object$fdataobj[["argvals"]],object$fdataobj[["rangeval"]],object$fdataobj[["names"]])
y=object$y
isfdata<-is.fdata(y)
if (!isfdata) {
   if (is.vector(y)) y.mat<-matrix(y,ncol=1)
   }
else y.mat<-y$data
gg<-1:nrow(new.fdataobj)
nas<-apply(new.fdataobj$data,1,count.na)
if (any(nas)) {
   bb<-!nas
   cat("Warning: ",sum(nas)," curves with NA are omited\n")
   new.fdataobj$data<-new.fdataobj$data[bb,]
   gg<-gg[bb]
   }
newx<-new.fdataobj[["data"]]
tt<-new.fdataobj[["argvals"]]
rtt<-new.fdataobj[["rangeval"]]
nn <- nrow(new.fdataobj)
 if (is.null(rownames(newx)))         rownames(newx) <- 1:nn
 if (object$call[[1]]=="fregre.pc") {
 a1<-object$coefficients[1]*rep(1,len=nrow(newx))
 object$beta.est$data<-matrix(object$beta.est$data,nrow=1)
# b2<-new.fdataobj$data%*%t(object$beta.est$data)
 b1<-inprod.fdata(fdata.cen(new.fdataobj,object$fdata.comp$mean)[[1]],object$beta.est)#/(ncol(newx)-1)
 yp<- a1+b1
 }
 else {
  if (object$call[[1]]=="fregre.pls") {
  a1<-object$coefficients[1]*rep(1,len=nrow(newx))
  object$beta.est$data<-matrix(object$beta.est$data,nrow=1)
#  b2<-new.fdataobj$data%*%t(object$beta.est$data)
  b1<-inprod.fdata(fdata.cen(new.fdataobj,object$fdata.comp$mean)[[1]],object$beta.est)#/(ncol(newx)-1)
  yp<- a1+b1
 }
 else {
 if (object$call[[1]]=="fregre.basis" || object$call[[1]]=="fregre.basis.cv"){
  x=newx
  basis.x=object$basis.x.opt             #
# 	xmean=apply(x,2,mean)
#  xcen=sweep(x,2,xmean,FUN="-")
  xcen<-fdata.cen(new.fdataobj,object$mean)[[1]]
	x.fd=Data2fd(argvals=tt,y=t(xcen$data),basisobj=basis.x)
  C=t(x.fd$coefs)
  if (is.vector(object$b.est)) object$b.est<-matrix(object$b.est,ncol=1,nrow=length(object$b.est))
  yp=object$a.est* rep(1,len=nn) + C%*%object$J%*%object$b.est
  yp <- matrix(yp,ncol=1,nrow=nn)
  }
 else {
 if (object$call[[1]]=="fregre.np" || object$call[[1]]=="fregre.np.cv"){
# if (is.vector(newx))  newx <- t(as.matrix(newx))
 x=object$fdataobj
 h=object$h.opt
 n = nrow(x)
 nn = nrow(newx)
 np <- ncol(x)
# if (n != (length(y)))         stop("ERROR IN THE DATA DIMENSIONS")
 if (is.null(rownames(newx)))         rownames(newx) <- 1:nn
 par.S<-object$par.S
 bs=as<-list()
 Ker=object$Ker
 #   par.metric<-list()
   par.metric<-attr(object$mdist,"par.metric")
   par.metric[["fdata1"]]<-new.fdataobj
   par.metric[["fdata2"]]<-x
#   parm<-attr(object$mdist,"par.metric")
#   lenpar<-length(parm)
#names(par.metric[3:(2+lenpar)])<-names(attr(object$mdist,"par.metric"))
  a1<-attr(object$mdist,"call")
  a2<-attr(object$par.S,"call")
  nmdist <- do.call(a1,par.metric)
  par.S$tt<-nmdist
  par.S$cv=FALSE
  H<-do.call(a2,par.S)
  yp=H%*%y.mat
 if (isfdata) {
    return(fdata(yp,y$argvals,y$rtt,y$names))
    }
 }     }
 }     }
#rownames(yp)=rownames(newx)
yp<-drop(yp)
names(yp)<-gg
return(yp)
}


