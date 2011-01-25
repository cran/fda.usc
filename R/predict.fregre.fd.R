predict.fregre.fd<-function(object,new.fdataobj=NULL,...){
if (is.null(object)) stop("No fregre.fd object entered")
if (is.null(new.fdataobj)) stop("No newx entered")
if (!is.fdata(new.fdataobj)) new.fdataobj=fdata(new.fdataobj,object$fdataobj[["argvals"]],object$fdataobj[["rangeval"]],object$fdataobj[["names"]])

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
# b2<-newx%*%(object$beta.est[["data"]])/(ncol(newx)-1)
 object$beta.est$data<-matrix(object$beta.est$data,nrow=1)
 b1<-inprod.fdata(new.fdataobj,object$beta.est)#/(ncol(newx)-1)
 yp<- a1+b1
 }
 else {
 if (object$call[[1]]=="fregre.basis" || object$call[[1]]=="fregre.basis.cv"){
  x=newx
  basis.x=object$basis.x.opt             #
 	xmean=apply(x,2,mean)
  xcen=sweep(x,2,xmean,FUN="-")
  class(xcen)="matrix"
	x.fd=Data2fd(argvals=tt,y=t(xcen),basisobj=basis.x)
  C=t(x.fd$coefs)
  if (is.vector(object$b.est)) object$b.est<-matrix(object$b.est,ncol=1,nrow=length(object$b.est))
  yp=object$a.est* rep(1,len=nn) + C%*%object$J%*%object$b.est
  yp <- matrix(yp,ncol=1,nrow=nn)
  }
 else {
 if (object$call[[1]]=="fregre.np" || object$call[[1]]=="fregre.np.cv"){
# if (is.vector(newx))  newx <- t(as.matrix(newx))
 x=object$fdataobj
 y=object$y
 h=object$h.opt
 n = nrow(x)
 nn = nrow(newx)
 np <- ncol(x)
# if (n != (length(y)))         stop("ERROR IN THE DATA DIMENSIONS")
 if (is.null(rownames(newx)))         rownames(newx) <- 1:nn
 bs=as<-list()
 C<-object$call
 m<-object$m
 Ker=object$Ker
 imetric<-m[5]
 if (imetric==0) {
      a1<-metric.lp
      len.metricc<-length(formals(a1))
      vv<-array(0,dim=c(len.metricc))      }
 else {
   a1<-match.fun(C[[imetric]])
   len.metricc<-length(formals(a1))
   vv<-array(0,dim=c(len.metricc))   }
   ii<-imetric+1
if (C[ii]!='NULL()')     {
    ind.m<-3
    while (C[ii]!='NULL()'&& ind.m<=len.metricc) {
    aa<-any(names(C)==names(formals(a1))[ind.m])
    if (aa) {
           vv[ind.m]<-which(names(C)==names(formals(a1)[ind.m]))
           ii<-ii+1
           as[[ind.m]]<-C[[vv[ind.m]]]           }
     else {as[[ind.m]]<-formals(a1)[[ind.m]]}
       ind.m<-ind.m+1
 } }
 as[[2]]<-new.fdataobj
 as[[1]]<-x
 mdist<-do.call(a1,as,...)
 kmdist=Ker(mdist/h)
 a=t(kmdist)%*%matrix(y,ncol=1)
 b=apply(kmdist,2,sum)
 yp=a/b
#   yp=t(object$type.S(mdist,h,object$Ker,cv=FALSE))%*%matrix(y,ncol=1)
 }}
 }
#rownames(yp)=rownames(newx)
yp<-drop(yp)
names(yp)<-gg
return(yp)
}

