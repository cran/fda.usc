predict.fregre.fd<-function(object,new.fdataobj=NULL,se.fit=FALSE,...){
if (is.null(object)) stop("No fregre.fd object entered")
#if (is.null(new.fdataobj)) stop("No newx entered")
if (is.null(new.fdataobj)) return(object$fitted.values)
if (object$call[[1]]=="gam")  return(predict(object,new.fdataobj,...))
if (object$call[[1]]=="lm")  return(predict(object,new.fdataobj,...)  )
if (object$call[[1]]=="glm")  return(predict(object,new.fdataobj,...))
#falta np (no funcional) y basis (no funcional)
#print("bbbbb")
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
np <- ncol(new.fdataobj)
 if (is.null(rownames(newx)))         rownames(newx) <- 1:nn
 if (object$call[[1]]=="fregre.pc" || object$call[[1]]=="fregre.ppc") {
  Z<- inprod.fdata(fdata.cen(new.fdataobj,object$fdata.comp$mean)[[1]],object$fdata.comp$rotation)
  colnames(Z)<-names(object$lm$coefficients[-1])
  XX<-data.frame(Z)   
   if (object$rn==0) return(predict.lm(object=object$lm,newdata=XX,se.fit=se.fit,x=TRUE,y=TRUE,...))
   else {
    a1<-object$coefficients[1]*rep(1,len=nrow(newx))
    object$beta.est$data<-matrix(object$beta.est$data,nrow=1)
    b1<-inprod.fdata(fdata.cen(new.fdataobj,object$fdata.comp$mean)[[1]],object$beta.est)#/(ncol(newx)-1)
    yp<- a1+b1      
    yp<-drop(yp)
    XX2<-cbind(rep(1,len=nn),Z)
#    yp.lm<-predict.lm(object=object$lm,newdata=XX,se.fit=se.fit,...)  
    names(yp) <-gg
    if (se.fit) {
     se.fit<-sqrt(rowSums((XX2 %*%object$Vp*XX2)))
     names(se.fit)<-gg
     return(list("fit"=yp,"se.fit"=se.fit))
    }
    else      return(yp)   
    }
 }
else{
  if (object$call[[1]]=="fregre.pls" || object$call[[1]]=="fregre.ppls") {
  newXcen<-fdata.cen(new.fdataobj,object$fdata.comp$mean)[[1]]
  if (object$fdata.comp$norm)  {
    sd.X <- sqrt(apply(object$fdataobj$data, 2, var))
    newXcen$data<- newXcen$data/(rep(1, nn) %*% t(sd.X))
  }
  Z<- inprod.fdata(newXcen,object$fdata.comp$rotation) 
  colnames(Z)<-names(object$lm$coefficients[-1])  
  XX<-data.frame(Z)  
  return(predict.lm(object=object$lm,newdata=XX,se.fit=se.fit,...))
#  a1<-object$coefficients[1]*rep(1,len=nrow(newx))  
#  object$beta.est$data<-matrix(object$beta.est$data,nrow=1)
#  b1<-inprod.fdata(newXcen,object$beta.est)#/(ncol(newx)-1)
#  yp<- a1+b1      
#  yp<-drop(yp)        
#  names(yp) <-gg  
#  if (se.fit) { 
#     p1<-predict.lm(object=object$lm,newdata=XX,se.fit=se.fit,...)
#     XX2<-cbind(rep(1,len=nn),Z) 
#     se.fit<-sqrt(rowSums((XX2 %*%object$Vp*XX2)))       
#       names(se.fit)<-gg    
# return(list(p1$fit,yp=yp,p1$se.fit,se.fit)) 
#       }
#  else {       return(yp) } 
 } 
else {
 if (object$call[[1]]=="fregre.basis" || object$call[[1]]=="fregre.basis.cv"){
  x=newx
  basis.x=object$basis.x.opt             #
  xcen<-fdata.cen(new.fdataobj,object$mean)[[1]]
	x.fd=Data2fd(argvals=tt,y=t(xcen$data),basisobj=basis.x)
  C=t(x.fd$coefs)
  if (is.vector(object$b.est)) object$b.est<-matrix(object$b.est,ncol=1,nrow=length(object$b.est))
  Z<-C%*%object$J
#  colnames(Z)<-names(object$lm$coefficients[-1])  
#  p1<-predict.lm(object=object$lm,newdata=data.frame(Z),se.fit=se.fit,...)  
  yp=object$a.est* rep(1,len=nn) + Z%*%object$b.est
   if (isfdata) {
    return(fdata(yp,y$argvals,y$rtt,y$names))
    }
  if (se.fit) {
     XX2<-cbind(rep(1,len=nn),Z) 
     se.fit<-sqrt(rowSums((XX2 %*%object$Vp*XX2)))       
 #    names(se.fit)<-gg    
     return(list("fit"=yp,"se.fit"=se.fit)) 
 }
 else {       return(yp) } 
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
#print(y.mat)
  yp=H%*%y.mat    
 if (isfdata) {
    return(fdata(yp,y$argvals,y$rtt,y$names))
    }
 else{
  yp<-drop(yp)
  names(yp)<-gg
  if (se.fit) {
   se.fit<- H%*%y.mat^2-yp^2
   return(list("fit"=yp,"se.fit"=sqrt(se.fit)))
  } 
  else (return(yp))
 }   
 }     }
}     }
#rownames(yp)=rownames(newx)
return(yp)
}


