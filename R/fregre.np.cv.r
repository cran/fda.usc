fregre.np.cv=function(fdataobj,y,h=NULL,Ker=AKer.norm,metric=metric.lp,
type.CV = GCV.S,type.S=S.NW,par.CV=list(trim=0),...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
x<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
   C<-match.call()
   mf <- match.call(expand.dots = FALSE)
   m<-match(c("x", "y","h","Ker","metric","type.CV","type.S","par.CV"),names(mf),0L)
#   if (is.vector(x))         x <- t(as.matrix(x))
   n = nrow(x)
   np <- ncol(x)
   if (n != (length(y)))         stop("ERROR IN THE DATA DIMENSIONS")
   if (is.null(rownames(x)))         rownames(x) <- 1:n
   if (is.null(colnames(x)))         colnames(x) <- 1:np
   mdist=metric(x,x,...)
types=FALSE
if (is.null(h)) {
if (m[7]==0){
     mdist2=mdist
     diag(mdist2)=Inf
     h0<-apply(mdist2,1,min,na.rm=TRUE)
     h.max=max(h0)
     h.med=median(h0)
     q.min=quantile(mdist2,probs=0.025,na.rm=TRUE)
     q.max=quantile(mdist2,probs=0.2,na.rm=TRUE)
     h.min=max(q.min,h.med)
     h.max=max(q.max,h.max)
     if (m[4]!=0) {
       i=1
       lenC=length(C)
    while (i<=lenC) {
        if (C[[i]]!="AKer.norm") {
                h.max=min(q.max,h.max)
                h.min=min(q.min,h.med)
                i=lenC+1
                 }
        else i=i+1        }
}
else {h.max=min(q.max,h.max)
                h.min=min(q.min,h.med)}
h = seq(h.min, h.max, len = 51)
}
else if (C[[m[7]]]=="S.KNN") {
          h.min=2
          h.max=floor(quantile(1:n,probs=0.10,na.rm=TRUE,type=4))
          h=h.min:h.max
     }
else {
     mdist2=mdist
     diag(mdist2)=Inf
     h0<-apply(mdist2,1,min,na.rm=TRUE)
     h.max=max(h0)
     h.med=median(h0)
     q.min=quantile(mdist2,probs=0.025,na.rm=TRUE)
     q.max=quantile(mdist2,probs=0.2,na.rm=TRUE)
     h.min=max(q.min,h.med)
     h.max=max(q.max,h.max)
     if (m[4]!=0) {
       i=1
       lenC=length(C)
  while (i<=lenC) {
        if (C[[i]]!="AKer.norm") {

                h.max=min(q.max,h.max)
                h.min=min(q.min,h.med)
                i=lenC+1
                 }
        else i=i+1        }
}
else {h.max=min(q.max,h.max)
                h.min=min(q.min,h.med)}
h = seq(h.min, h.max, len = 51)
}}
else {if   (any(h<=0)) stop("Error: Invalid range for h")}
  lenh <- length(h)
gcv=cv.error <- array(NA, dim = c(lenh))
  for (i in 1:lenh) {
     H =type.S(mdist,h[i],Ker,cv=FALSE)
#     gcv[i] <- type.CV(y, H,trim=par.CV$trim,draw=par.CV$draw,...)
     par.CV$S<-H
     par.CV$y<-y
     gcv[i]<- do.call(type.CV,par.CV)
#     e2=y-H%*%y
#     cv.error[i]=sum(e2^2)/(n-traza(H))
   }
if (all(is.infinite(gcv))) print(" Warning: Invalid range for h")
   l = which.min(gcv)
   h.opt <- h[l]
   if (h.opt==min(h)) cat(" Warning: h.opt is the minimum value of bandwidths
   provided, range(h)=",range(h),"\n")
   else if (h.opt==max(h)) cat(" Warning: h.opt is the maximum value of bandwidths
   provided, range(h)=",range(h),"\n")
   H=matrix(NA,ncol=n,nrow=n)
if (h.opt>0) {   H =type.S(mdist,h.opt,Ker,cv=FALSE)  }
else {
     h.opt=-h.opt
     H =type.S(mdist,h.opt,cv=FALSE)
     }
 	yp=H%*%y
  e=matrix(y,ncol=1)-yp
  df=traza(H)
  sr2=sum(e^2)/(n-df)
  ycen=y-mean(y)
	r2=1-sum(e^2)/sum(ycen^2)
	names(gcv)<-h
out<-list("call"=C,"fitted.values"=yp,"H"=H,"residuals"=e,"df"=df,"r2"=r2,"sr2"=sr2,"y"=y,"fdataobj"=fdataobj,"mdist"=mdist,"Ker"=Ker,"metric"=metric,"type.S"=type.S,"gcv"=gcv,"h.opt"=h.opt,"h"=h,"m"=m)
 class(out)="fregre.fd"
return(out)
}
