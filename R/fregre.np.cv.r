fregre.np.cv=function(fdataobj,y,h=NULL,Ker=AKer.norm,metric=metric.lp,
type.CV = GCV.S,type.S=S.NW,par.CV=list(trim=0),...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
nas<-apply(fdataobj$data,1,count.na)
nas.g<-is.na(y)
if (is.null(names(y))) names(y)<-1:length(y)
if (any(nas) & !any(nas.g)) {
   bb<-!nas
   cat("Warning: ",sum(nas)," curves with NA are omited\n")
   fdataobj$data<-fdataobj$data[bb,]
  y<-y[bb]
   }
else {
if (!any(nas) & any(nas.g)) {
   cat("Warning: ",sum(nas.g)," values of group with NA are omited \n")
   bb<-!nas.g
   fdataobj$data<-fdataobj$data[bb,]
     y<-y[bb]
   }
else {
if (any(nas) & any(nas.g))  {
   bb<-!nas & !nas.g
   cat("Warning: ",sum(!bb)," curves  and values of group with NA are omited \n")
   fdataobj$data<-fdataobj$data[bb,]
   y<-y[bb]
   }
}}

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
types=FALSE
mdist=metric(fdataobj,fdataobj,...)
ke<-deparse(substitute(Ker))
ty<-deparse(substitute(type.S))
if (is.null(h)) h=h.default(fdataobj,probs=c(0.025,0.25),len=25,metric = mdist,Ker =ke,
 type.S =ty,...)
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

