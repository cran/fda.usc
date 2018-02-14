fdata.deriv<-function(fdataobj,nderiv=1,method="bspline",class.out='fdata'
,nbasis=NULL,...) {
 if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
 nas1<-is.na.fdata(fdataobj)
 if (any(nas1))  stop("fdataobj contain ",sum(nas1)," curves with some NA value \n")
 DATA<-fdataobj[["data"]]
 tt=fdataobj[["argvals"]]
 rtt=fdataobj[["rangeval"]]
 labels=fdataobj[["names"]]
 nr <- nrow(DATA)
 nc <- ncol(DATA)
 if (method=="diff") {
  res=matrix(NA, nrow=nr, ncol = nc)
  for (i in 1:nr) {
   a=diff(DATA[i,],differences=nderiv)/(tt[2:nc]-tt[1:(nc-1)])^nderiv
   ab=matrix(NA,ncol=nc,nrow=2)
   ab[1,2:nc]=a
   ab[2,1:(nc-1)]=a
   res[i,]=colMeans(ab,na.rm=TRUE)
  }
  labels$main<-paste("d(",labels$main,",",nderiv,")",sep="")
  res<-fdata(res,tt,rtt,names=labels)
 }
  else {
      if (any(method==c("fmm", "periodic", "natural", "monoH.FC"))) {
       res=matrix(NA,nrow=nrow(DATA),ncol=ncol(DATA))
       for (i in 1:nrow(DATA)) {
         f1<-splinefun(x=tt,y=DATA[i,],method=method)
         res[i,]=f1(tt,deriv=nderiv)
        }
          labels$main<-paste("d(",labels$main,",",nderiv,")",sep="")
         res<-fdata(res,tt,rtt,names=labels)
       }
  else{
    if (any(method==c("bspline","exponential", "fourier",
      "monomial","polynomial"))) {
      #no run  "constant","polygonal","power"
      res=fdata2fd(fdataobj=fdataobj,type.basis=method,nbasis=nbasis,nderiv=nderiv,...)
      if (class.out=='fdata') {
         ff<-eval.fd(tt,res)
         labels$ylab<-paste("d(",labels$ylab,",",nderiv,")",sep="")
         res=fdata(t(ff),tt,rtt,names=labels)
         }
      }
  }
}
res
}