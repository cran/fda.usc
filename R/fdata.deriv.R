fdata.deriv<-function(fdataobj,nderiv=1,method="bspline",class.out='fdata'
,nbasis=NULL,...) {
 if (!is.fdata(fdataobj))  fdataobj=fdata(fdataobj)
 DATA<-fdataobj[["data"]]
 tt=fdataobj[["argvals"]]
 labels=fdataobj[["labels"]]
 ndist=ncol(DATA)
 if (method=="diff") {
  res=matrix(NA,nrow=nrow(DATA),ncol=ncol(DATA))
  for (i in 1:nrow(DATA)) {
   a=diff(DATA[i,],differences=nderiv)/(tt[2:ndist]-tt[1:(ndist-1)])
   ab=matrix(NA,ncol=ndist,nrow=2)
   ab[1,2:ndist]=a
   ab[2,1:(ndist-1)]=a
   res[i,]=apply(ab,2,mean,na.rm=TRUE)
  }
  res<-fdata(res,tt)
 }
  else {
      if (any(method==c("fmm", "periodic", "natural", "monoH.FC"))) {
       print(method)
       res=matrix(NA,nrow=nrow(DATA),ncol=ncol(DATA))
       for (i in 1:nrow(DATA)) {
         f1<-splinefun(x=tt,y=DATA[i,],method=method)
         res[i,]=f1(tt,deriv=nderiv)
        }
         res<-fdata(res,tt)
       }
      else{
      if (any(method==c("bspline","exponential", "fourier",
      "monomial","polynomial"))) {
#no run  "constant","polygonal","power"
      res<-fdata(DATA,tt)
      res=fdata2fd(res,method,nbasis,nderiv,...)
      if (class.out=='fdata') res=fdata(res,tt,labels)
      }
  }
  }
res
}
