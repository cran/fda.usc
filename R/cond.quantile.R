cond.quantile<-function(qua=0.5,fdata0,fdataobj,y,fn,
a=min(y),b=max(y),tol=10^floor(log10(max(y)-min(y))-3),iter.max=100,...){
  i<-0
  tol.up=qua+tol
  tol.lo=qua-tol
  medio<-(a+b)/2
  rmed=fn(fdata0,medio, fdataobj, y,...)$Fc
  rlo=fn(fdata0,a, fdataobj, y, ...)$Fc
  rup=fn(fdata0,b, fdataobj, y, ...)$Fc
  res<-c(i,medio,rmed)
  if ((rup*rlo>qua)|(a>b)) {return('Error in input data') }
  else {
    while (abs(rmed-qua)>tol & i<iter.max) {
      i<-i+1
      if (rmed<qua) {
        a<-medio
        medio<-(a+b)/2
        rlo=rmed
      }
      else {
        b<-medio
        medio<-(a+b)/2
        rup=rmed
      }
      res<-rbind(res,c(i,medio,rmed))
      rmed=fn(fdata0,medio, fdataobj, y,...)$Fc
  }
  }
print(paste("Fc=",medio,"with Conditional quantile=",qua,", tol=",tol," and iter=",i) )
return(medio)
}

