cond.mode<-function(Fc,method="monoH.FC") {
   x=Fc$y0
   Fc=Fc$Fc
  ndist=length(x)
  par(mfrow=c(2,1))
  if (method=="diff") {
     plot(x,Fc,type="l",main="Conditional Distribution Function",ylab="Fc")
     if (is.vector(Fc)) Fc=matrix(Fc,ncol=1)
     ind3=1
	   der = (Fc[2:ndist,ind3] - Fc[1:(ndist-1),ind3]) /( x[2:ndist]-x[1:(ndist-1)])
     names(der) = round(x[2:(ndist)],2)
   }
  else {
    f<-splinefun(x,Fc,method=method)
    curve(f(x),x[1],x[length(Fc)],main="Conditional Distribution Function",ylab="Fc")
    der=f(x,deriv=1)
  }
  max.der=max(der)
  ind.mode.cond=which.max(der)
  mode.cond=x[ind.mode.cond]
 plot(der,type="l",main=c("Conditional mode=",round(mode.cond,3)),
 xlab=c("Index conditional mode",ind.mode.cond),ylab="Derivative")
 abline(v=ind.mode.cond,col=2,lty=2)
  return(list("mode.cond"=mode.cond,"ind.mode.cond"=ind.mode.cond))
}

