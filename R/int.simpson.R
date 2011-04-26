int.simpson=function(fdataobj,equi=TRUE,method="CSR"){
 if (!inherits(fdataobj, "fdata"))  fdataobj<-fdata(fdataobj)
 n<-nrow(fdataobj)
 out<-rep(NA,n)
 for (i in 1:n) {
   out[i]<-int.simpson2(fdataobj$argvals,fdataobj$data[i,],equi,method)
   }
	return(out)
}
int.simpson2=function(x,y,equi=TRUE,method="CSR"){
  n=length(x);ny=length(y)
	if (n!=ny) stop("Different length in the input data")
  out <- switch(method,
    "CSR" = {
     if (!equi){app=approx(x,y,n=2*length(x)-1);x=app$x;y=app$y}
     	h=(max(x)-min(x))/(n-1)
	   value=(h/3)*(y[n]+y[1]+2*sum(y[2*(1:((n-1)/2))+1])+4*sum(y[2*(1:((n-1)/2))]))
    },
    "ESR" = {
     if (!equi){app=approx(x,y,n=2*length(x)-1);x=app$x;y=app$y}
     h=(max(x)-min(x))/(n-1)
	   if (n<=4) stop("This method needs n>4")
     value=17*(y[1]+y[n])+59*(y[2]+y[n-1])+43*(y[3]+y[n-2])+49*(y[4]+y[n-3])
	   value=value+48*sum(y[5:(n-4)])
  	 value=(h/48)*value
    },
    "TRAPZ" = {
   		idx=2:length(x)
	    value<-as.double((x[idx]-x[idx-1])%*%(y[idx]+y[idx-1]))/2
	  }
   )
	return(out)
}

