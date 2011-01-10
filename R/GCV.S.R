GCV.S=function(y,S,criteria="Rice",W=diag(1,ncol=ncol(S),nrow=nrow(S)),trim=0,draw=FALSE,...){
    fdat=TRUE
    if (is.fdata(y))    y2=t(y$data)
    else {
        if (is.matrix(y)&&(ncol(y)==1) ){y2<-y;fdat<-draw<-FALSE}
        else if (is.vector(y)){y2<-y;fdat<-draw<-FALSE}
        else stop("y is not a fdata,  vector or matrix")
    }
    n=ncol(S)
    y.est=S%*%y2
    tab=list("GCV","AIC","FPE","Shibata","Rice")
    type.i=pmatch(criteria,tab)
    e=y2-y.est
    if (trim>0) {
      if (fdat) {
         e = fdata(t(e),y$argvals,y$rangeval,y$names)
         ee<-norm.fdata(e)
         e.trunc=quantile(ee,probs=(1-trim),na.rm=TRUE,type=4)
         ind<-ee<=e.trunc
         if (draw)  plot(y,col=(2-ind))
         l<-which(abs(ee)<=e.trunc)
         res=traza((e[["data"]][l,])%*%W%*%t(e[["data"]][l,]))
        }
        else {
             ee = t(e)
             e.trunc=quantile(abs(ee),probs=(1-trim),na.rm=TRUE,type=4)
             l<-which(abs(ee)<=e.trunc)
             res=traza(t(e[l])%*%W[l,l]%*%e[l])   }
    }
    else      {l=1:n
               res=traza(t(e)%*%W%*%e)             }
    d<-diag(S)[l]
    if (is.na(type.i))   {
                   if (mean(d,na.rm=TRUE)>0.5) vv=Inf
                   else   vv= 1/(1-2*mean(d,na.rm=TRUE))        }
    else {
      if (type.i==1)  vv=(1-mean(d,na.rm=TRUE))^(-2)
      if (type.i==2)  vv= exp(2*mean(d))
      if (type.i==3)  vv=(1+mean(d,na.rm=TRUE))/(1-mean(d,na.rm=TRUE))
      if (type.i==4)  vv=(1+2*mean(d,na.rm=TRUE))
      if (type.i==5)  {
      if (mean(d,na.rm=TRUE)>0.5) vv=Inf
      else   vv= 1/(1-2*mean(d,na.rm=TRUE))   }          }
return(res*vv/n) }

