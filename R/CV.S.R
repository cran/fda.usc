CV.S=function (y, S, W = diag(ncol(S)),trim=0,draw=FALSE,...)
{
    n = ncol(S)
    fdat=TRUE
    if (is.fdata(y))    y2=t(y$data)
    else {
        if (is.matrix(y)&&(ncol(y)==1) ){y2<-y;fdat<-draw<-FALSE}
        else if (is.vector(y)){y2<-y;fdat<-draw<-FALSE}
        else stop("y is not a fdata,  vector or matrix")
    }
    y.est = S %*% y2
    I=diag(n)/(1-diag(S))^2
    W = W *I
    e<-y2 - y.est
    if (trim>0) {
        if (fdat) {
                e = fdata(t(e),y$argvals,y$rangeval,y$names)
                ee<-norm.fdata(e)
                e.trunc=quantile(ee,probs=(1-trim),na.rm=TRUE,type=4)
                ind<-ee<=e.trunc
#print("Trimmed curves")
#print(rownames(y[["data"]])[ind==F])
                if (draw)  plot(y,col=(2-ind))
                res = mean(diag(W)[ind]*e[ind,][["data"]]^2,na.rm=TRUE)    }
        else {
             ee = t(e)
             e.trunc=quantile(abs(ee),probs=(1-trim),na.rm=TRUE,type=4)
             l<-which(abs(ee)<=e.trunc)
             res = mean(diag(W)[l]*e[l]^2,na.rm=TRUE)            }
     }
     else        res = mean(diag(W)*e^2)
    if (is.nan(res)) res=Inf
return(res)
}

 
#CV.S_SinTruncar=function (y, S, W = diag(ncol(S)),...){
#    n = ncol(S)
#    if (is.matrix(y) ) {if (ncol(y)>1) y=t(y)}
#    y.est = S %*% y
#    e = y - y.est
#    I=diag(n)/(1-diag(S))^2
#    W = W *I
#    res = mean(diag(W)*e^2)
#    if (is.nan(res)) res=Inf
#    return(res)
#}



